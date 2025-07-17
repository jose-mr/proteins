import gzip
from urllib.request import urlretrieve
from itertools import islice
import urllib.request
import hashlib

from Bio import SeqIO, SwissProt
from Bio.SeqFeature import UnknownPosition

from django.db import models, transaction
from django.db.models import Q, Count, Func
from django.db.models.functions import Length
import go.models as go
import taxonomy.models as taxonomy

class Sequence(models.Model):
    seq = models.TextField()
    seq_hash = models.CharField(max_length=32, unique=True, editable=False)

def get_seq_hash(seq):
    return hashlib.md5(seq.encode("utf-8")).hexdigest()

class EntryQuerySet(models.QuerySet):

    def reviewed(self):
        return self.filter(reviewed=True)

    def peptides(self, max_length=50):
        qs = self.annotate(seq_length=Length("seq__seq"))
        return qs.filter(seq_length__lte=max_length)

    def enzymes_ec(self, catalytic=True):
        return self.filter(ec_entries__isnull=not(catalytic))

    def enzymes_kw(self, catalytic=True):
        enzymes = self.filter(keywords__in=Keyword.objects.enzymatic())
        if catalytic:
            return enzymes
        else:
            return self.exclude(ac__in=enzymes)

    def enzymes_go(self, catalytic=True):
        # enzymes = self.filter(go_associations__in=go.TermUniProtEntry.objects.catalytic())

        catalytic_association = go.TermUniProtEntry.objects.filter(
                term__in=go.Term.objects.catalytic(),
                qualifier="enables"
                )
        enzymes = self.filter(go_associations__in=catalytic_association)
        if catalytic:
            return enzymes
        else:
            return self.exclude(ac__in=enzymes)

    def catalytic_activity(self):
        return self.filter(comment__contains="CATALYTIC ACTIVITY:")

    def caution(self):
        return self.filter(comment__contains="CAUTION:")

    def inactive(self):
        return self.filter(
                Q(name__istartswith="inactive")|
                Q(name__istartswith="probable inactive")|
                Q(name__istartswith="probably inactive")|
                Q(name__istartswith="putative inactive")
                )

    def single_domain(self):
        """return uniprot entries that are associated with a single CATH superfamily"""
        return self.annotate(domain_count=Count("cath_superfamilies")).filter(domain_count=1)\
                .prefetch_related("cath_superfamilies")

class Entry(models.Model):
    ac = models.CharField(
            primary_key=True,
            max_length=10,
            verbose_name="primary accession number",
    )
    secondary_ac = models.JSONField(
            verbose_name="secondary accession numbers",
            default=list
    )
    name = models.TextField()
    comment = models.TextField()
    keywords = models.ManyToManyField("Keyword", related_name="entries")
    reviewed = models.BooleanField(default=True)
    seq = models.ForeignKey("sequence", related_name="uniprot_entries", on_delete=models.PROTECT)
    species = models.ForeignKey(
            "taxonomy.Taxon",
            related_name="uniprot_entries",
            null=True,
            on_delete=models.SET_NULL
    )
    features = models.JSONField()

    objects = EntryQuerySet.as_manager()

    def __str__(self):
        return self.ac
    
    @classmethod
    def create_from_dat_file(cls, filename, skip_first=0):
        """Create all UniProt entries from a uniprot dat_gz file does not update old entries"""
        batch_size = 100000

        taxid_map = taxonomy.Taxon.objects.all().old_to_new()
        with gzip.open(filename, "rb") as dat_file:
            records = []
            for no, record in enumerate(islice(SwissProt.parse(dat_file), skip_first, None)):
                if no > 0 and (no % batch_size) == 0:
                    cls.create_from_records(records, taxid_map)
                    records = []
                records.append(record)
            cls.create_from_records(records, taxid_map)

    @classmethod
    @transaction.atomic
    def create_from_records(cls, records,taxid_map):
        "create UniProt objects biopython record objects"
        print(f"Creating entries from {len(records)} records")

        # checking for existing objects is done inside the function to save memory
        # skipping existing entries
        ac2rec = {rec.accessions[0]: rec for rec in records}
        existing_acs = set(cls.objects.filter(ac__in=ac2rec).values_list("ac", flat=True))
        records = [rec for ac, rec in ac2rec.items() if ac not in existing_acs]
        print(f"Records not in db: {len(records)}")
   
        # adding missing sequence objects
        hex2rec = {get_seq_hash(rec.sequence): rec for rec in records}
        hex2id = {t[0]: t[1] for t in Sequence.objects.filter(seq_hash__in=hex2rec)\
            .values_list("seq_hash","id")}
        seq_created = Sequence.objects.bulk_create(
                [Sequence(seq=rec.sequence, seq_hash=seq_hash) for seq_hash, rec in
                 hex2rec.items() if seq_hash not in hex2id],
                batch_size=10000
        )
        hex2id.update({ seq.seq_hash: seq.pk for seq in seq_created})
        print(f"{len(seq_created)} new sequence objects created")

        # now adding the uniprot entries
        to_create = []
        kw_recs = set()
        entry_keywords = []
        print("looping over records")
        for record in records:
            ac = record.accessions[0]
            keywords = clean_kws_from_record(record)
            kw_recs.update(keywords)

            # prepare keyword through objects
            for kw in keywords:
                entry_keywords.append(cls.keywords.through(entry_id=ac, keyword_id=kw))
            
            features = [
                    {
                        "location": [
                            f.location.start if not isinstance(f.location.start, UnknownPosition) else None, 
                            f.location.end if not isinstance(f.location.end, UnknownPosition) else None, 
                            ],
                        "type": f.type,
                        "qualifiers": f.qualifiers
                    }
                    for f in record.features]

            obj = cls(
                ac=ac,
                name=record.entry_name,
                seq_id=hex2id[get_seq_hash(record.sequence)],
                species_id=taxid_map.get(int(record.taxonomy_id[0]), None),
                secondary_ac=record.accessions[1:],
                comment=record.comments,
                features=features,
                reviewed=record.data_class=="Reviewed",
            )
            to_create.append(obj)

        print(f"Creating {len(to_create)} new UniProt entries")
        cls.objects.bulk_create(to_create, batch_size=5000)

        # Add new keywords
        db_kws = set(Keyword.objects.filter(name__in=kw_recs).values_list("name", flat=True))
        Keyword.objects.bulk_create([Keyword(name=kw) for kw in kw_recs if kw not in db_kws])

        print(f"Creating {len(entry_keywords)} new UniProt - Keywords relations")
        cls.keywords.through.objects.bulk_create(entry_keywords, batch_size=10000)
        print(f"Done creating {len(to_create)} records and related annotations")


    @property
    def catalytic_activities(self):
        return [line for line in self.comment.split("\n") if line.startswith("CATALYTIC ACTIVITY:")]

    @property
    def cautions(self):
        return [line for line in self.comment.split("\n") if line.startswith("CAUTION:")]

class KeywordQuerySet(models.QuerySet):

    def enzymatic(self):
        # see https://www.uniprot.org/keywords/KW-9992
        return self.filter(name__in=Keyword.TOP_LEVEL_EZ_KWS)

class Keyword(models.Model):
    name = models.TextField(primary_key=True)

    TOP_LEVEL_EZ_KWS = {
          "oxidoreductase",
          "transferase",
          "hydrolase",
          "lyase",
          "isomerase",
          "ligase",
          "translocase",
          "allosteric enzyme",
          "dna invertase",
          "excision nuclease",
        }

    objects = KeywordQuerySet.as_manager()

class MD5(Func):
    function = "MD5"
    template = "%(function)s(%(expressions)s)"

class MultiSequenceAlignment(models.Model):
    name = models.TextField()
    family = models.TextField()
    seq = models.ForeignKey(
            "Sequence",
            related_name="multi_sequence_alignments",
            on_delete=models.CASCADE
            )
    alignment = models.TextField()

    class Meta:
        unique_together = ['name', 'family', 'seq']


def clean_kws_from_record(record):
    """get clean keywords from sequence record object
    
    some cleanup is necessary for the txt records coming from the rest api
    biopython does not divide the keywords correctly if they contain braces
    """
    clean_kws = []
    for kw in record.keywords:
        kw = kw.lower()
        if kw.startswith("eco:") or kw.startswith("prorule:") or kw.startswith("rule:"):
            continue
        else:
            if "eco:" in kw:
                kw = kw[:kw.find(" {eco:")]
            if kw and not kw.endswith("}"):
                clean_kws.append(kw)
    return clean_kws


