import gzip
from urllib.request import urlretrieve
import urllib.request
import hashlib

from Bio import SeqIO, SwissProt

from django.db.models import Q, Count, Func
from django.db.models.functions import Length
from django.db import models, transaction
from pseudoenzymes.settings import SWISSPROT_DAT_FILE, SWISSPROT_ACS_FILE, PDB_UNIPROT_DAT_FILE
import go.models as go
import taxonomy.models as taxonomy

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

    objects = EntryQuerySet.as_manager()

    def __str__(self):
        return self.ac
    
    @classmethod
    def create_from_dat_file(cls, filename):
        """Create all UniProt entries from a uniprot dat_gz file"""
        batch_size = 100000
        with gzip.open(filename, "rb") as dat_file:
            records = []
            for record in SwissProt.parse(dat_file):
                records.append(record)
                if len(records) == batch_size:
                    cls.create_from_records(records)
                    records = []
            cls.create_from_records(records)

    @classmethod
    def create_from_ac_list(cls, acs):
        """create entries from a list of accessions.
        This is very slow use only for a small number of entries"""
        records = []
        batch_size = 100
        for i, ac in enumerate(acs):
            url = f"https://rest.uniprot.org/uniprotkb/{ac}.txt"
            with urllib.request.urlopen(url) as result:
                biorecord = SwissProt.parse(result)
                records.extend([r for r in biorecord])
            if len(records) == batch_size:
                cls.create_from_records(records)
                records = []
                print(".", end="")
        cls.create_from_records(records)

    @classmethod
    def create_from_records(cls, records):
        "create UniProt objects biopython record objects"
        print(f"Creating entries from {len(records)} records")
   
        # Add new Sequence objects
        existing_seqs = set(Sequence.objects.values_list("seq_hash", flat=True))
        to_create = []
        for record in records:
            if (seq_hash := get_seq_hash(record.sequence)) not in existing_seqs:
                to_create.append(Sequence(seq=record.sequence, seq_hash=seq_hash))
                existing_seqs.add(seq_hash)
        print(f"Creating {len(to_create)} new sequence objects")
        Sequence.objects.bulk_create(to_create)
        seq2id = {seq.seq: seq.id for seq in Sequence.objects.all()}
        # this is needed to include alternative tax ids
        old_to_new_taxid = taxonomy.Taxon.objects.all().old_to_new()

        to_create = []
        to_update = []
        kws_records = set()
        entry_keywords = []
        existing_entry_kws = set(cls.keywords.through.objects.all().values_list("entry_id", "keyword_id"))
        existing_entries = set(cls.objects.all().values_list("ac", flat=True))
        for record in records:
            ac = record.accessions[0]
            keywords = clean_kws_from_record(record)
            kws_records.update(keywords)

            # prepare keyword through objects
            for kw in keywords:
                if (record.accessions[0], kw) not in existing_entry_kws:
                    entry_keywords.append(cls.keywords.through(entry_id=ac, keyword_id=kw))

            obj = cls(
                ac=ac,
                name=record.entry_name,
                seq_id=seq2id[record.sequence],
                species=old_to_new_taxid.get(int(record.taxonomy_id[0]), None),
                secondary_ac=record.accessions[1:],
                comment=record.comments,
                reviewed=record.data_class=="Reviewed",
            )
            if obj.ac not in existing_entries:
                to_create.append(obj)
            else:
                to_update.append(obj)

        print(f"Creating {len(to_create)} new UniProt entries")
        created = cls.objects.bulk_create(to_create)

        print(f"Updating {len(to_update)} UniProt Entries")
        # TODO this is to slow, need to compare first and only update what is necessary
        # cls.objects.bulk_update(to_update, fields=["comment", "reviewed", "secondary_ac"])

        # Add new keywords
        existing_kws = set(Keyword.objects.all().values_list("name", flat=True))
        kw_to_create = kws_records - existing_kws
        print(f"Creating {len(kw_to_create)} new Keywords")
        Keyword.objects.bulk_create([Keyword(name=kw) for kw in kw_to_create])

        print(f"Creating {len(entry_keywords)} new UniProt - Keywords relations")
        cls.keywords.through.objects.bulk_create(entry_keywords)
        return created

    @classmethod
    def dump_all_acs(cls):
        "write all primary acs and secondary acs into a file"
        acs = cls.objects.values_list("ac", flat=True)
        with open(SWISSPROT_ACS_FILE, "w") as acs_file:
            acs_file.writelines(acs)

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

class Sequence(models.Model):
    seq = models.TextField()
    seq_hash = models.CharField(max_length=32, unique=True, editable=False)

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


def get_seq_hash(seq):
    return hashlib.md5(seq.encode("utf-8")).hexdigest()
