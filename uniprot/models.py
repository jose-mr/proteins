import gzip
from urllib.request import urlretrieve
import urllib.request

from Bio import SeqIO

from django.db.models import Q, Count
from django.db import models, transaction
from pseudoenzymes.settings import SWISSPROT_DAT_FILE, SWISSPROT_ACS_FILE, PDB_UNIPROT_DAT_FILE
import go.models as go

class EntryQuerySet(models.QuerySet):

    def reviewed(self):
        return self.filter(reviewed=True)

    def enzymes_ec(self, catalytic=True):
        return self.filter(ec_entries__isnull=not(catalytic))

    def enzymes_kw(self, catalytic=True):
        enzymes = self.filter(keywords__in=Keyword.objects.enzymatic())
        if catalytic:
            return enzymes
        else:
            return self.exclude(ac__in=enzymes)

    def enzymes_go(self, catalytic=True):
        enzymes = self.filter(go_associations__in=go.TermUniProtEntry.objects.catalytic())

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
    taxid = models.IntegerField(
            verbose_name="NCBI TaxId"
    )
    keywords = models.ManyToManyField("Keyword", related_name="entries")
    seq = models.ForeignKey("sequence", related_name="uniprot_entries", on_delete=models.PROTECT)
    reviewed = models.BooleanField(default=True)

    objects = EntryQuerySet.as_manager()

    def __str__(self):
        return self.ac
    
    @classmethod
    def download_pdb_dat_file(cls):
        """get the dat file for the accession codes in the list
        it is used to get entries associated with PDBs that are not in uniprot"""

        print("downloading uniprot data for entries associated with pdb")
        url = ("https://rest.uniprot.org/uniprotkb/search?query=((structure_3d:true)"
               "+AND+(reviewed:false))&compressed=true&format=txt&size=500")
        
        with open(PDB_UNIPROT_DAT_FILE, "wb") as out_file:
            while True:
                with urllib.request.urlopen(url) as result:
                    out_file.write(result.read())
                    link = result.getheader("Link")
                    if link is None:
                        break
                    url = link.split(";")[0][1:-1]

    @classmethod
    def create_from_dat_file(cls, filename=SWISSPROT_DAT_FILE, swissprot=True):
        """Create all UniProt entries from a uniprot dat_gz file"""
        batch_size = 100000
        with gzip.open(filename, "rb") as dat_file:
            records = []
            for record in SeqIO.parse(dat_file, "swiss"):
                records.append(record)
                if len(records) == batch_size:
                    cls.create_from_records(records, swissprot=swissprot)
                    records = []
            cls.create_from_records(records, swissprot=swissprot)

    @classmethod
    def create_from_ac_list(cls, acs, swissprot=False):
        """create entries from a list of accessions.
        This is very slow use only for a small number of entries"""
        records = []
        for i, ac in enumerate(acs):
            url = f"https://rest.uniprot.org/uniprotkb/{ac}.txt"
            with urllib.request.urlopen(url) as result:
                biorecord = SeqIO.parse(result, format="swiss")
                records.extend([r for r in biorecord])
        cls.create_from_records(records, swissprot=swissprot)

    @classmethod
    def create_from_pdb_dat_file(cls):
        return cls.create_from_dat_file(filename=PDB_UNIPROT_DAT_FILE, swissprot=False)

    @classmethod
    def create_from_records(cls, records, swissprot=True):
        "create UniProt objects biopython record objects"
        print(f"Creating entries from {len(records)} records")
        to_create = []
        kw_to_create = set()
        entry_keyword_through_objs = []
        # sequence objects
        seqs_to_create = []
        existing_seqs = set(Sequence.objects.values_list("seq", flat=True))
        for record in records:
            if record.seq not in existing_seqs:
                seqs_to_create.append(Sequence(seq=record.seq))
                existing_seqs.add(record.seq)
        Sequence.objects.bulk_create(seqs_to_create) 
        seq2id = {seq.seq: seq.id for seq in Sequence.objects.all()}

        for record in records:
            # prepare keyword through objects
            keywords = record.annotations.get("keywords", [])
            keywords = [kw.lower() for kw in keywords]

            # this is necessary because keywords sometimes are not cleaned up properly
            # this happens when the dat file comes from the api
            keywords = [k[:k.find("{eco")+1].strip("{").strip() for k in keywords
                        if not (k.startswith("eco:") or k.startswith("rule:") or k.startswith("prorule:"))]
            keywords = [k for k in keywords if k]
            # above should not be necessary if biopython fixes this
            print(keywords)

            entry_keyword_through_objs.extend([
                cls.keywords.through(entry_id=record.id, keyword_id=kw)
                for kw in keywords
            ])
            kw_to_create.update(keywords)

            secondary_ac = record.annotations["accessions"]
            secondary_ac.remove(record.id)

            to_create.append(cls(
                ac=record.id,
                name=record.description.split("=")[1].split(";")[0].split(" {")[0],
                seq_id=seq2id[record.seq],
                taxid=int(record.annotations['ncbi_taxid'][0]),
                comment=record.annotations.get("comment", ""),
                secondary_ac=secondary_ac,
                reviewed=swissprot,
            ))
        created = cls.objects.bulk_create(to_create)
        print(f"Created {len(created)} UniProt entries")

        kw_existing = set(Keyword.objects.all().values_list("name", flat=True))
        kw_to_create = kw_to_create - kw_existing
        kw_to_create = [Keyword(kw) for kw in kw_to_create]
        kw_created = Keyword.objects.bulk_create(kw_to_create)
        print(f"Created {len(kw_created)} Keywords")

        rel_created = cls.keywords.through.objects.bulk_create(entry_keyword_through_objs)
        print(f"Created {len(rel_created)} UniProt - Keywords relations")
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

class Sequence(models.Model):
    seq = models.TextField(unique=True)

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
