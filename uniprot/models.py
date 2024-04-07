import gzip

from Bio import SeqIO

from django.db import models
from pseudoenzymes.settings import DATA_FOLDER
from django.db import connection

DAT_FILE = DATA_FOLDER / "uniprot/uniprot_sprot.dat.gz"
CURRENT_RELEASE_FILE = DATA_FOLDER / "uniprot/reldate.txt"
INDEX_FILE = DATA_FOLDER / "uniprot/index"
INDEX_RELEASE_FILE = DATA_FOLDER / "uniprot/index_release"

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
    seq = models.TextField()


    @classmethod
    def create_from_dat_file(cls):
        """Create all UniProt entries from a uniprot dat_gz file"""
        batch_size = 100000
        with gzip.open(DAT_FILE, "rb") as dat_file:
            records = []
            for record in SeqIO.parse(dat_file, "swiss"):
                records.append(record)
                if len(records) == batch_size:
                    cls.create_from_records(records)
                    records = []
            cls.create_from_records(records)

    @classmethod
    def create_from_records(cls, records):
        "create UniProt objects biopython record objects"
        print(f"Creating entries from {len(records)} records")
        to_create = []
        kw_to_create = set()
        entry_keyword_through_objs = []
        for record in records:
            # prepare keyword through objects
            keywords = record.annotations.get("keywords", [])
            keywords = [kw.lower() for kw in keywords]
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
                seq=record.seq,
                taxid=int(record.annotations['ncbi_taxid'][0]),
                comment=record.annotations.get("comment", ""),
                secondary_ac=secondary_ac,
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


class Keyword(models.Model):
    name = models.TextField(primary_key=True)
