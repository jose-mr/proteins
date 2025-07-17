# python imports
from collections import defaultdict
from urllib.request import urlretrieve
import gzip
import tarfile

# django imports
from django.db import models
from django.contrib.postgres.fields import ArrayField

# library imports
from pseudoenzymes.settings import (NCBI_TAXDUMP_FILE, NCBI_FOLDER, NCBI_NAMES_FILE,
                                    NCBI_NODES_FILE, NCBI_MERGED_FILE)

class TaxonQuerySet(models.QuerySet):

    def old_to_new(self) -> dict[int, int]:
        out = {}
        for taxid, old_taxids in self.values_list("taxid", "old_taxids"):
            out[taxid] = taxid
            for old_taxid in  old_taxids:
                out[old_taxid] = taxid
        return out

    def children_of(self, parents):
        """recursively find all taxa that are children of these parent taxa

        parents are included in the resulting query
        """
        ids = set(parents.values_list("taxid", flat=True))
        no_previous_ids = 0
        while no_previous_ids != len(ids):
            no_previous_ids = len(ids)
            ids.update(Taxon.objects.filter(parent__in=ids).values_list("taxid", flat=True))
        return self.filter(pk__in=ids)

class Taxon(models.Model):
    """A biological species"""
    taxid = models.IntegerField(verbose_name="NCBI TaxID", primary_key=True)
    parent = models.ForeignKey(
            "Taxon",
            null=True,
            on_delete=models.SET_NULL,
            related_name="children"
    )
    old_taxids = ArrayField(models.IntegerField())
    rank = models.CharField(max_length=30)
    name = models.CharField(max_length=255, blank=False)
    common_name = models.CharField(max_length=255, blank=True)

    objects = TaxonQuerySet.as_manager()

    def __str__(self):
        output = f"txid{self.taxid}: {self.name}"
        if self.common_name:
            output += " ({})".format(self.common_name)
        return output

    class Meta:
        verbose_name_plural = "Taxa"

    @classmethod
    def download_taxdump_from_ncbi(cls):
        """Downloads and extracts the ncbi taxonomy files"""
        urlretrieve("https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", NCBI_TAXDUMP_FILE)
        with gzip.open(NCBI_TAXDUMP_FILE, 'rb') as gz_file:
            with tarfile.open(fileobj=gz_file) as tar:
                tar.extractall(path=NCBI_FOLDER)

    @classmethod
    def create_from_ncbi_files(cls):
        """Creates new taxa objects"""
        existing = set(cls.objects.values_list("taxid", flat=True))

        info = defaultdict(dict)
        with open(NCBI_NAMES_FILE, "rt") as names_file:
            for line in names_file:
                words = [word.strip() for word in line.split("|")]
                taxid, name, _, name_class, _ = words
                taxid = int(taxid)
                if name_class == "scientific name":
                    info[taxid]["name"] = name
                    info[taxid]["old_taxids"] = []
                elif name_class == "genbank common name":
                    info[taxid]["common_name"] = name

        with open(NCBI_NODES_FILE, "rt") as nodes_file:
            for line in nodes_file:
                words = [word.strip() for word in line.split("|")[:3]]
                taxid, parent_taxid, rank = words
                taxid = int(taxid)
                info[taxid]["rank"] = rank.lower()
                info[taxid]["parent_id"] = parent_taxid

        with open(NCBI_MERGED_FILE, "rt") as merged_file:
            for line in merged_file:
                old, new = [int(word.strip()) for word in line.split("|")[:2]]
                info[new]["old_taxids"].append(old)

        to_create = [cls(taxid=taxid,**v) for taxid, v in info.items()
                     if taxid not in existing]

        print(f"Creating {len(to_create)} new Taxon objects")
        cls.objects.bulk_create(to_create, batch_size=10000)
