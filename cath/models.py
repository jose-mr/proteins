"""Includes models related with the CATH database"""

# python standard imports
import gzip
from urllib.request import urlretrieve

# django imports
from django.db import models
from django.db.models import Count, Q

# pseudoenzymes imports
import uniprot.models as uniprot
from pseudoenzymes.settings import CATH_NAMES_FILE, INTERPRO_ONLY_G3_SP_DAT_FILE

class SuperfamilyQuerySet(models.QuerySet):

    def superfamilies(self):
        superfamilies = [obj.number for obj in self if
                         obj.number.count(".") == 3]
        return self.filter(number__in=superfamilies)

    def annotate_uniprot_count(self):
        return self.annotate(uniprot_entries_count=Count('uniprot_entries'))

    def annotate_uniprot_enzyme_ec_count(self):
        with_ec = Count("uniprot_entries", filter=Q(uniprot_entries__ec_entries__isnull=False))
        # with_ec = Count("uniprot_entries", filter=Q(
            # uniprot_entries__in=uniprot.Entry.objects.enzymes_ec()))
        no_ec = Count("uniprot_entries", filter=Q(uniprot_entries__ec_entries__isnull=True))
        return self.annotate(uniprot_entries_ec_count=with_ec,
                             uniprot_entries_no_ec_count=no_ec)


class Superfamily(models.Model):
    """Model for the CATH superfamilies"""
    number = models.TextField(primary_key=True, verbose_name="CATH number")
    name = models.TextField(verbose_name="CATH name")
    uniprot_entries = models.ManyToManyField(
            "uniprot.Entry",
            related_name="cath_superfamilies",
            through="SuperfamilyUniprotEntry"
            )

    objects = SuperfamilyQuerySet.as_manager()

    def __repr__(self):
        return str(self.number)

    def __str__(self):
        return str(self.number)

    @classmethod
    def download_names_file(cls):
        """download the cath-names.txt file from the CATH website"""
        urlretrieve(
            "http://download.cathdb.info/cath/releases/daily-release/"
            "newest/cath-b-newest-names.gz",
            CATH_NAMES_FILE
        )

    @classmethod
    def create_from_names_file(cls):
        """create CATH objects by reading the cath-names.txt file"""
        objs = []
        with gzip.open(CATH_NAMES_FILE, 'rt', encoding='utf8') as names_file:
            for line in names_file:
                number, name = line.split(' ', 1)
                objs.append(cls(number=number, name=name.strip()))
        print(f"Creating {len(objs)} CATH superfamilies and their parents")
        cls.objects.bulk_create(objs)


class SuperfamilyUniprotEntry(models.Model):
    """through table that links uniprot entries to CATH superfamilies"""
    uniprot_entry = models.ForeignKey('uniprot.Entry', on_delete=models.CASCADE)
    superfamily = models.ForeignKey("Superfamily", on_delete=models.CASCADE)
    start = models.IntegerField()
    end = models.IntegerField()
    # evalue = models.FloatField()

    class Meta:
        unique_together = ("uniprot_entry", "superfamily", "start", "end")

    @classmethod
    def create_from_interpro_file(cls):
        """Update table using Gene3D file"""
        objs = []
        uniprot_entries = set(uniprot.Entry.objects.values_list("ac", flat=True))
        with open(INTERPRO_ONLY_G3_SP_DAT_FILE, 'r') as interpro_file:
            for line in interpro_file:
                ac, _, _, g3d, start, end = line.strip().split("\t")
                superfamily = g3d.split(":")[1]
                # file contain some uniprot entries that do not belong to swissprot
                if ac in uniprot_entries:
                    objs.append(cls(
                        uniprot_entry_id=ac,
                        superfamily_id=superfamily,
                        start=start,
                        end=end,
                        ))
        print(f"Creating {len(objs)} cath <-> uniprot associations")
        cls.objects.bulk_create(objs, batch_size=100000)
