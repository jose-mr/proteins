"""Includes models related with the CATH database"""

# python standard imports
import gzip
from urllib.request import urlretrieve

# django imports
from django.db import models

# pseudoenzymes imports
from pseudoenzymes.settings import CATH_NAMES_FILE, INTERPRO_ONLY_G3_SP_DAT_FILE


class Superfamily(models.Model):
    """Model for the CATH superfamilies"""
    number = models.TextField(primary_key=True)
    name = models.TextField()

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
    evalue = models.FloatField()

    class Meta:
        unique_together = ("uniprot_entry", "superfamily", "start", "end")

    @classmethod
    def create_from_interpro_file(cls):
        """Update table using Gene3D file"""
        objs = []
        with open(INTERPRO_ONLY_G3_SP_DAT_FILE, 'r') as interpro_file:
            for line in interpro_file:
                print(line)
            # csv_reader = csv.reader(cath_file)
            # for no, words in enumerate(csv_reader):
                # uniprot, _, _, _, _, cath, all_limits, evalue, _ = words
                # for limits in all_limits.split(","):
                    # start, end = limits.split("-")
                    # try:
                        # info_gene3d.add(
                            # (uniprot_to_id_dict[uniprot], cath_to_id_dict[cath], int(start), int(end), float(evalue)))
                    # except KeyError:
                        # continue
        # existing_pairs = set(cls.objects.values_list('sequence_id', 'cath_id', 'start', 'end', 'evalue'))
        # to_delete = existing_pairs - info_gene3d
        # print("Deleting {} cath<->Uniprot links".format(len(to_delete)))
        # for t in to_delete:
            # cls.objects.get(sequence_id=t[0], cath_id=t[1], start=t[2], end=t[3]).delete()
        # to_create = info_gene3d - existing_pairs
        # print("Adding {} new cath<->UniProt links".format(len(to_create)))
        # cls.objects.bulk_create([cls(sequence_id=t[0], cath_id=t[1], start=t[2], end=t[3], evalue=t[4])
                                 # for t in to_create], batch_size=5000)

