"""Includes models related with the CATH database"""

# python standard imports
import gzip
from urllib.request import urlretrieve

# django imports
from django.db import models
from django.db.models import Count, Q, Prefetch

# pseudoenzymes imports
import uniprot.models as uniprot
from pseudoenzymes.settings import CATH_NAMES_FILE, INTERPRO_ONLY_G3_SP_DAT_FILE

class SuperfamilyQuerySet(models.QuerySet):

    def superfamilies(self) -> models.QuerySet:
        superfamilies = [obj.number for obj in self if
                         obj.number.count(".") == 3]
        return self.filter(number__in=superfamilies)

    def annotate_uniprot_count(self) -> models.QuerySet:
        return self.annotate(uniprot_entries_count=Count('uniprot_entries'))

    def annotate_uniprot_enzyme_ec_count(self) -> models.QuerySet:
        with_ec = Count("uniprot_entries", filter=Q(uniprot_entries__ec_entries__isnull=False))
        # with_ec = Count("uniprot_entries", filter=Q(
            # uniprot_entries__in=uniprot.Entry.objects.enzymes_ec()))
        no_ec = Count("uniprot_entries", filter=Q(uniprot_entries__ec_entries__isnull=True))
        return self.annotate(uniprot_entries_ec_count=with_ec,
                             uniprot_entries_no_ec_count=no_ec)

    def annotate_single_domain_sequences(self):
        """return the with a list of sequences that are comprised by this single domain"""
        return self.filter(uniprot_entries__in=uniprot.Entry.objects.single_domain())\
            .prefetch_related(
                Prefetch(
                    "uniprot_entries",
                    queryset=uniprot.Entry.objects.single_domain(),
                    to_attr="single_domain_sequences"
                    )
                )


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

class SuperfamilyUniprotEntryQuerySet(models.QuerySet):

    def single_domain_sequences(self):
        single_domain_sequences = SuperfamilyUniprotEntry.objects\
                .values_list("uniprot_entry_id", flat=True)\
                .annotate(Count('id'))\
                .order_by()\
                .filter(id__count=1)
        return self.filter(uniprot_entry_id__in=set(single_domain_sequences))

class SuperfamilyUniprotEntry(models.Model):
    """through table that links uniprot entries to CATH superfamilies"""
    uniprot_entry = models.ForeignKey('uniprot.Entry', on_delete=models.CASCADE)
    superfamily = models.ForeignKey("Superfamily", on_delete=models.CASCADE)
    start = models.IntegerField()
    end = models.IntegerField()
    seq = models.ForeignKey('uniprot.Sequence', on_delete=models.PROTECT)
    # evalue = models.FloatField()

    objects = SuperfamilyUniprotEntryQuerySet.as_manager()

    class Meta:
        unique_together = ("uniprot_entry", "superfamily", "start", "end")

    @classmethod
    def create_from_interpro_file(cls):
        """Update table using Gene3D file"""
        objs = []
        uniprot_entries_to_seq = dict(uniprot.Entry.objects.values_list("ac", "seq__seq"))
        existing_seqs = set(uniprot.Sequence.objects.values_list("seq", flat=True))
        ac_sup_start_end_to_seq = {}
        seq_to_create = []
        with open(INTERPRO_ONLY_G3_SP_DAT_FILE, 'r') as interpro_file:
            for line in interpro_file:
                ac, _, _, g3d, start, end = line.strip().split("\t")
                superfamily = g3d.split(":")[1]
                # file contain some uniprot entries that do not belong to swissprot
                if ac in uniprot_entries_to_seq:
                    seq = uniprot_entries_to_seq[ac][int(start)-1:int(end)]
                    if seq not in existing_seqs:
                        existing_seqs.add(seq)
                        seq_to_create.append(uniprot.Sequence(seq=seq))

                    ac_sup_start_end_to_seq[(ac, superfamily, start, end)] = seq
                    objs.append(cls(
                        uniprot_entry_id=ac,
                        superfamily_id=superfamily,
                        start=start,
                        end=end,
                        ))
        uniprot.Sequence.objects.bulk_create(seq_to_create)
        seq_to_seq_obj = {seq.seq: seq for seq in uniprot.Sequence.objects.all()}
        for obj in objs:
            obj.seq = seq_to_seq_obj[
                    ac_sup_start_end_to_seq[(obj.uniprot_entry_id, obj.superfamily_id, obj.start, obj.end)]]

        print(f"Creating {len(objs)} cath <-> uniprot associations")
        cls.objects.bulk_create(objs, batch_size=100000)
