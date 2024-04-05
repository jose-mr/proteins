import io
import gzip
import subprocess
import shutil

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
    keywords = models.JSONField(
            default=list
    )
    seq = models.TextField()

    @classmethod
    def create_index(cls, overwrite=False):
        """Creates a index for the dat gz file of UniProt"""
        print(f"Creating index for {DAT_FILE}")
        assert CURRENT_RELEASE_FILE.exists(), f"Cannot access {CURRENT_RELEASE_FILE}"
        if INDEX_RELEASE_FILE.exists() and INDEX_FILE.exists() and not overwrite:
            if INDEX_RELEASE_FILE.read_text() == CURRENT_RELEASE_FILE.read_text():
                print("Updated index file already exists")
        else:
            assert DAT_FILE.exists(), f"Cannot access {DAT_FILE}"
            command = (
                f"""zgrep -b -A1 '^ID' {DAT_FILE} | sed 's/:ID .*/:/' | sed '/^--/d' |"""
                f"""sed 's/^.*-AC\\s*/ /' | awk '/:$/ {{ printf("%s", $0); next }}1' |tr -d " " |"""
                f"""tr  ";" ":" | sed 's/:$//' > {INDEX_FILE}""")
            subprocess.call(command, shell=True)
            shutil.copyfile(CURRENT_RELEASE_FILE, INDEX_RELEASE_FILE)
            print("Index created")

    @classmethod
    def create_and_update_all(cls):
        """Add new UniProt entries"""
        with open(INDEX_FILE, 'r') as ids_file:
            uniprot_ids = set([line.split(":")[1].strip() for line in ids_file])
        # sanity check
        assert len(uniprot_ids) > 550000, f"index file only contains{len(uniprot_ids)} ids"
        # removing obsolete
        existing = set(Entry.objects.all().values_list("ac", flat=True))
        obsolete = existing - uniprot_ids
        if obsolete:
            if input("Found {len(obsolete)} ids. Delete? (y/n") != "y":
                return
            else:
                print(Entry.objects.filter(ac__in=obsolete).delete())
        print(f"Adding or updating {len(uniprot_ids)} UniProt entries")
        cls.create_and_update_list(uniprot_ids)

    @classmethod
    def create_and_update_list(cls, uniprot_ids):
        """Adds the uniprot_ids in the list to the database or updates existing ones

        :param uniprot_ids: a container of uniprot accession numbers
        # :param in_swissprot: if these sequences belong to swissprot or trembl
        """
        # TODO divide in batches to use less ram
        records = cls.read_sequence_records(list(uniprot_ids), INDEX_FILE, DAT_FILE)
        new_sequences = cls.create_and_update_from_records(records)
        print(f"{len(new_sequences)} sequences added to database")
        return new_sequences

    @classmethod
    def read_sequence_records(cls, uniprot_ids, index_file, dat):
        """
        Reads a list of Uniprot ids and returns Bio Seq records

        :param uniprot_ids: the input container with the uniprot ids
        :param index_file: the index file to use
        :param dat: the data file to use
        :return: Bio Seq records
        """
        print("Reading Index")
        index = cls.read_index(index_file, uniprot_ids)
        print("Index Read")
        to_read_tuples = []
        for uniprot_id in uniprot_ids:
            if uniprot_id in index:
                to_read_tuples.append(index[uniprot_id])
        print(f"{len(to_read_tuples)} sequences found in index")
        to_read_tuples = sorted(to_read_tuples, key=lambda tup: tup[0])
        with gzip.open(dat, 'rb') as dat_file:
            previous_end = 0
            data = []
            for no, pair in enumerate(to_read_tuples):
                jump = pair[0]-previous_end
                dat_file.seek(jump, 1)
                data.append(dat_file.read(pair[1]))
                previous_end = sum(pair)
        return SeqIO.parse(io.StringIO(b''.join(data).decode("utf-8")), "swiss")

    @classmethod
    def find_new_sequences(cls, uniprot_ids) -> set:
        """
        Return a set of uniprot ids that are not on the database yet

        :param uniprot_ids: container of accession numbers to be tested
        :return: a set of the accession numbers that are not on the database yet
        """
        return set(uniprot_ids) - set(cls.objects.values_list("ac", flat=True))

    @classmethod
    def create_and_update_from_records(cls, records):
        """
        Creates and updates sequence objects from the record objects
        :param records: Bio Seq record with several sequences
        :param in_swissprot: if the sequences are in swissprot
        :return: the Sequence objects that were added to the database
        """
        to_create = []
        fields_to_check = [
                "name",
                "taxid",
                "comment",
                "seq",
                "keywords",
                "secondary_ac",
                ]

        existing = {entry.ac: entry for entry in cls.objects.all()}

        for record in records:
            name = record.description
            name = name[name.find("Full=")+5:]
            name = name.split(";")[0]
            name = name.split(" {")[0]

            taxid = int(record.annotations['ncbi_taxid'][0])
            comment = record.annotations.get("comment", "")
            seq = record.seq
            keywords = record.annotations.get("keywords", [])
            secondary_ac = record.annotations["accessions"]

            if record.id in existing:
                current_record = existing[record.id]
                needs_update = False
                for field in fields_to_check:
                    new_value = eval(field)
                    if getattr(current_record, field) != new_value:
                        setattr(current_record, field, new_value)
                        needs_update = True
                if needs_update:
                    print("updating ", current_record)
                    current_record.save()
            else:
                to_create.append(cls(
                    ac=record.id,
                    name=name,
                    seq=seq,
                    taxid=taxid,
                    comment=comment,
                    secondary_ac=secondary_ac,
                    keywords=keywords
                ))
        return cls.objects.bulk_create(to_create, batch_size=1000)

    @classmethod
    def read_index(cls, filename, uniprot_ids=None):
        """Read index file with the byte offsets of the uniprot accession numbers

        return only the byte offset of the requested uniprot ids or all by default
        """
        index_dict = {}
        uniprot_ids = set(uniprot_ids) if uniprot_ids is not None else None
        with open(filename, 'r') as index_file:
            previous_entries = []
            for line in index_file:
                words = line.strip().split(":")
                index = int(words[0])
                these_uniprot_ids = set(words[1:])
                for entry in previous_entries:
                    index_dict[entry][1] = index - index_dict[entry][0]
                previous_entries = []
                for uniprot_id in these_uniprot_ids:
                    if uniprot_ids is None or uniprot_id in uniprot_ids:
                        index_dict[uniprot_id] = [index, -1]
                        previous_entries.append(uniprot_id)
        return index_dict
