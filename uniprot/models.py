import subprocess
import shutil

from django.db import models
from pseudoenzymes.settings import DATA_FOLDER

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

    @classmethod
    def create_dat_gz_index(cls, overwrite=False):
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
    def create_and_update_list(cls, uniprot_ids, in_swissprot=True):
        """Adds the uniprot_ids in the list to the database or updates existing ones

        :param uniprot_ids: a container of uniprot accession numbers
        # :param in_swissprot: if these sequences belong to swissprot or trembl
        """
        records = cls.read_sequence_records(list(uniprot_ids), INDEX_FILE, DAT_FILE)
        new_sequences = cls.create_and_update_sequences_from_records(records, in_swissprot=in_swissprot)
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
        index = cls.read_uniprot_index(index_file, uniprot_ids)
        print("Index Read")
        to_read_tuples = []
        for uniprot_id in uniprot_ids:
            if uniprot_id in index:
                to_read_tuples.append(index[uniprot_id])
        print("{} sequences found in index".format(len(to_read_tuples)))
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
    def create_and_update_sequences_from_records(cls, records, in_swissprot: bool):
        """
        Creates and updates sequence objects from the record objects
        :param records: Bio Seq record with several sequences
        :param in_swissprot: if the sequences are in swissprot
        :return: the Sequence objects that were added to the database
        """
        species_ncbi_to_dbid = dict(common.models.Species.objects.get_idx2id())
        to_create = []
        exp_eco = set(common.models.EcoTerm.objects.experimental().values_list("idx", flat=True))
        fields_to_check = ["name", "species_id", "in_reference_proteome", "function_experimental_evidence",
                           "description", "caution", "misc"]
        if in_swissprot:
            existing = cls.objects.reviewed().get_idx2obj()
        else:
            existing = cls.objects.get_idx2obj()

        for no, record in enumerate(records):
            name = record.description
            name = name[name.find("Full=")+5:]
            name = name.split(";")[0]
            name = name.split(" {")[0]
            function_experimental_evidence = False
            description = caution = misc = ""
            if 'comment' in record.annotations:
                for line in record.annotations['comment'].splitlines():
                    if line.startswith('FUNCTION'):
                        description = line[10:line.find("{")]
                        start = 0
                        while True:
                            eco_position = line.find("ECO:", start)
                            if eco_position == -1:
                                break
                            elif line[eco_position+4:eco_position+11] in exp_eco:
                                function_experimental_evidence = True
                                break
                            else:
                                start = eco_position + 1
                    if line.startswith('CAUTION'):
                        caution = line[9:line.find("{")]
                        start = 0
                        while True:
                            eco_position = line.find("ECO:", start)
                            if eco_position == -1:
                                break
                            elif line[eco_position+4:eco_position+11] in exp_eco:
                                function_experimental_evidence = True
                                break
                            else:
                                start = eco_position + 1
                    if line.startswith('MISCELLANEOUS'):
                        start = 0
                        while True:
                            eco_position = line.find("ECO:", start)
                            if eco_position == -1:
                                break
                            elif line[eco_position+4:eco_position+11] in exp_eco:
                                function_experimental_evidence = True
                                break
                            else:
                                start = eco_position + 1
                        misc = line[15:line.find("{")]
            # in reference proteome
            in_reference_proteome = 'Reference proteome' in record.annotations.get('keywords', '')
            taxid = int(record.annotations['ncbi_taxid'][0])
            try:
                species_id = species_ncbi_to_dbid[taxid]
            except KeyError:
                species = common.models.Species.objects.create(taxid=taxid)
                species_id = species.id
                species_ncbi_to_dbid[taxid] = species_id

            if record.id in existing:
                current_record = existing[record.id]
                needs_update = False
                for field in fields_to_check:
                    new_value = eval(field)
                    if getattr(current_record, field) != new_value:
                        print("update field", field, new_value, getattr(current_record, field))
                        setattr(current_record, field, new_value)
                        needs_update = True
                if needs_update:
                    print("updating ", current_record)
                    current_record.save()
            else:
                to_create.append(Sequence(
                    uniprot_id=record.id,
                    name=name,
                    sequence=record.seq,
                    species_id=species_id,
                    in_swissprot=in_swissprot,
                    in_reference_proteome=in_reference_proteome,
                    function_experimental_evidence=function_experimental_evidence,
                    description=description,
                    caution=caution,
                    misc=misc,
                ))
        return Sequence.objects.bulk_create(to_create, batch_size=1000)
