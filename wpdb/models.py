from django.db import models
import gzip
from urllib.request import urlretrieve

from pseudoenzymes.settings import PDB_ENTRIES_IDX, PDB_UNIPROT_SIFTS, PDB_UNIPROT_DAT_FILE
import uniprot.models as uniprot

class Entry(models.Model):
    pdb_id = models.CharField(primary_key=True, max_length=12)
    title = models.TextField()
    date = models.DateField()
    method = models.TextField()
    # resolution = models.CharField(max_length=256)
    # aa_sequences = models.ManyToManyField('AASequence', through='Chain', related_name="structures")
    uniprot_entries = models.ManyToManyField(
            "uniprot.Entry",
            related_name="pdb_entries",
            through="EntryUniProtEntry"
            )
    # biological_assembly = models.IntegerField(null=True)


    # objects = StructureQuerySet.as_manager()

    def __str__(self):
        return self.pdb_id

    @classmethod
    def download_entries_idx(cls):
        """download entries.idx file"""
        print("Start downloading PDB entries.idx file")
        urlretrieve("https://ftp.ebi.ac.uk/pub/databases/pdb/derived_data/index/entries.idx",
                    PDB_ENTRIES_IDX)

    @classmethod
    def create_from_entries_idx_file(cls):
        """Reads entries.idx file and adds new pdb entries to the table"""
        existing = set(cls.objects.values_list("pdb_id", flat=True))
        to_create = []
        with open(PDB_ENTRIES_IDX, 'r') as pdb_idx_file:
            next(pdb_idx_file)
            next(pdb_idx_file)
            for line in pdb_idx_file:
                words = line.strip().split("\t")
                pdb_id = words[0].lower()
                month, day, year = words[2].split("/")
                year_prefix = "19" if year[0] in ("7", "8", "9") else "20"
                date = f"{year_prefix}{year}-{month}-{day}"
                title = words[3].lower()
                method = words[7].lower()
                if pdb_id not in existing:
                    to_create.append(cls(
                        pdb_id=pdb_id, title=title, date=date, method=method))

        print(f"Addding {len(to_create)} new PDB entries")
        cls.objects.bulk_create(to_create)



    # @classmethod
    # def link_to_sequences(cls):
        # """Populate m2m structure_sequences table"""
        # existing = set(cls.objects.values_list('pdb_id', 'sequences__uniprot_id'))
        # nones = set(Structure.objects.values_list('pdb_id', 'sequences__uniprot_id').filter(sequences=None))
        # print('Structures that are not associated with any Uniprot Sequence: {}'.format(len(nones)))
        # in_file = set()
        # with gzip.open(files.SIFTS_UNIPROT_PDB, 'rt') as uniprot_pdb_file:
            # next(uniprot_pdb_file)
            # next(uniprot_pdb_file)
            # for line in uniprot_pdb_file:
                # uniprot_id, pdb_ids = line.split("\t")
                # pdb_ids = pdb_ids.strip('\n').split(";")
                # for pdb_id in pdb_ids:
                    # in_file.add((pdb_id, uniprot_id))
        # to_add = in_file - existing
        # # todo delete old ones
        # print("Adding {} Structure-Sequence links".format(len(to_add)))
        # # this could be made faster by adding to the intermediate table directly
        # # see Chain.link_to_sequences to see how
        # for (pdb_id, uniprot_id) in to_add:
            # try:
                # sequence = common.models.Sequence.objects.get(uniprot_id=uniprot_id)
            # except common.models.Sequence.DoesNotExist:
                # print("Sequence not found in database: {0} (maybe deleted from uniprot)".format(uniprot_id))
                # continue
            # structure = cls.objects.get(pdb_id=pdb_id)
            # sequence.structures.add(structure)

    # def get_sifts_file(self):
        # """Returns the pdb file name. Downloads it if necessary"""
        # sub_directory = files.SIFTS_FILES_PATH / "{}/".format(self.pdb_id[1:3])
        # sub_directory.mkdir(exist_ok=True)
        # sifts_filename = sub_directory / "{}.xml.gz".format(self.pdb_id)
        # if not sifts_filename.exists():
            # url = "ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/{}/{}.xml.gz".format(self.pdb_id[1:3], self.pdb_id)
            # print(url)
            # try:
                # urlretrieve(url, sifts_filename)
                # print("Downloaded pdb_id sifts file {} to {}.".format(url, sifts_filename))
            # except URLError:
                # print("Could not download sifts file.".format(url))
                # print("SIFTS XML file not found: {}. Obsolete entry or need to update the devserver files".format(sifts_filename))
                # return None
        # return sifts_filename

    # def get_sequence_resid_chains_dict(self, output_key="uniprot_resid", use_assembly_chains=True):
        # """Returns a dict where keys are Uniprot res and values are pdb res chains"""
        # sifts_filename = self.get_sifts_file()
        # if sifts_filename is None:
            # return None
        # with gzip.open(sifts_filename) as sifts_file:
            # xml_str = sifts_file.read()
        # xml_str = re.sub(b'\sxmlns="[^"]+"', b'', xml_str, count=1)
        # root = etree.fromstring(xml_str)
        # sequence_chain_dict = collections.defaultdict(list)
        # chain_resid_to_auth_dict = {}
        # for entity in root.findall(".//entity"):
            # assembly_chains = set()
            # for chain_group in self.assembly_chains_groups.values():
                # if entity.attrib['entityId'] in chain_group:
                    # assembly_chains.update(chain_group)
            # for residue in entity.findall(".//residue"):
                # uniprot_cross_rf = residue.find("./crossRefDb[@dbSource='UniProt']")
                # if uniprot_cross_rf is None:
                    # continue
                # seq_resid = int(residue.attrib['dbResNum'])
                # pdb_cross_ref = residue.find("./crossRefDb[@dbSource='PDB']")
                # auth_resid = pdb_cross_ref.get('dbResNum')
                # # in some rare cases auth_resid is not a number
                # auth_resid = int(re.sub('[^0-9]', '', auth_resid)) if auth_resid != 'null' else None
                # pdb_resname = pdb_cross_ref.get('dbResName').capitalize()
                # u_resid = int(uniprot_cross_rf.attrib['dbResNum'])
                # uniprot_id = uniprot_cross_rf.attrib['dbAccessionId']
                # if not use_assembly_chains:
                    # chains = {residue.find("./crossRefDb[@dbSource='PDB']").get('dbChainId'), }
                # else:
                    # chains = assembly_chains
                # if (uniprot_id, u_resid) in sequence_chain_dict:
                    # sequence_chain_dict[(uniprot_id, u_resid)][0].update(chains)
                # else:
                    # sequence_chain_dict[(uniprot_id, u_resid)] = [set(chains), seq_resid, auth_resid, pdb_resname]
                # for chain in chains:
                    # chain_resid_to_auth_dict[(chain, seq_resid)] = auth_resid

        # return sequence_chain_dict if output_key == "uniprot_resid" else chain_resid_to_auth_dict

    # @cached_property
    # def sequence_resid_chains_dict(self):
        # return self.get_sequence_resid_chains_dict(output_key="uniprot_resid", use_assembly_chains=True)

    # @cached_property
    # def assembly_chain_to_pdb_chain(self) -> Dict[str, str]:
        # """return map of the auth_asym_id in mmmcif files to the label_asym_id"""
        # data = self.mmcif_dict
        # # return {label: data["_atom_site.auth_asym_id"][i] for i, label in enumerate(data["_atom_site.label_asym_id"])}
        # return {label: data["_pdbe_chain_remapping.orig_auth_asym_id"][i]
                # for i, label in enumerate(data["_pdbe_chain_remapping.new_label_asym_id"])}

    # @property
    # def pdb_chain_to_assembly_chains(self):
        # pdb_to_assembly = {}
        # for assembly_chain, pdb_chain in self.assembly_chain_to_pdb_chain.items():
            # if pdb_chain in pdb_to_assembly:
                # pdb_to_assembly[pdb_chain].append(assembly_chain)
            # else:
                # pdb_to_assembly[pdb_chain] = [assembly_chain, ]
        # return pdb_to_assembly

    # @cached_property
    # def assembly_chains_groups(self):
        # entity_id_to_assembly_chains = collections.defaultdict(list)
        # for entity_id, asym_id in zip(
            # self.mmcif_dict["_pdbe_chain_remapping.orig_label_asym_id"],
            # self.mmcif_dict["_pdbe_chain_remapping.new_auth_asym_id"],
        # ):
            # entity_id_to_assembly_chains[entity_id].append(asym_id)
        # return dict(entity_id_to_assembly_chains)

    # def get_assembly_xml(self):
        # if self._assembly_xml_tree is None:
            # url = ("http://www.ebi.ac.uk/pdbe/static/entry/download/{}-assembly.xml".format(self.pdb_id))
            # xml_str = urlopen(url).read()
            # xml_str = re.sub(b'\sxmlns="[^"]+"', b'', xml_str, count=1)
            # self._assembly_xml_tree = etree.fromstring(xml_str)
        # return self._assembly_xml_tree

    # def get_default_biological_assembly(self):
        # return int(self.get_assembly_xml().find('./assembly[@prefered="True"]').attrib['id'])

    # # def get_pdb_file(self):
        # # """Returns the pdb file name. Downloads it if necessary"""
        # # pdb_filename = pathlib.Path("{0}/pdb{1.pdb_id}.ent".format(files.PDB_FILES_PATH, self))
        # # if not pdb_filename.exists():
            # # Bio.PDB.PDBList().retrieve_pdb_file(self.pdb_id, pdir=files.PDB_FILES_PATH, file_format="pdb")
        # # return pdb_filename

    # # def get_pdb_parser(self):
        # # """Get Biopython parser for pdb file"""
        # # parser = Bio.PDB.PDBParser()
        # # return parser.get_structure('pdb', self.get_pdb_file())

    # @property
    # def mmcif_dict(self):
        # """return mmcif data as a dict. Just a wrapper to the MMCIF2Dict biopython function"""
        # cif_file = gzip.open(self.get_assembly_cif_file(), mode="rt", encoding='utf-8')
        # mmcif_dict = MMCIF2Dict(cif_file)
        # cif_file.close()
        # return mmcif_dict

    # # this does not return the same objects as io
    # # wip
    # @cached_property
    # def assembly_parser(self):
        # """return the Biopython parser for the mmcif file of the biological assembly"""
        # parser = MMCIFParser( QUIET=True)
        # cif_file = gzip.open(self.get_assembly_cif_file(), mode="rt", encoding='utf-8')
        # print(self.get_assembly_cif_file())
        # assembly_parser = parser.get_structure('mmcif-{}'.format(self.pdb_id), cif_file)
        # cif_file.close()
        # return assembly_parser
 
    # def residue_min_distance(self, res_a, res_b):
        # """returns the min distance between the two residues. each residue is defined as (resid, chain) tuple"""
        # minimum = 999
        # for atom_a in res_a:
            # for atom_b in res_b:
                # minimum = min(minimum, atom_a-atom_b)
        # return minimum

    # @cached_property
    # def assembly_het_residues(self):
        # """return list of Biopython residue objects relative to HET residues in the structure"""
        # het_residues = []
        # for model in self.assembly_parser:
            # for chain in model:
                # for residue in chain:
                    # if residue.get_id()[0].startswith('H') or residue.get_id()[0] == "W":
                        # het_residues.append(residue)
        # return het_residues

    # @cached_property
    # def assembly_het_labels_to_residues(self):
        # return {self.get_label_from_residue(res): res for res in self.assembly_het_residues}

    # @cached_property
    # def assembly_het_labels_to_orig_chain_het_labels(self):
        # return {self.get_label_from_residue(res, orig_chain=False): self.get_label_from_residue(res, orig_chain=True)
                # for res in self.assembly_het_residues}

    # @cached_property
    # def assembly_orig_chain_het_labels_to_het_labels(self):
        # return {v: k for k, v in self.assembly_het_labels_to_orig_chain_het_labels.items()}

    # @cached_property
    # def assembly_het_labels_to_mols(self):
        # labels_to_mols = {}
        # for label, residue in self.assembly_het_labels_to_residues.items():
            # mol = self.get_mol_from_pdb_residue(residue)
            # # Chem.rdDepictor.Compute2DCoords(mol)
            # labels_to_mols[label] = mol
        # return labels_to_mols

    # def get_label_from_residue(self, residue, orig_chain=False):
        # """a label usefull to identify residues"""
        # parent_id = self.assembly_chain_to_pdb_chain[residue.parent.id] if orig_chain else residue.parent.id
        # return "{},{},{},{}".format(parent_id, *residue.id)

    # def get_residue_from_label(self, pdb_residue_label):
        # """get the pdb residue for this label"""
        # words = re.split(",", pdb_residue_label)
        # print(words)
        # chain = words[0]
        # resid = (words[1], int(words[2]), words[3])
        # return self.assembly_parser[0][chain][resid]

    # def get_residue_center_of_geometry(self, residue):
        # """return the coordinates corresponding to the center of geometry"""
        # return sum([atom.coord for atom in residue])/len(residue)

    # def get_same_assembly_chains(self, chain_name):
        # """returns a list of the assembly chains that are based on the same original pdb chain"""
        # for assembly_chain_group in self.assembly_chains_groups.values():
            # if chain_name in assembly_chain_group:
                # return assembly_chain_group
        # return None

    # def get_assembly_cif_file(self, assembly=None, overwrite=False):
        # """Returns the mmcif filename and downloads it if necessary"""
        # if assembly is None:
            # if self.biological_assembly is None:
                # self.biological_assembly = self.get_default_biological_assembly()
                # self.save()
            # assembly = self.biological_assembly

        # directory = files.MMCIF_FILES_PATH / "{}/{}/".format(self.pdb_id[1:3], self.pdb_id)
        # directory.mkdir(exist_ok=True, parents=True)
        # filename = directory / "{}-assembly-{}.cif.gz".format(self.pdb_id, assembly)
        # if not filename.exists() or overwrite: # this should only happen in the dev server (run setup_dev_files)
            # url = "https://www.ebi.ac.uk/pdbe/static/entry/download/{}-assembly-{}.cif.gz".format(self.pdb_id, assembly)
            # try:
                # urlretrieve(url, filename)
                # print("Downloaded pdb_id: {}, assembly: {} to {}.".format(self.pdb_id, assembly, filename))
            # except HTTPError:
                # print("Could not download pdb_id: {}, assembly: {} from {}.".format(self.pdb_id, assembly, url))
        # if filename.exists():
            # return filename
    
    # def move_assembly_file_to_media(self):
        # mmcif_media_dir = pathlib.Path(settings.MEDIA_ROOT) / "mmcif/"
        # mmcif_filename = mmcif_media_dir / "{0.pdb_id}-assembly-{0.biological_assembly}.cif".format(self)
        # if not mmcif_filename.exists():
            # try:
                # os.system("cd {0}; cp {1} .; gzip -d {1.stem}".format(mmcif_media_dir, self.get_assembly_cif_file()))
            # except AttributeError as e:
                # print(e)

    # @property
    # def cif_file_url(self):
        # self.move_assembly_file_to_media()
        # return "{0}mmcif/{1}-assembly-{2}.cif".format(settings.MEDIA_URL, self.pdb_id, self.biological_assembly)

    # @cached_property
    # def biopython_pdb_io(self):
        # pdb_io = PDBIO()
        # pdb_io.set_structure(self.assembly_parser)
        # return pdb_io


    # # def save_pdb_with_residue_selection(self, pdb_residues, new_filename):
        # # tmp_folder = pathlib.Path("/tmp/")
        # # tmp_filename = tmp_folder / 'test.pdb' #  .format(pdb_residue.id, pdb_residue.parent.id)

        # # residues_id_chain = set()
        # # old_chains = {}
        # # # to avoid problems with chains with dashes
        # # if (z_chains := [chain for chain in self.assembly_parser[0] if chain.get_id() == "z"]):
            # # new_chain = z_chains[0]
        # # else:
            # # new_chain = BioChain("z")
            # # self.assembly_parser[0].add(new_chain)

        # # for pdb_residue in pdb_residues:
            # # print("a", pdb_residue, id(pdb_residue))
            # # old_chains[id(pdb_residue)] = pdb_residue.parent
            # # old_chains[id(pdb_residue)].detach_child(pdb_residue.id)
            # # new_chain.add(pdb_residue)
            # # residues_id_chain.add((pdb_residue.id, pdb_residue.parent.id))

        # # for pdb_residue, chain in old_chains.items():
            # # print(pdb_residue, id(pdb_residue), chain)

        # # self.biopython_pdb_io.save(str(tmp_filename), ResidueSelect(residues_id_chain))

        # # for pdb_residue in pdb_residues:
            # # print("b", pdb_residue, id(pdb_residue))
            # # new_chain.detach_child(pdb_residue.id)
            # # old_chains[id(pdb_residue)].add(pdb_residue)

        # # self.assembly_parser[0].detach_child(new_chain.id)


    # _mol_from_pdb_residue_data = {}
    # def get_mol_from_pdb_residue(self, pdb_residue: PdbResidue, atom_selection="") -> Chem.Mol:
        # """return a rdkit mol file for this pdb_residue"""
        # if (pdb_residue.get_full_id(), atom_selection) not in self._mol_from_pdb_residue_data:
            # tmp_folder = pathlib.Path("/tmp/")
            # tmp_filename = tmp_folder / '{}_{}.pdb'.format(pdb_residue.id, pdb_residue.parent.id)
            # # to avoid problems with chains with dashes
            # if (z_chains := [chain for chain in pdb_residue.parent.parent if chain.get_id() == "z"]):
                # new_chain = z_chains[0]
            # else:
                # new_chain = BioChain("z")
                # pdb_residue.parent.parent.add(new_chain)
            # old_chain = pdb_residue.parent
            # old_chain.detach_child(pdb_residue.id)
            # new_chain.add(pdb_residue)

            # # print(pdb_residue.id, pdb_residue.parent.id, old_chain.id, atom_selection)
            # self.biopython_pdb_io.save(str(tmp_filename), ResidueSelect(
                # {(pdb_residue.id, pdb_residue.parent.id)}, atom_selection=atom_selection))

            # new_chain.detach_child(pdb_residue.id)
            # old_chain.add(pdb_residue)
            # # ###

            # pdb_residue_mol = Chem.MolFromPDBFile(str(tmp_filename), sanitize=False)
            # self._mol_from_pdb_residue_data[(pdb_residue.get_full_id(), atom_selection)] = pdb_residue_mol

        # return self._mol_from_pdb_residue_data[(pdb_residue.get_full_id(), atom_selection)]

class EntryUniProtEntry(models.Model):
    """through table that links pdb to uniprot entries"""
    uniprot_entry = models.ForeignKey('uniprot.Entry', on_delete=models.CASCADE, 
                                      related_name="pdb_entries_through")
    pdb_entry = models.ForeignKey("Entry", on_delete=models.CASCADE,
                                  related_name="uniprot_entries_through")

    class Meta:
        unique_together = ("uniprot_entry", "pdb_entry")

    @classmethod
    def download_uniprot_pdb_sifts(cls):
        """download uniprot_pdb.tsv.gz from sifts"""
        urlretrieve("https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/"
                    "uniprot_pdb.tsv.gz", PDB_UNIPROT_SIFTS)

    @classmethod
    def create_from_uniprot_pdb_sifts(cls):
        """Add uniprot - pdb associations"""
        existing = set(cls.objects.values_list("uniprot_entry_id", "pdb_entry_id"))
        ac_to_pdb = {}
        with gzip.open(PDB_UNIPROT_SIFTS, 'rt') as uniprot_pdb_file:
            next(uniprot_pdb_file)
            next(uniprot_pdb_file)
            for line in uniprot_pdb_file:
                words = line.strip().split("\t")
                ac_to_pdb[words[0]] = words[1].split(";")


        db_acs = set(uniprot.Entry.objects.values_list("pk", flat=True))
        in_sifts = set()
        missing_uniprot = ac_to_pdb.keys() - db_acs
        # newer uniprot ids might not be in the uniprot dat file yet
        # ading them here
        print(f"{len(missing_uniprot)} missing UniProt")
        uniprot.Entry.create_from_ac_list(missing_uniprot)

        for ac, pdbs in ac_to_pdb.items():
            if ac in db_acs:
                for pdb_id in pdbs:
                    in_sifts.add((ac, pdb_id))
        to_create = in_sifts - existing
        objs_to_create = [cls(uniprot_entry_id=ac, pdb_entry_id=pdb_id) for ac, pdb_id in to_create]
        print(f"Creating {len(to_create)} pdb-uniprot associations")
        cls.objects.bulk_create(objs_to_create)
