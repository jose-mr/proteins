
from django.db import transaction
from django.db.models import F
from uniprot.models import Entry, Keyword, Sequence
import go.models as go
import eco.models as eco
import ec.models as ec
import cath.models as cath
import wpdb.models as wpdb
import taxonomy.models as taxonomy
from pseudoenzymes.settings import UNIPROT_DAT_FILE, SWISSPROT_DAT_FILE


def run():
    """run this script to generate the pseudoenzyme datasets"""
    # taxonomy
    print("Adding taxonony data")
    print("downloading")
    # taxonomy.Taxon.download_taxdump_from_ncbi()
    print("add taxonomy data to DB")
    # taxonomy.Taxon.create_from_ncbi_files()

    print("Adding Swiss-prot data")
    # Entry.create_from_dat_file(SWISSPROT_DAT_FILE)

    print("Adding trembl UniProt data")
    # use skip_first if restarting
    # Entry.create_from_dat_file(UNIPROT_DAT_FILE, skip_first=0)

    # PDB
    # wpdb.Entry.download_entries_idx()
    # wpdb.Entry.create_from_entries_idx_file()

    # link pdb to uniprot
    # wpdb.EntryUniProtEntry.download_uniprot_pdb_sifts()
    # wpdb.EntryUniProtEntry.create_from_uniprot_pdb_sifts()

    # Gene Ontology
    # go.Term.download_ontology()
    # go.Term.create_from_ontology_file()
    # go.Relation.create_from_ontology_file()

    # # eco ontology
    # eco.Term.download_ontology()
    # eco.Term.create_from_ontology_file()
    # eco.Relation.create_from_ontology_file()

    # link ontologies to uniprot entries
    # go.TermUniProtEntry.prepare_csv_file_for_ingestion()
    go.TermUniProtEntry.ingest_csv_file()
    # go.TermUniProtEntry.create_from_gpa_file()


    # TODO go here




    # # ec numbers
    # ec.Entry.download_classes_file()
    # ec.Entry.download_ec_dat_file()
    # ec.Entry.objects.all().delete()
    # ec.Entry.create_classes_classes_files()
    # ec.Entry.create_entries_from_dat_file()

    # ec.Entry.download_intenz_xml_file()
    # ec.Synonym.objects.all().delete()
    # ec.Entry.create_synonyms_from_intenz_file()

    # # ec <-> uniprot associations
    # ec.EntryUniProtEntry.objects.all().delete()
    # ec.EntryUniProtEntry.create_from_uniprot_dat_file()
    # # TODO add ecs for nonswissprot files
    # ec.EntryUniProtEntry.create_from_pdb_uniprot_dat_file()


    # # cath
    # cath.Superfamily.objects.all().delete()
    # cath.Superfamily.download_names_file()
    # cath.Superfamily.create_from_names_file()

    # # cath <-> uniprot associations
    # cath.SuperfamilyUniprotEntry.objects.all().delete()
    # cath.SuperfamilyUniprotEntry.create_from_interpro_file()




