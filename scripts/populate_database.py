
from uniprot.models import Entry, Keyword, Sequence
import go.models as go
import eco.models as eco
import ec.models as ec
import cath.models as cath
import wpdb.models as wpdb

def run():
    """run this script to generate the pseudoenzyme datasets"""

    # # Delete all UniProt previous entries
    Entry.objects.all().delete()
    Keyword.objects.all().delete()
    
    # SwissProt
    # Add all swissprot entries and keyword associations
    Entry.create_from_dat_file()

    # PDB
    wpdb.Entry.objects.all().delete()
    wpdb.Entry.download_entries_idx()
    wpdb.Entry.create_from_entries_idx_file()

    # add UniProt in PDBs
    Entry.objects.filter(reviewed=False).delete()
    # Add information about these uniprot entries to the database
    Entry.create_from_pdb_dat_file()
    #  add all uniprot entries associated with a PDB and their keywords
    wpdb.EntryUniProtEntry.download_uniprot_pdb_sifts()
    wpdb.EntryUniProtEntry.objects.all().delete()
    wpdb.EntryUniProtEntry.create_from_uniprot_pdb_sifts()

    # Gene Ontology
    go.Term.download_ontology()
    go.Term.objects.all().delete()
    go.Term.create_from_ontology_file()

    go.Relation.objects.all().delete()
    go.Relation.create_from_ontology_file()


    # eco ontology
    eco.Term.download_ontology()
    eco.Term.objects.all().delete()
    eco.Term.create_from_ontology_file()
    eco.Relation.create_from_ontology_file()


    # link ontologies to uniprot entries
    go.TermUniProtEntry.objects.all().delete()
    go.TermUniProtEntry.create_from_gpa_file()


    # ec numbers
    ec.Entry.download_classes_file()
    ec.Entry.download_ec_dat_file()
    ec.Entry.objects.all().delete()
    ec.Entry.create_classes_classes_files()
    ec.Entry.create_entries_from_dat_file()

    ec.Entry.download_intenz_xml_file()
    ec.Synonym.objects.all().delete()
    ec.Entry.create_synonyms_from_intenz_file()

    # ec <-> uniprot associations
    ec.EntryUniProtEntry.objects.all().delete()
    ec.EntryUniProtEntry.create_from_uniprot_dat_file()
    # TODO add ecs for nonswissprot files
    ec.EntryUniProtEntry.create_from_pdb_uniprot_dat_file()


    # cath
    cath.Superfamily.objects.all().delete()
    cath.Superfamily.download_names_file()
    cath.Superfamily.create_from_names_file()

    # cath <-> uniprot associations
    cath.SuperfamilyUniprotEntry.objects.all().delete()
    cath.SuperfamilyUniprotEntry.create_from_interpro_file()




