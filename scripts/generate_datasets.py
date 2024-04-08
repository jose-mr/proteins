
from uniprot.models import Entry, Keyword
import go.models as go

def run():
    """run this script to generate the pseudoenzyme datasets"""

    # Delete all UniProt previous entries and add new entries
    # Entry.objects.all().delete()
    # Keyword.objects.all().delete()
    # Entry.create_from_dat_file()

    # Gene Ontology
    go.Term.download_ontology()
    go.Term.objects.all().delete()
    go.Term.create_from_ontology_file()

    go.Relation.objects.all().delete()
    go.Relation.create_from_ontology_file()

    

    



