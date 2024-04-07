
from uniprot.models import Entry, Keyword

def run():
    """run this script to generate the pseudoenzyme datasets"""

    # Delete all previous entries and add new entries
    Entry.objects.all().delete()
    Keyword.objects.all().delete()
    Entry.create_from_dat_file()


