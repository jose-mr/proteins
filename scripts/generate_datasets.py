
from uniprot.models import Entry

def run():
    """run this script to generate the pseudoenzyme datasets"""

    # update uniprot index file
    Entry.create_dat_gz_index()

