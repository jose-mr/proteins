
from uniprot.models import Entry

def run():
    """run this script to generate the pseudoenzyme datasets"""

    # update uniprot index file
    ## create index
    Entry.create_dat_gz_index()
    
    # Add and update entries
    Entry.create_and_update_all()

