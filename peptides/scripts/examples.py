
from uniprot.models import Entry



def run():

    peptides = Entry.objects.peptides()
    print("Number of Peptides:", peptides.count())

    # some filters
    # filter by sequence length
    print(peptides.filter(seq_length__gte=10, seq_length__lte=20).count())

    # get a specific peptide
    print(peptides.get(ac="A0A0C5B5G6"))

    # get all reviewed peptides
    # this is how to filter in general
    print(peptides.filter(reviewed=True).count())
    # but in this case I have already made a function for this
    print(peptides.reviewed().count())

    # example peptide info
    example_peptide = peptides.first()
    print("Example peptide:", example_peptide.ac)
    print("Example peptide seq:", example_peptide.seq.seq)
    print("Example peptide sequence length:", example_peptide.seq_length)
    print("Example peptide in swissprot?", example_peptide.reviewed)

    # get all the peptides with a pdb structure
    # important to use distinct when filtering with other tables data
    # otherwise you will get duplicates
    with_structure = peptides.filter(pdb_entries__isnull=False).distinct()
    print("number of peptides with structure: ", with_structure.count())

    # pdbs associated with a peptide
    peptide_with_structure = with_structure.first()
    for pdb in peptide_with_structure.pdb_entries.all():
        print(peptide_with_structure.ac, pdb.pdb_id)







