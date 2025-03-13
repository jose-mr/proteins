
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






