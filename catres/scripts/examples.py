
from uniprot.models import Entry
from taxonomy.models import Taxon



def run():

    protein_list = ["A1L3X0", "A2RUC4", "O00165"]

    proteins = Entry.objects.filter(ac__in=protein_list)
    print("Number of proteins:", proteins.count())
    for p in proteins:
        print(p)
        # print(p.comment)
        print(p.seq.seq)
        print(p.species)
        print(p.species.taxid)
        print(p.cath_superfamilies.all())

    print(Entry.objects.reviewed().count())
    return








