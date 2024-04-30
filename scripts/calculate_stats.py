
import uniprot.models as uniprot
import go.models as go
from collections import Counter
from pseudoenzymes.settings import OUT_FOLDER

def run():

    go_stats()
    venn_enzyme_annotation()

def go_stats():
    # try to understand why go has so many catalytic associations more 
    # than ec and go

    go_terms = go.Term.objects.all()
    go_uniprot = go.TermUniProtEntry.objects.all()
    print("go uniprot associations", go_uniprot.distinct().count())
    # check enable and not|enable
    catalytic_go_uniprot = go_uniprot.filter(term__in=go_terms.catalytic())
    print("catalytic", catalytic_go_uniprot.distinct().count())
    print("catalytic qualifiers",
          Counter(catalytic_go_uniprot.values_list("qualifier")))
    catalytic_association = go.TermUniProtEntry.objects.filter(
                term__in=go.Term.objects.catalytic(),
                qualifier="enables"
                )
    print("catalytic enable",
          Counter(catalytic_association.values_list("qualifier")))


    # check experimental and not experimental

def venn_enzyme_annotation():
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn3_unweighted as venn3

    proteins = uniprot.Entry.objects.all()
    ec = set(proteins.enzymes_ec().values_list("ac", flat=True))
    kw = set(proteins.enzymes_kw().values_list("ac", flat=True))
    go = set(proteins.enzymes_go().values_list("ac", flat=True))
    no_ez_count = proteins.count() - len(ec|kw|go)
    plt.text(0.5, -0.41, f"No catalytic\nannotation\n{no_ez_count}")
    plt.text(0.5, -0.61, f"Total\n{proteins.count()}")
    
    venn3([ec, kw, go], ('EC Number', 'Uniprot Keyword', 'Go Term'))
    plt.savefig(OUT_FOLDER/"venn_annotation.svg")
