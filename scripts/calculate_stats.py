
import uniprot.models as uniprot
import go.models as go
import eco.models as eco
from collections import Counter
from pseudoenzymes.settings import OUT_FOLDER

def run():
    venn_enzyme_annotation()
    go_stats()
    uniprot_stats()

def uniprot_stats():
    uniprot_entries = uniprot.Entry.objects.all()
    # print("all swissprot entries", uniprot_entries.distinct().count())

    # # entries with a FUNCTION: section in the description
    # functional = uniprot_entries.filter(comment__contains="FUNCTION:")
    # print("with FUNCTION: description", functional.distinct().count())
    # catalytic = uniprot_entries.catalytic_activity()
    # print("with CATALYTIC ACTIVITY: description", catalytic.distinct().count())
    # print("with CATALYTIC ACTIVITY: and FUNCTION: description", (catalytic & functional).distinct().count())
    # print("with kw", uniprot_entries.enzymes_kw().distinct().count())
    # print("functional with kw", functional.enzymes_kw().distinct().count())
    print("inactive", uniprot_entries.inactive().distinct().count())
    print("caution", uniprot_entries.caution().distinct().count())

    # for entry in uniprot_entries.caution():
        # for caution in entry.cautions:
            # print(entry, caution)
    # for entry in catalytic:
        # for cat in entry.catalytic_activities:
            # print(entry, cat)
    # do all with keyword have a FUNCTION section?

def go_stats():
    go_uniprot = go.TermUniProtEntry.objects.all()
    print("all go-uniprot associations: {} for {} swissprot sequences".format(
          go_uniprot.distinct().count(),
          go_uniprot.values("uniprot_entry_id").distinct().count(),
          ))
    functional = go_uniprot.functional()
    print("functional go-uniprot associations: {} for {} swissprot sequences".format(
          functional.distinct().count(),
          functional.values("uniprot_entry_id").distinct().count(),
          ))
    experimental = go_uniprot.experimental()
    print("experimental go-uniprot associations: {} for {} swissprot sequences".format(
          experimental.distinct().count(),
          experimental.values("uniprot_entry_id").distinct().count(),
          ))
    functional_exp = functional.experimental()
    print("experimental functional go-uniprot associations: {} for {} swissprot sequences".format(
          functional_exp.distinct().count(),
          functional_exp.values("uniprot_entry_id").distinct().count(),
          ))


def venn_enzyme_annotation():
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn3_unweighted as venn3

    proteins = uniprot.Entry.objects.reviewed()
    print("number of proteins in venn diagram", proteins.count())
    ec = set(proteins.enzymes_ec().values_list("ac", flat=True))
    kw = set(proteins.enzymes_kw().values_list("ac", flat=True))
    gos = set(proteins.enzymes_go().values_list("ac", flat=True))
    no_ez_count = proteins.count() - len(ec|kw|gos)
    plt.text(0.5, -0.41, f"No catalytic\nannotation\n{no_ez_count}")
    plt.text(0.5, -0.61, f"Total\n{proteins.count()}")
    
    venn3([ec, kw, gos], ('EC Number', 'UniProt Keyword', 'GO Term'))
    plt.savefig(OUT_FOLDER/"venn_annotation.svg")

    # looking into sets
    ec_kw_not_go = (ec & kw) - gos
    entries = uniprot.Entry.objects.filter(ac__in=ec_kw_not_go)
    # what ecs are these?)
    ec_classes = [n.split(".")[0] for n in  entries.values_list("ec_entries__number", flat=True)]
    print("ec classes of ec and kw but not go", Counter(ec_classes))


    kw_go_not_ec = (gos & kw) - ec
    entries = uniprot.Entry.objects.filter(ac__in=kw_go_not_ec)
    catalytic_kw = [kw for kw in entries.values_list("keywords", flat=True)
                    if kw in uniprot.Keyword.TOP_LEVEL_EZ_KWS]
    go_terms = list(entries.values_list("go_terms__name", flat=True))
    # print("go terms", go_terms)
    all_catalytic_go = set(go.Term.objects.catalytic().values_list("name", flat=True))
    # print(all_catalytic_go)
    catalytic_go = [ term for term in go_terms if term in all_catalytic_go ]
    print("go terms in kw and go but not ec", Counter(catalytic_go))
    # print(list(entries))


    go_only = gos - ec - kw
    entries = uniprot.Entry.objects.filter(ac__in=go_only)
    go_terms = list(entries.values_list("go_terms__name", flat=True))
    catalytic_go = [ term for term in go_terms if term in all_catalytic_go ]
    print("go terms go only", Counter(catalytic_go))
    print(list(entries))



