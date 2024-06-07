
from collections import defaultdict, Counter
import math
import numpy as np
from scipy.stats import truncnorm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import networkx as nx


import cath.models as cath
import uniprot.models as uniprot
import go.models as go

def run():
    # plot_nseqs_by_pc_single_domain()

    plot_nseqs_by_pc_enzymes(enzyme_list="ec", single_domain_only=True)
    non_catalytic_functions()

def get_cath_to_uniprot(single_domain_only=True):
    if single_domain_only:
        sd_acs = set(cath.SuperfamilyUniprotEntry.objects.get_single_domain_sequence_acs())
    else:
        sd_acs = set()
    cath_to_uniprot = defaultdict(set)
    for cath_number, ac in cath.SuperfamilyUniprotEntry.objects.values_list("superfamily_id", "uniprot_entry_id"):
        if single_domain_only and ac not in sd_acs:
            continue
        cath_to_uniprot[cath_number].add(ac)
    return cath_to_uniprot

def get_cath_to_ecs3():
    cath_to_ecs = defaultdict(set)
    for cath_number, ec in cath.SuperfamilyUniprotEntry.objects.values_list("superfamily_id", "uniprot_entry__ec_entries__number"):
        if ec is not None:
            cath_to_ecs[cath_number].add(ec)

    cath_to_ecs3 = {}
    # remove incomplete ec numbers if there are similar complete ec numbers
    # done at subsubclass level
    for cath_number, ecs in cath_to_ecs.items():
        # ignore anything after the third point to consider only subsubclasses
        ecs = {".".join(ec.split(".")[:3]) for ec in ecs}
        ecs_to_exclude = set()
        for ec in ecs:
            if "-" in ec:
                for other_ec in ecs:
                    if ec != other_ec:
                        incomplete_ec = ec.split("-", 1)[0]
                        if other_ec.startswith(incomplete_ec):
                            ecs_to_exclude.add(ec)
        cath_to_ecs3[cath_number] = ecs - ecs_to_exclude
    return cath_to_ecs3

def plot_nseqs_by_pc_enzymes(enzyme_list, single_domain_only=True):
    """plot the number of sequences in each family and % of enzymes"""
    plt.style.use('_mpl-gallery')

    enzymes = set()
    if enzyme_list == "ec":
        enzymes = set(uniprot.Entry.objects.enzymes_ec().values_list("ac", flat=True))

    cath_to_ec3_count = {k: len(v) for k, v in get_cath_to_ecs3().items()}

    single_domain_acs = cath.SuperfamilyUniprotEntry.objects.get_single_domain_sequence_acs()

    x = []
    y = []
    x0 = []
    y0 = []
    x1 = []
    y1 = []
    ec_color = []
    # ec color for families that only have enzymes
    ec_color1 = []
    number_of_pseudoenzymes = 0
    for cath_number, acs in get_cath_to_uniprot().items():
        # decide color
        if (ec3_count := cath_to_ec3_count.get(cath_number, 0)) == 1:
            color = "grey"
        elif  2 <= (ec3_count) <= 4:
            color = "green"
        elif  5 <= (ec3_count) <= 9:
            color = "orange"
        else:
            color = "red"

        if single_domain_only:
            acs &= single_domain_acs
        if len(acs) > 20:
            number_of_sequences = len(acs)
            pc_enzymes = 100 * (len(acs & enzymes) / number_of_sequences)

            if math.isclose(0, pc_enzymes):
                x0.append(number_of_sequences)
                y0.append(min(-1, np.random.normal(loc=-10, scale=2.5)))
                continue
            if math.isclose(100, pc_enzymes):
                x1.append(number_of_sequences)
                y1.append(max(101, np.random.normal(loc=110, scale=2.5)))
                ec_color1.append(color)
                continue
            if number_of_sequences > 1000:
                print(number_of_sequences, pc_enzymes, cath_number)
            x.append(number_of_sequences)
            number_of_pseudoenzymes += number_of_sequences - len(acs&enzymes)
            y.append(pc_enzymes)
            ec_color.append(color)

    print("number of pseudoenzymes",  number_of_pseudoenzymes)

    print(len(x), len(x0), len(x1))
    fig, ax = plt.subplots()
    fig.set_size_inches(7,5)
    ax.scatter(x, y, alpha=0.5, color=ec_color, s=70)
    ax.scatter(x0, y0, alpha=0.5, color="white", edgecolors="black", s=70)
    ax.scatter(x1, y1, alpha=0.5, color=ec_color1, edgecolors="black", s=70)
    ax.set_yticks(np.arange(0, 100.1, 100/5), [">0%", "20%", "40%", "60%", "80%", "<100%"])
    ax.set_xscale("log")
    ax.set_xlabel("Number of proteins in family")
    ax.set_ylabel("Percentage of enzymes in family")
    plt.tight_layout()
    plt.savefig("a.png", dpi=300)
    plt.show()

def get_cath_family_type(min_proteins=20, single_domain_only=True, enzyme_set="ec"):
    """return a dict with the 3 types of cath families: enzyme, nonenzyme, mixed"""
    cath_type = {
            "enzymes": set(),
            "nonenzymes": set(),
            "mixed": set()
            }
    single_domain_acs = cath.SuperfamilyUniprotEntry.objects.get_single_domain_sequence_acs()

    if enzyme_set == "ec":
        enzymes = set(uniprot.Entry.objects.enzymes_ec().values_list("ac", flat=True))
    elif enzyme_set == "go":
        enzymes = set(uniprot.Entry.objects.enzymes_go().values_list("ac", flat=True))

    for cath_number, acs in get_cath_to_uniprot(single_domain_only=single_domain_only).items():
        if single_domain_only:
            acs &= single_domain_acs
        if len(acs) >= min_proteins:
            number_of_sequences = len(acs)
            pc_enzymes = 100 * (len(acs & enzymes) / number_of_sequences)
            if math.isclose(0, pc_enzymes):
                cath_type["nonenzymes"].add(cath_number)
                for ac in acs:
                    if ac in enzymes:
                        print(ac)
            elif math.isclose(100, pc_enzymes):
                cath_type["enzymes"].add(cath_number)
            else:
                cath_type["mixed"].add(cath_number)
    return cath_type



def plot_nseqs_by_pc_single_domain():
    """plot the number of sequences in each family and % of single domain sequences"""
    import matplotlib.pyplot as plt
    plt.style.use('_mpl-gallery')

    single_domain_acs = cath.SuperfamilyUniprotEntry.objects.get_single_domain_sequence_acs()

    x = []
    y = []
    x0 = []
    y0 = []
    x1 = []
    y1 = []

    for _, acs in get_cath_to_uniprot().items():
        if len(acs) > 100:
            number_of_sequences = len(acs)
            pc_single_domain_sequences = len(acs & single_domain_acs) / number_of_sequences

            if math.isclose(0, pc_single_domain_sequences):
                x0.append(number_of_sequences)
                y0.append(min(-0.01, np.random.normal(loc=-0.1, scale=0.03)))
                continue
            if math.isclose(1, pc_single_domain_sequences):
                x1.append(number_of_sequences)
                y1.append(max(1.01, np.random.normal(loc=1.1, scale=0.03)))
                continue
            x.append(number_of_sequences)
            y.append(pc_single_domain_sequences)

    fig, ax = plt.subplots()
    ax.scatter(x, y, alpha=0.2)
    ax.scatter(x0, y0, alpha=0.2, color="red")
    ax.scatter(x1, y1, alpha=0.2, color="red")
    ax.set_xscale("log")
    plt.show()




def non_catalytic_functions():
    # TODO needs refactoring, was rushing to finish ppt

    go_id_to_term = {t.id: t for t in go.Term.objects.all()}

    ac_to_gos = defaultdict(set)
    values = go.TermUniProtEntry.objects.functional().filter(qualifier="enables").values_list("uniprot_entry", "term__id")
    for ac, term_name in values:
        ac_to_gos[ac].add(term_name)

    go_to_ancestors = go.Term.objects.go_to_ancestors()

    protein_to_gos = defaultdict(set)
    for cath, proteins in get_cath_to_uniprot().items():
        if len(proteins) > 10:
            for protein in proteins:
                go_terms = ac_to_gos[protein]
                for go_term in go_terms:
                    protein_to_gos[protein].update(go_to_ancestors[go_term])

    # gos_to_proteins = defaultdict(set)
    # for protein, gos in protein_to_gos.items():
        # for go_id in gos:
            # gos_to_proteins[go_id].add(protein)

    # gos_count = {go_id: len(proteins) for go_id, proteins in gos_to_proteins.items()}
    # gos_count = {go_id: gos_count[go_id] for go_id in sorted(gos_count, key= lambda x: gos_count[x])}

    cath_to_uniprot = get_cath_to_uniprot()

    enzymes = set(uniprot.Entry.objects.enzymes_go().values_list("ac", flat=True))
    fig, axs = plt.subplots(1, 4)
    fig.set_size_inches(9,5)

    gos_to_cath = {}
    cath_family_type = get_cath_family_type(enzyme_set="go")
    cath_family_type["mixed_enzymes"] = cath_family_type["mixed"]
    cath_family_type["mixed_nonenzymes"] = cath_family_type["mixed"]

    catalytic_go_ids = set(go.Term.objects.catalytic().values_list("id", flat=True))

    go_to_type_to_pc = defaultdict(dict)
    for plot, (cath_type, cath_numbers) in enumerate(cath_family_type.items()):
        gos_to_cath[cath_type] = defaultdict(set)

        for cath_number in cath_numbers:
            for protein in cath_to_uniprot[cath_number]:
                if protein in enzymes and cath_type == "mixed_nonenzymes":
                    continue
                if protein not in enzymes and cath_type == "mixed_enzymes":
                    continue
                # if protein in enzymes and cath_type == "nonenzymes":
                    # print(protein, cath_number)
                go_terms = protein_to_gos[protein]
                for go_term in go_terms:
                    gos_to_cath[cath_type][go_term].add(cath_number)

    
        for go_id, caths in gos_to_cath[cath_type].items():
            pc = len(caths) / len(cath_numbers)
            go_to_type_to_pc[go_id][cath_type] = pc

    gos_to_cath.pop("mixed")

    max_gos = {}
    for go_id, inner in go_to_type_to_pc.items():
        for cath_type, pc in inner.items():
            max_gos[go_id] = max(max_gos.get(go_id, 0), pc)

    graph = go.Term.objects.functional().go_graph()
    
    dfs = nx.dfs_tree(graph, 3674, depth_limit=2,
                      sort_neighbors=lambda x: sorted(x, key=lambda i: max_gos.get(i, 0), reverse=True))
    # print(dfs)

    xs = []
    ys = [[], [], [], []]

    cath_types = list(gos_to_cath.keys())
    # print(cath_types)

    for go_id in dfs:
        if go_id == 3674:
            continue
        if max_gos.get(go_id, 0) > 0.2:
            xs.insert(0, go_id_to_term[go_id].name)
            for i, cath_type in enumerate(cath_types):
                ys[i].insert(0, go_to_type_to_pc.get(go_id, {}).get(cath_type, 0))

    for i, cath_type in enumerate(cath_types):
        axs[i].barh(xs, ys[i])
        axs[i].set_xlabel(cath_type)

        axs[i].label_outer()
        axs[i].set_xticks(np.arange(0, 1.001, 0.5))

    # yax = axs[0].get_yaxis()
    # yax.set_tick_params(pad=50)

    # fig.align_ylabels()
    # fig.align_xlabels()
    plt.tight_layout()
    plt.show()
