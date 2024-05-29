
from collections import defaultdict, Counter
import math
import numpy as np
from scipy.stats import truncnorm
import matplotlib.colors as mcolors

import cath.models as cath
import uniprot.models as uniprot

def run():
    # plot_nseqs_by_pc_single_domain()

    plot_nseqs_by_pc_enzymes(enzyme_list="ec", single_domain_only=True)

def plot_nseqs_by_pc_enzymes(enzyme_list, single_domain_only=True):
    """plot the number of sequences in each family and % of enzymes"""
    import matplotlib.pyplot as plt
    plt.style.use('_mpl-gallery')

    enzymes = set()
    if enzyme_list == "ec":
        enzymes = set(uniprot.Entry.objects.enzymes_ec().values_list("ac", flat=True))

    cath_to_uniprot = defaultdict(set)
    for cath_number, ac in cath.SuperfamilyUniprotEntry.objects.values_list("superfamily_id", "uniprot_entry_id"):
        cath_to_uniprot[cath_number].add(ac)

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

    cath_to_ec3_count = {k: len(v) for k, v in cath_to_ecs3.items()}

    # print(sorted(cath_to_ecs, key = lambda x: len(cath_to_ecs[x])))
    # print(cath_to_ecs["3.40.50.150"])
    print(Counter(cath_to_ec3_count.values()))

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
    for cath_number, acs in cath_to_uniprot.items():
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
            x.append(number_of_sequences)
            y.append(pc_enzymes)
            ec_color.append(color)

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



def plot_nseqs_by_pc_single_domain():
    """plot the number of sequences in each family and % of single domain sequences"""
    import matplotlib.pyplot as plt
    plt.style.use('_mpl-gallery')

    cath_to_uniprot = defaultdict(set)
    for cath_number, ac in cath.SuperfamilyUniprotEntry.objects.values_list("superfamily_id", "uniprot_entry_id"):
        cath_to_uniprot[cath_number].add(ac)
    single_domain_acs = cath.SuperfamilyUniprotEntry.objects.get_single_domain_sequence_acs()

    x = []
    y = []
    x0 = []
    y0 = []
    x1 = []
    y1 = []
    for cath_number, acs in cath_to_uniprot.items():
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




