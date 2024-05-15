from collections import defaultdict
from multiprocessing.pool import ThreadPool
import os
import subprocess
from pathlib import Path
from textwrap import wrap

from ete3 import NCBITaxa
ncbi = NCBITaxa()
from ete3 import PhyloTree as etetree

import cath.models as cath
import uniprot.models as uniprot
from pseudoenzymes.settings import (MSA_FOLDER, MSA_BY_DOMAIN_STRIP,
                                    TREE_OUT, NOTUNG_FOLDER, NOTUNG_JAR)

def run():
    print("creating fasta files for multi-sequence-alignment")
    # create_msa_fasta_single_domain_strip()
    print("Running Mafft")
    # run_mafft(folder=MSA_BY_DOMAIN_STRIP)
    print("Loading msas into database")
    # load_msa_single_domain()
    print("Running fasttree")
    # run_fasttree(MSA_BY_DOMAIN_STRIP, TREE_OUT)
    print("Runing notung")
    run_notung(TREE_OUT, NOTUNG_FOLDER)


def load_msa_single_domain():
    """
    load the multisequencealignment result into the db for sequence with single domain only

    it uses the MSA_BY_DOMAIN_STRIP which also includes domains from multidomain sequences
    these improve the quality of the alignment
    """

    name = f"single_domain_strip"
    deleted = uniprot.MultiSequenceAlignment.objects.filter(name=name).delete()
    print(f"Deleted {deleted} with name: {name}")

    to_create = []
    for out_filename in MSA_BY_DOMAIN_STRIP.glob("*out.fasta"):
        family = out_filename.stem.replace("_out", "")
        print(".",end="", flush=True)
        with open(out_filename, "r") as out_file:
            seq_id_to_align = defaultdict(list)
            for line in out_file:
                if line.startswith(">"):
                    seq_id = line[1:].strip()
                else:
                    seq_id_to_align[seq_id].append(line.strip())
            for seq_id, lines in seq_id_to_align.items():
                to_create.append(uniprot.MultiSequenceAlignment(
                    name=name,
                    family=family,
                    seq_id=seq_id,
                    alignment="".join(lines)
                ))
    print(f"\nCreating {len(to_create)} rows")
    uniprot.MultiSequenceAlignment.objects.bulk_create(to_create)
    print("Done")


def create_msa_fasta_single_domain_strip():
    """create fasta files for each domain, containing only the relevant sequence portion

    includes domains from sequences with several domains"""
    superfamily_to_sequences = defaultdict(set)
    for cath_uniprot in cath.SuperfamilyUniprotEntry.objects.all().select_related("seq"):
        superfamily_to_sequences[cath_uniprot.superfamily_id].add(cath_uniprot.seq)

    for superfamily, sequences in superfamily_to_sequences.items():
        fasta_name = f"{superfamily}_in.fasta"
        with open(MSA_BY_DOMAIN_STRIP / fasta_name, "w") as fasta_file:
            for sequence in sequences:
                fasta_file.write(f">{sequence.id}\n")
                for line in wrap(sequence.seq, 80):
                    fasta_file.write(line + "\n")


def run_mafft(folder=MSA_FOLDER):
    """run mafft for every fasta file in folder, overwrites existing files"""
    n_cpu = os.cpu_count()
    commands = []
    for path in folder.iterdir(): 
        if path.stem.endswith("_in"):
            commands.append(f"mafft --anysymbol --thread {n_cpu} --auto --leavegappyregion {path} "
                            f"> {folder}/{path.stem.replace("in", "out")}.fasta "
                            f"2> {folder}/{path.stem}.out")
    with ThreadPool() as pool:
        pool.map(lambda p: subprocess.run(p, shell=True), commands)
  

def run_fasttree(in_folder: Path, out_folder: Path) -> None:
    """
    Runs fasttree on every multi-sequence alignment fasta file in the input folder

    :param input_folder: will run fasttree for every msa fasta file in this folder
    :param output_folder: puts output tree files here. Overwrites existing files
    :return:
    """
    # fasttree runs only on 1 cpu so starting with the slowest and changing map to imap
    # to use cores as they are freed. this avoids having only one process using one core at the end
    commands = []
    for path in sorted(in_folder.iterdir(), key=os.path.getsize, reverse=True): 
        if path.stem.endswith("_out"):
            family = path.stem.replace("_out", "")
            commands.append(f"trimal -in {path} -gappyout | "
                            f"fasttree > {out_folder / family}.tree 2> {out_folder / family}.out")
    with ThreadPool() as pool:
        list(pool.imap(lambda p: subprocess.run(p, shell=True), commands, chunksize=1))


def run_notung(tree_folder: Path, notung_folder: Path):
    """run notun to take into account species informartion in the trees"""
    # this relationship is not one-to-one since some proteins from different
    # species can have the same sequence. this only happens in close species,
    # though, so it is not relevant for the rerooting
    taxid_old2new = {t[0]: t[1] for t in ncbi.db.execute("SELECT * FROM merged").fetchall()}

    seqid_to_taxid = {str(s): str(taxid_old2new.get(t, t)) 
                      for s, t in cath.SuperfamilyUniprotEntry.objects.values_list(
                          "seq_id", "uniprot_entry__taxid")}
    ncbi_tree = ncbi.get_topology(seqid_to_taxid.values())

    # do not overwrite or redo existing output files
    done = set([path.stem.split("_")[0] for path in notung_folder.glob("*rooting.ntglog")])
    tree_files = [path for path in tree_folder.glob("*.tree") if path.stem not in done]

    # tree_files = [path for path in tree_files if path.stem == "1.10.1740.80"]

    # for tree in tree_files:
        # run_notung_tree(tree, ncbi_tree, seqid_to_taxid, notung_folder)
    # return
    with ThreadPool(16) as pool:
        pool.map(lambda tree: run_notung_tree(
            tree, ncbi_tree, seqid_to_taxid, notung_folder), tree_files, chunksize=1)
        # list(pool.imap(lambda tree: run_notung_tree(
            # tree, ncbi_tree, seqid_to_taxid, notung_folder), tree_files, chunksize=1))

def run_notung_tree(tree_file, ncbi_tree, seqid_to_taxid, notung_folder):
    family_id = tree_file.stem
    nw_string = tree_file.read_text()
    nw_tree = etetree(nw_string, sp_naming_function=lambda s: seqid_to_taxid[s])
    taxids = {seqid_to_taxid[n] for n in nw_tree.get_leaf_names()}

    species_name_change, pruned_taxids = write_pruned_ncbi_tree(
            ncbi_tree,
            taxids,
            notung_folder / f"{family_id}.ncbi.nw")

    if len(nw_tree.get_leaves()) < 3 or len(pruned_taxids) < 2:
        with open(notung_folder/ f"{family_id}_empty_rooting.ntglog", "w") as empty_log_file:
            empty_log_file.write("TREE CONTAINS LESS THAN TWO PROTEINS OR SINGLE SPECIES")
        return

    print(family_id)

    midpoint = nw_tree.get_midpoint_outgroup()
    nw_tree.set_outgroup(midpoint)

    for node in nw_tree.iter_leaves():
        taxid = str(seqid_to_taxid[node.name])
        taxid = species_name_change.get(taxid, taxid)
        node.name = f"{taxid}_{node.name}"

    nw_tree.standardize()
    nw_tree.write(outfile=NOTUNG_FOLDER / f"{family_id}.nw", format=2)

    # notung reconcile and rearrange
    batch_filename = NOTUNG_FOLDER / f"{family_id}_batch_rearrange"
    with open(batch_filename, "w") as batch_file:
        batch_file.write("{0}.ncbi.nw\n{0}.nw\n".format(family_id))
    command = (f"java -jar {NOTUNG_JAR} -b {batch_filename} --rearrange "
               f" --costloss 0.01 --prune --speciestag prefix --edgeweights "
               f" name --threshold 0.9 --outputdir {NOTUNG_FOLDER} --silent --log")
    subprocess.run(command, shell=True, stdout=subprocess.DEVNULL)

    batch_root_filename = NOTUNG_FOLDER / f"{family_id}_batch_root"
    with open(batch_root_filename, "w") as batch_file:
        batch_file.write("{0}.ncbi.nw\n{0}.nw.rearrange.0\n".format(family_id))
    command = (f"java -jar {NOTUNG_JAR} -b {batch_root_filename} --root --nolosses "
               f"--speciestag prefix --outputdir {NOTUNG_FOLDER} --silent --log")
    subprocess.run(command, shell=True)  # , stderr=subprocess.DEVNULL)
    return

def write_pruned_ncbi_tree(ncbi_tree, taxids, path):
    """writes the ncbi newick tree here"""
    ncbi_tree = ncbi_tree.copy()
    ncbi_tree.prune(taxids, preserve_branch_length=True)

    # intermediate nodes that appear in taxids will not be recognized by notung 
    # (species when strain are also present). in this case, the function merges all the strains
    # into the same species
    name_change = {}
    for node in ncbi_tree.traverse("preorder"):
        if (not node.is_leaf()) and node.name in taxids:
            for desc in node.get_descendants():
                name_change[desc.name] = node.name
                taxids.discard(desc.name)

    ncbi_tree.prune(taxids, preserve_branch_length=True)
    ncbi_tree.standardize()
    # ncbi_tree.show()
    ncbi_tree.write(outfile=path, format=8)
    return name_change, taxids
