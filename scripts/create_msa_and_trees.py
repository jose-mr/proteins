from textwrap import wrap
import subprocess
from multiprocessing.pool import ThreadPool
import os
from collections import defaultdict

import cath.models as cath
import uniprot.models as uniprot
from pseudoenzymes.settings import MSA_FOLDER, MSA_BY_DOMAIN_ARCHITECTURE, MSA_BY_DOMAIN_STRIP

def run():
    print("creating fasta files for multi-sequence-alignment")
    # create_msa_fa# sta_files()
    # create_msa_fasta_same_domain_architecture()
    # run_mafft(folder=MSA_BY_DOMAIN_ARCHITECTURE)

    # create_msa_fasta_single_domain_strip()
    run_mafft(folder=MSA_BY_DOMAIN_STRIP)


def create_msa_fasta_single_domain_strip():
    """create fasta files for each domain, containing only the relevant sequence portion"""
    cath_uniprots = cath.SuperfamilyUniprotEntry.objects.all()\
            .order_by("uniprot_entry_id", "start")\
            .select_related("uniprot_entry")

    domain_to_seqs = defaultdict(list)
    for cath_uniprot in cath_uniprots:
        uniprot = cath_uniprot.uniprot_entry
        superfamily = cath_uniprot.superfamily_id
        domain_to_seqs[superfamily].append(
                (uniprot, int(cath_uniprot.start), int(cath_uniprot.end)) )

    for superfamily, sequences_start in domain_to_seqs.items():
        if len(sequences_start) >= 10:
            fasta_name = f"{superfamily}_in.fasta"
            print(fasta_name)
            with open(MSA_BY_DOMAIN_STRIP / fasta_name, "w") as fasta_file:
                lines = []
                for sequence, start, end in sequences_start:
                    print(sequence,start, end)
                    lines.append(f"> {sequence.ac} {start} {end}")
                    stripped_seq = sequence.seq[start-1:end]
                    lines.extend(wrap(stripped_seq, 80))
                fasta_file.writelines("\n".join(lines))

def create_msa_fasta_same_domain_architecture():
    """create fasta files with sequences that have the same domain architecture"""
    cath_uniprots = cath.SuperfamilyUniprotEntry.objects.all()\
            .order_by("uniprot_entry_id", "start")\
            .select_related("uniprot_entry")

    sequence_to_domain_architecture = defaultdict(list)
    for cath_uniprot in cath_uniprots:
        uniprot = cath_uniprot.uniprot_entry
        superfamily = cath_uniprot.superfamily_id
        sequence_to_domain_architecture[uniprot].append(superfamily)

    domain_architecture_to_sequences = defaultdict(list)
    for sequence, domain_architecture in sequence_to_domain_architecture.items():
        domain_architecture_to_sequences[tuple(domain_architecture)].append(sequence)

    for domain_architecture, sequences in domain_architecture_to_sequences.items():
        if len(sequences) >= 10:
            fasta_name = f"{'_'.join(domain_architecture)}_in.fasta"
            print(fasta_name)
            with open(MSA_BY_DOMAIN_ARCHITECTURE / fasta_name, "w") as fasta_file:
                lines = []
                for sequence in sequences:
                    lines.append(f"> {sequence.ac}")
                    lines.extend(wrap(sequence.seq, 80))
                fasta_file.writelines("\n".join(lines))

def create_msa_fasta_files(min_number_seqs=10):
    # create fasta files with all the sequences belonging to each cath superfamily
    superfamilies = cath.Superfamily.objects.superfamilies().prefetch_related("uniprot_entries")

    for superfamily in superfamilies:
        if (uniprot_entries := superfamily.uniprot_entries.all()).count() >= min_number_seqs:
            # with open(MSA_FOLDER / f"{str(superfamily).replace(".", "_")}.fasta", "w") as fasta_file:
            with open(MSA_FOLDER / f"{superfamily}_in.fasta", "w") as fasta_file:
                lines = []
                for uniprot_entry in uniprot_entries:
                    lines.append(f"> {uniprot_entry.ac}")
                    lines.extend(wrap(uniprot_entry.seq, 80))
                fasta_file.writelines("\n".join(lines))




def run_mafft(delete_old=False, folder=MSA_FOLDER):
    """
    Runs mafft for every fasta file in folder
    """
    fasta_filenames = [path for path in folder.iterdir() if path.stem.endswith("_in")]

    n_cpu = 16
    commands = []
    for filename in fasta_filenames:
        # command = "{} --thread {} --auto --leavegappyregion {} --anysymbol {} > {}/{}.out.fasta".format(
        #     MAFFT_EXE, 1 if turbo else 1, "--memsave " if turbo else "", filename, output_folder, uniprot_id)
        command = f"mafft --thread {n_cpu} --auto --leavegappyregion {filename} > {folder}/{filename.stem.replace("in", "out")}.fasta 2> {folder}/{filename.stem}.out"
        commands.append(command)
        # subprocess.run(command, shell=True)

    with ThreadPool(n_cpu) as pool:
        pool.map(lambda p: subprocess.run(p, shell=True), commands)

