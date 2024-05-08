from textwrap import wrap
import subprocess
from multiprocessing.pool import ThreadPool
import os
from collections import defaultdict
from django.db import connection
from django.db.models import Count
from Bio import Align

import cath.models as cath
import uniprot.models as uniprot
from pseudoenzymes.settings import MSA_FOLDER, MSA_BY_DOMAIN_ARCHITECTURE, MSA_BY_DOMAIN_STRIP

def run():
    print("creating fasta files for multi-sequence-alignment")
    # create_msa_fasta_same_domain_architecture()
    # run_mafft(folder=MSA_BY_DOMAIN_ARCHITECTURE)

    create_msa_fasta_single_domain_strip()
    run_mafft(folder=MSA_BY_DOMAIN_STRIP)
    load_msa_single_domain()


def load_msa_single_domain():
    """
    load the multisequencealignment result into the db for sequence with single domain only

    it uses the MSA_BY_DOMAIN_STRIP which also includes domains from multidomain sequences
    these improve the quality of the alignment but are not saved in the db
    """

    name = f"single_domain_strip"

    for out_filename in MSA_BY_DOMAIN_STRIP.glob("*out.fasta"):
        family = out_filename.stem.replace("_out", "")
        deleted = uniprot.MultiSequenceAlignment.objects.filter(family=family, name=name).delete()
        print("Deleted", deleted)
        alignment = Align.read(out_filename, "fasta")
        print(f"Reading {len(alignment)} aligned sequences")
        to_create = []
        with open(out_filename, "r") as out_file:
            seq_id_to_align = defaultdict(list)
            for line in out_file:
                if line.startswith(">"):
                    seq_id = line.split()[1]
                else:
                    seq_id_to_align[seq_id].append(line.strip())
            for seq_id, lines in seq_id_to_align.items():
                to_create.append(uniprot.MultiSequenceAlignment(
                    name=name,
                    family=family,
                    seq_id=seq_id,
                    alignment="".join(lines)
                ))
        print(f"Creating {len(to_create)} msa rows")
        uniprot.MultiSequenceAlignment.objects.bulk_create(to_create)


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
                fasta_file.write(f"> {sequence.id}\n")
                for line in wrap(sequence.seq, 80):
                    fasta_file.write(line + "\n")

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


def run_mafft(folder=MSA_FOLDER):
    """run mafft for every fasta file in folder, overwrites existing files"""
    fasta_filenames = [path for path in folder.iterdir() if path.stem.endswith("_in")]
    n_cpu = 16
    commands = []
    for filename in fasta_filenames:
        command = (f"mafft --anysymbol --thread {n_cpu} --auto --leavegappyregion {filename}"
                   f"> {folder}/{filename.stem.replace("in", "out")}.fasta"
                   f"2> {folder}/{filename.stem}.out")
        commands.append(command)
    with ThreadPool(n_cpu) as pool:
        pool.map(lambda p: subprocess.run(p, shell=True), commands)

