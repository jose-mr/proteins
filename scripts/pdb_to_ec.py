
import csv
import datetime
from collections import defaultdict

from wpdb.models import Entry as PdbEntry
from pseudoenzymes.settings import OUT_FOLDER


def run():
    print("oi")
    
    date_2021 = datetime.date(2021, 1, 1)
    all_pdbs = PdbEntry.objects.all().prefetch_related("uniprot_entries__ec_entries")
    all_pdbs_ecs = all_pdbs.values_list("pdb_id", "uniprot_entries__ec_entries").order_by("date")
    all_ecs = defaultdict(set)
    for pdb_id, ec in all_pdbs_ecs:
        if ec is not None:
            all_ecs[ec].add(pdb_id)

    previous_pdbs = all_pdbs.filter(date__lt=date_2021)
    previous_pdbs_ec = previous_pdbs.values_list("pdb_id", "uniprot_entries__ec_entries").order_by("date")
    previous_ecs = defaultdict(set)
    for pdb_id, ec in previous_pdbs_ec:
        if ec is not None:
            previous_ecs[ec].add(pdb_id)

    recent_pdbs = all_pdbs.filter(date__gte=date_2021)

    recent_pdbs_ec = recent_pdbs.values_list("pdb_id", "uniprot_entries__ec_entries").order_by("date")
    recent_ecs = defaultdict(set)
    recent_only_ecs = defaultdict(set)

    for pdb_id, ec in recent_pdbs_ec:
        if ec is not None:
            recent_ecs[ec].add(pdb_id)
            if ec not in previous_ecs:
                recent_only_ecs[ec].add(pdb_id)

    print("ALL PDBs: ", all_pdbs.count())
    print("PDBS before 2021: ",  previous_pdbs.count())
    print("PDBS on and after 2021: ",  recent_pdbs.count())

    print("All ECs represented in PDB: ", len(previous_ecs.keys()|recent_ecs.keys()))
    print("ECs represented in PDB before 2021: ", len(previous_ecs))
    print("ECs represented in PDB on an after 2021: ", len(recent_ecs))
    print("ECs represented in PDB exclusively on an after 2021 (new ecs): ", len(recent_only_ecs))

    pdb_id_to_entry = {e.pdb_id: e for e in PdbEntry.objects.all()}

    fieldnames = ["PDB_ID", "NAME", "YEAR", "DATE", "EC", "UniProt IDS", "CATH_IDS"]
    with open(OUT_FOLDER/"pdb2ec_recent_only.csv", "w") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(fieldnames)
        for ec, pdbs in recent_only_ecs.items():
            for pdb_id in pdbs:
                pdb = pdb_id_to_entry[pdb_id]
                uniprot_ids = ";".join(list(pdb.uniprot_entries.values_list("ac", flat=True)))
                cath_superfamilies = set()
                for uniprot in pdb.uniprot_entries.all():
                    cath_superfamilies.update(uniprot.cath_superfamilies.values_list("number", flat=True))
                cath_superfamilies = ";".join(cath_superfamilies)
                csv_writer.writerow([pdb_id, pdb.title, pdb.date.year, pdb.date, ec, uniprot_ids, cath_superfamilies])

    # TODO check why too few associations
    with open(OUT_FOLDER/"pdb2ec_all.csv", "w") as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(fieldnames)
        for ec, pdbs in (all_ecs).items():
            for pdb_id in pdbs:
                pdb = pdb_id_to_entry[pdb_id]
                uniprot_ids = ";".join(list(pdb.uniprot_entries.values_list("ac", flat=True)))
                cath_superfamilies = set()
                for uniprot in pdb.uniprot_entries.all():
                    cath_superfamilies.update(uniprot.cath_superfamilies.values_list("number", flat=True))
                cath_superfamilies = ";".join(cath_superfamilies)
                csv_writer.writerow([pdb_id, pdb.title, pdb.date.year, pdb.date, ec, uniprot_ids, cath_superfamilies])

