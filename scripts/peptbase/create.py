#!/usr/bin/env python3
"""
This scripts scraps the Binding MOAD looking for protein-peptide complexes. It goes further and also look for unbound
95% similar structures.
"""
import pandas as pd
import requests_cache
from pymol import cmd as pm

requests_cache.install_cache("requests_cache")

from pymol_labimm.fetch_similar.blast import find_similar_chain_ids, get_resolution


def convert_unit(value, unit):
    mult = None
    if unit == "uM":
        mult = 1000
    if unit == "mM":
        mult = 1e6
    if unit == "pM":
        mult = 0.0001
    if unit == "nM":
        mult = 1
    if mult:
        return value * mult


#
# MAIN FUNCTIONS
#

AMINOACIDS = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL PYL SEC"


def parse_binding_moad(moad_csv_file):
    entries = []

    for i, line in enumerate(open(moad_csv_file)):
        parts = line.split(",")

        # EC number (ignore)
        if parts[0]:
            continue

        # PDB
        if parts[2]:
            pdb_id = parts[2].upper()

        # There is a ligand record
        if parts[3]:

            # Ignorable molecule
            if parts[4] == "invalid":
                continue

            # Valid ligands
            elif parts[4] == "valid":
                ligand_id, chain, resid = parts[3].split(":")
                ligand_smiles = parts[9]
                entry = {
                    "type": "ligand",
                    "pdb_id": pdb_id,
                    "ligand_id": ligand_id,
                    "ligand_smiles": ligand_smiles,
                    "ligand_resid": resid,
                    "ligand_chain": chain,
                }

                # Check for peptid ligands.
                if " " in ligand_id:
                    if not all(l in AMINOACIDS.split() for l in ligand_id.split()):
                        continue

                    # Affinity data found.
                    if parts[5]:
                        affinity_measure = parts[5]
                        affinity_value = float(parts[7])
                        affinity_unit = parts[8]

                        affinity_value = convert_unit(affinity_value, affinity_unit)

                        # Ligand with incorrect unit (eg. fM or M)
                        if not affinity_value:
                            continue

                        entry.update(
                            {
                                "affinity_measure": affinity_measure,
                                "affinity_value": affinity_value,
                            }
                        )

                        entries.append(
                            {**entry,}
                        )

    # Structures with no valid ligand must be removed
    entries = pd.DataFrame(entries)
    counts = entries[entries.type == "ligand"].groupby("pdb_id").count().type
    entries = entries[entries.pdb_id.isin(counts[counts > 0].index)]
    entries = entries.groupby(["pdb_id"]).first().reset_index()


    #
    # Find apo/holo pairs
    #
    new_entries = []
    for i, entry in entries.iterrows():

        # if entry.pdb_id != "1OBX":
        #     continue
        # breakpoint()

        pm.reinitialize()

        pm.fetch(entry.pdb_id)
        if entry.pdb_id not in pm.get_object_list():
            continue

        # Store the fasta string length
        ref_len = len(pm.get_fastastr(entry.pdb_id))

        # Look for similar chains using the Blast clusterization algorithm
        similars = find_similar_chain_ids(entry.pdb_id.upper() + "_A", 95)
        for sim_pdb, sim_chain in similars:

            # if sim_pdb != "1R6J":
            #     continue
            # breakpoint()

            if sim_pdb == entry.pdb_id.upper():
                continue

            pm.fetch(sim_pdb)
            if sim_pdb not in pm.get_object_list():
                continue

            # Remove structures with different number of aminoacids. I'm
            # considering that peptides impact less than 50 characters on
            # the fasta string.
            cur_len = len(pm.get_fastastr(sim_pdb))
            if abs(ref_len - cur_len) >= 50:
                pm.delete(sim_pdb)
                continue

            # Remove alternative conformations, which may be useful for docking.
            pm.remove('not (alt "" or alt A)')

            # Check the resolution
            resol = get_resolution(sim_pdb)
            if resol is None or resol > 2.5:
                pm.delete(sim_pdb)
                continue

            # Align the chains
            pm.align(sim_pdb, entry.pdb_id)

            # Check nearby ligands
            is_apo = True
            model = pm.get_model(
                f"({sim_pdb} and (organic)) within 5 of ({entry.pdb_id} and (bysegi resid {entry.ligand_resid} and chain {entry.ligand_chain}))"
            )
            resns = set(a.resn for a in model.atom)
            for resn in resns:
                # ignore if is a prosthetic group
                prosthetic_groups = "HEM FAD NAP NDP ADP FMN"
                if resn not in prosthetic_groups.split():
                    pm.delete(sim_pdb)
                    is_apo = False
                    break
            if not is_apo:
                continue

            # Check nearby peptides
            model = pm.get_model(
                f"({sim_pdb} and polymer) within 5 of ({entry.pdb_id} and (bysegi resid {entry.ligand_resid} and chain {entry.ligand_chain}))"
            )
            is_apo = True
            for at in model.atom:
                fasta = pm.get_fastastr(
                    f"bs. {sim_pdb} and chain {at.chain} and resid {at.resi}"
                )
                for j, fast in enumerate(fasta.split("\n>")):
                    # Peptides have a fasta string smaller than 25 characters
                    if len(fast) <= 25:
                        pm.delete(sim_pdb)
                        is_apo = False
                        break
                if not is_apo:
                    break
            if not is_apo:
                continue

            # Found a pair apo/holo
            new_entries.append({**entry, "apo95_pdb": sim_pdb})
            pm.save(f"{entry.pdb_id}.pdb", entry.pdb_id)
            pm.save(f"{sim_pdb}.pdb", sim_pdb)
            break
        else:
            continue
    return pd.DataFrame(new_entries)


#
# Entrada e Saída de Arquivos
#

PEPTBASE = parse_binding_moad("nr.csv")
PEPTBASE.to_csv("peptbase.csv")
