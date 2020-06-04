import requests
from lxml import etree
from pymol import cmd as pm


def find_similar_chain_ids(chain_id, threshold):
    """Fetch structure similar chain ids from RCSB PDB.
    chain_id - reference chain id
    threshold - similarity threshold
    """
    resp = requests.get(
        f"https://www.rcsb.org/pdb/rest/sequenceCluster"
        f"?cluster={threshold}"
        f"&structureId={chain_id}"
    )
    tree = etree.fromstring(resp.content)
    for pdb_chain in tree.xpath("//pdbChain"):
        sim_chain_id = pdb_chain.get("name")
        rank = pdb_chain.get("rank")
        sim_pdb, sim_chain = sim_chain_id.split(".")
        yield sim_pdb, sim_chain_id, int(rank)


def get_resolution(pdb_id):
    """
    Get the resolution for a PDB id, or None in case it isn't of an X-ray experiment.
    """
    ret = requests.get(
        f"https://www.rcsb.org/pdb/rest/getEntityInfo?structureId={pdb_id}"
    )
    tree = etree.fromstring(ret.content)
    if tree.xpath("//Method")[0].get("name") == "xray":
        resol = tree.xpath("//PDB")[0].get("resolution")
        if not resol:
            return float("nan")
        return float(resol)
    else:
        return None


@pm.extend
def fetch_similar(
    chain_id,
    similarity=95,
    ligand=None,
    dist=5,
    compounds="organic or inorganic",  # pep_compounds=None,
    max_resolution=None,
    max_structures=50,
):
    """
Fetch sequence similar structures from RCSB PDB and optionally keep only
apo structures. Apo are evaluated respective to a choosen ligand on the
reference chain.
OPTIONS:
    chain_id        Reference structure chain id.
    similarity      Sequence similarity threshold (one of the available
                    from RCSB PDB).
    ligand          Reference ligand PDB id.
    dist            Distance cut-off around reference ligand for apo
                    evaluation.
    compounds       Selection of atoms that should be considered ligands
                    upon apo computation. Only used when ligand is given.
    max_resolution  Fetch only X-ray structures with up to such
                    resolution.
    max_structures  Fetch at most n structures. 0 for all structures.
EXAMPLES:
    fetch_similar 2XY9.A, 100
    fetch_similar 2XY9.A, 95, 3ES, 3, organic
    fetch_similar 6Y2F.A, max_structures=0
    """

    max_structures = int(max_structures)
    obj = f"{chain_id}.0"
    pm.fetch(chain_id, obj)

    sims = []
    similars = find_similar_chain_ids(chain_id, similarity)
    cont = 0
    for sim_pdb, sim_chain_id, rank in similars:

        if max_structures != 0 and cont >= max_structures:
            break

        if sim_chain_id.upper() == chain_id.upper():
            continue

        sim_obj = f"{sim_chain_id}.{rank}"
        pm.fetch(sim_chain_id, sim_obj, **{"async": 0})
        pm.align(sim_obj, obj)

        resol = None
        if max_resolution:
            resol = get_resolution(sim_pdb)
            if not resol or resol > max_resolution:
                pm.delete(sim_obj)
                continue
        if ligand:
            model = pm.get_model(
                f"({sim_obj} and ({compounds}))"
                f" within {dist} of"
                f"({obj} and (resn {ligand}))"
            )
            if len(model.atom) != 0:
                pm.delete(sim_obj)
                continue

            # for lig_obj in pm.get_object_list(sim_obj):
            #     if len(pm.get_fastastr(obj)) < 15:
            #         pm.get_model(
            #             f'%{lig_obj} within {dist} of '
            #         )

        cont += 1
        sims.append((sim_obj, sim_chain_id, sim_pdb, rank, resol))
    return sims
