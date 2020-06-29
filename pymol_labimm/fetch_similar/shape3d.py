import requests
from pymol import cmd as pm
from dataclasses import dataclass
from .blast import get_resolution
from typing import Optional


@dataclass
class SimilarStructure:
    pdb_id: str
    score: float
    resolution: Optional[float]


@pm.extend
def fetch_similar_shape3d(pdb_id: str,
                            min_similarity: float = 0,
                            max_resolution: Optional[float] = None,
                            ligand: Optional[str] = None,
                            dist: float = 5.0,
                            compounds: str = "organic or inorganic",
                            prosthetic_groups: str = "HEM FAD NAP NDP ADP FMN",
                            ) -> [SimilarStructure]:
    ret = requests.post(
        url="https://search.rcsb.org/rcsbsearch/v1/query",
        headers={"Content-Type": "application/json"},
        json={
              "query": {
                "type": "terminal",
                "service": "structure",
                "parameters": {
                  "value": {
                    "entry_id": pdb_id,
                    "assembly_id": "1"
                  },
                  "operator": "strict_shape_match"
                },
                "node_id": 1
              },
              "return_type": "entry"
        }
    )

    similars = []
    for i, result in enumerate(ret.json()["result_set"]):
        sim = SimilarStructure(
            pdb_id=result["identifier"],
            score=result["score"],
            resolution=None,
        )

        # Chek the similarity threshold
        if sim.score < min_similarity:
            continue

        # Check the resolution
        if i >= 1 and max_resolution:
            resol = get_resolution(sim.pdb_id)
            if resol and resol > max_resolution:
                continue
            sim.resolution = resol

        # Fetch the structure
        pm.fetch(sim.pdb_id)

        if i >= 1:
            first_pdb_id = similars[0].pdb_id

            # Align
            pm.align(sim.pdb_id, first_pdb_id)

            # Check nearby non-prosthetic ligands
            # Apo detection
            if ligand:
                model = pm.get_model(
                    f"({sim.pdb_id} and ({compounds}))"
                    f" within {dist} of"
                    f"({first_pdb_id} and (resn {ligand}))"
                )
                resns = set(a.resn for a in model.atom)
                is_apo = True
                for resn in resns:
                    if resn not in prosthetic_groups.split():
                        is_apo = False
                        break
                if not is_apo:
                    pm.delete(sim.pdb_id)
                    continue

        # Set the score property
        if hasattr(pm, 'set_property'):
            pm.set_property("shape3d_score", sim.score, sim.pdb_id)

        # Go on
        similars.append(sim)
    return similars
