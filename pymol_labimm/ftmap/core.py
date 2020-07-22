import fnmatch
import shutil
import tempfile
from dataclasses import dataclass, field
from glob import glob
from itertools import combinations
from os.path import basename, splitext

import numpy as np
import requests
import scipy.spatial
from cached_property import cached_property
from pymol import CmdException
from pymol import cmd as pm
from pymol import stored

from ..commons import (disable_feedback, fractional_overlap, get_atoms,
                       pairwise, settings, count_molecules)


@dataclass
class Cluster:
    """A consensus site."""

    id: str
    sele: str
    coords: np.ndarray = field(repr=False)

    @cached_property
    def selection(self):
        return self.sele

    @cached_property
    def strength(self):
        return count_molecules(self.sele)

    @cached_property
    def max_dist(self):
        return scipy.spatial.distance_matrix(self.coords, self.coords).max()

    @cached_property
    def polar_count(self):
        stored.resns = []
        pm.iterate_state(0, self.selection, "stored.resns.append(resn)")
        count = 0
        for resn in stored.resns:
            if resn in [
                "ACD",
                "ACN",
                "ACT",
                "ADY",
                "AMN",
                "BDY",
                "BUT",
                "DFO",
                "EOL",
                "PHN",
                "THS",
                "URE",
            ]:
                count += 1
        del stored.resns
        return count

    @cached_property
    def apolar_count(self):
        stored.resns = []
        pm.iterate_state(0, self.selection, "stored.resns.append(resn)")
        count = 0
        for resn in stored.resns:
            if resn in ["BEN", "CHX", "DME", "ETH"]:
                count += 1
        del stored.resns
        return count

    def __hash__(self):
        return hash(self.selection)

    @classmethod
    def collect(cls):
        for id in range(0, 50):
            obj_mask = "consensus.{id:03}.*".format(id=id)
            objs = fnmatch.filter(pm.get_object_list(), obj_mask)
            if len(objs) == 0:
                break
            if len(objs) > 1:
                raise Exception(f"Too much objects found: {', '.join(objs)}")
            obj = objs[0]
            pm.flag("ignore", obj, "clear", True)
            coords = pm.get_coords(obj)
            yield Cluster(id, obj, coords)

    @classmethod
    def collect_atlas(cls):
        """Collect cluster objects from an Atlas session."""
        yield from cls.collect()

    @classmethod
    def collect_ftmap(cls):
        """Collect cluster objects from an FTMap session."""
        for obj in pm.get_object_list():
            if obj.endswith(".pdb"):
                new_name = obj[:-4].replace("crosscluster", "consensus")
                pm.set_name(obj, new_name)
        yield from cls.collect()


@dataclass
class Ensemble:
    """A combination of clusters."""

    clusters: [Cluster]

    @cached_property
    def selection(self):
        return " or ".join(f"({c.selection})" for c in self.clusters)

    @cached_property
    def strength(self):
        return sum(cs.strength for cs in self.clusters)

    @cached_property
    def strength0(self):
        return self.clusters[0].strength

    @cached_property
    def max_dist(self):
        return max(
            scipy.spatial.distance_matrix(cs1.coords, cs2.coords).max()
            for (cs1, cs2) in pairwise(self.clusters, self.clusters)
        )

    @cached_property
    def center_to_center(self):
        cd = {}
        for (cs1, cs2) in pairwise(self.clusters, self.clusters):
            if cs1 == cs2:
                cd[(cs1, cs2)] = 0
            else:
                cd[(cs1, cs2)] = scipy.spatial.distance.euclidean(
                    np.average(cs1.coords, axis=0), np.average(cs2.coords, axis=0)
                )
        return cd

    @cached_property
    def max_center_to_center(self):
        return max(self.center_to_center.values())

    @cached_property
    def continuity(self):
        dist = 0
        for (cs1, cs2) in pairwise(self.clusters, self.clusters):
            if cs1.id == cs2.id:
                continue
            dist = max(
                dist,
                scipy.spatial.distance.euclidean(
                    np.average(cs1.coords, axis=0), np.average(cs2.coords, axis=0)
                ),
            )
        return dist

    @cached_property
    def polar_count(self):
        return sum(c.polar_count for c in self.clusters)

    @cached_property
    def apolar_count(self):
        return sum(c.apolar_count for c in self.clusters)

    def __hash__(self):
        return hash(self.selection)

    @classmethod
    def collect(cls, max_size, cluster_collector):
        clusters = list(cluster_collector())
        for sz in range(1, max_size + 1):
            for cs_comb in combinations(clusters, sz):
                yield cls(cs_comb)

    @classmethod
    def collect_atlas(cls, max_size):
        """Collect possible ensembles from an Atlas session."""
        yield from cls.collect(max_size, Cluster.collect_atlas)

    @classmethod
    def collect_ftmap(cls, max_size):
        """Collect possible ensembles from an FTMap session."""
        yield from cls.collect(max_size, Cluster.collect_ftmap)


class Kozakov2015Ensemble(Ensemble):
    """A combination of clusters as described by Kozakov (2015).
    DOI: 10.1021/acs.jmedchem.5b00586
    """

    @cached_property
    def klass(self):
        return (
            "D"
            if self.is_druggable
            else "Dl"
            if self.is_druggable_large
            else "Ds"
            if self.is_druggable_small
            else "B"
            if self.is_borderline
            else "Bl"
            if self.is_borderline_large
            else "Bs"
            if self.is_borderline_small
            else None
        )

    @cached_property
    def is_druggable(self):
        return (
            self.clusters[0].strength >= 16
            and all(cd < 8 for cd in self.center_to_center.values())
            and self.max_dist >= 10
        )

    @cached_property
    def is_druggable_large(self):
        return (
            self.clusters[0].strength >= 16
            and all(cd >= 8 for cd in self.center_to_center.values())
            and self.max_dist >= 10
        )

    @cached_property
    def is_druggable_small(self):
        return (
            self.clusters[0].strength >= 16
            and all(cd < 8 for cd in self.center_to_center.values())
            and 7 <= self.max_dist < 10
        )

    @cached_property
    def is_borderline(self):
        return (
            13 <= self.clusters[0].strength < 16
            and all(cd < 8 for cd in self.center_to_center.values())
            and self.max_dist >= 10
        )

    @cached_property
    def is_borderline_large(self):
        return (
            13 <= self.clusters[0].strength < 16
            and all(cd >= 8 for cd in self.center_to_center.values())
            and self.max_dist >= 10
        )

    @cached_property
    def is_borderline_small(self):
        return (
            13 <= self.clusters[0].strength < 16
            and all(cd < 8 for cd in self.center_to_center.values())
            and 7 <= self.max_dist < 10
        )


def process_session(ensemble_collector, pattern, group, max_size, base_root=None):
    """Main plugin code."""

    results = {}

    for path in sorted(glob(pattern)):
        if base_root is None:
            root = pm.get_legal_name(splitext(basename(path))[0])
        else:
            root = base_root
        results[root] = [], []
        ensembles, clusters = results[root]

        if group:
            root = f"{group}.{root}"
        else:
            root = root

        with disable_feedback("all", "warnings"):
            with settings(group_auto_mode=1):
                pm.load(path)

        try:
            collected_ensembles = list(ensemble_collector(max_size))
            if len(collected_ensembles) == 0:
                raise Exception()
        except:
            raise CmdException(f"File {path} is invalid.")

        with settings(group_auto_mode=2):
            pm.create(f"{root}.protein", "protein")
            pm.delete("protein")

            i = 0

            for ensemble in collected_ensembles:
                klass = ensemble.klass
                if klass:
                    i += 1
                    pm.hide("sticks", ensemble.selection)
                    pm.show("line", ensemble.selection)
                    pm.util.cbas(ensemble.selection)

                    obj = f"{root}.{klass}.{i:03}"
                    pm.create(obj, ensemble.selection)

                    if hasattr(pm, 'set_property'):
                        pm.set_property("Class", ensemble.klass, obj)
                        pm.set_property("S", ensemble.strength, obj)
                        pm.set_property("S (CS0)", ensemble.clusters[0].strength, obj)
                        pm.set_property("CD", ensemble.max_center_to_center, obj)
                        pm.set_property("MD", ensemble.max_dist, obj)

                    ensemble.selection = obj
                    ensembles.append(ensemble)

            for i, cluster in enumerate(Cluster.collect_atlas()):
                pm.hide("sticks", cluster.selection)
                pm.show("line", cluster.selection)
                pm.util.cbay(cluster.selection)
                obj = f"{root}.CS.{i:03}_{cluster.strength:03}"
                pm.create(obj, cluster.selection)
                pm.delete(cluster.selection)
                cluster.selection = obj
                clusters.append(cluster)

        pm.color("yellow", f"{root}.CS.*")
        pm.color("salmon", f"{root}.B.* or {root}.Bs.* {root}.Bl.*")
        pm.color("red", f"{root}.D.* or {root}.Ds.* {root}.Dl.*")

        pm.hide("lines", f"{root}.*")
        pm.disable(f"{root}.CS")

        pm.show("mesh", f"{root}.B.* or {root}.Bs.* {root}.Bl.*")
        pm.show("mesh", f"{root}.D.* or {root}.Ds.* {root}.Dl.*")

        pm.show("mesh", f"{root}.CS.*")
        pm.hide("nb_spheres", "*label")

        pm.orient(root)
    return results


#
# Plugin commands
#


def load_ftmap(path, group=None, max_cs=3):
    """
    Load a FTMap PDB file and classify hotspot ensembles in accordance to
    Kozakov et al. (2015).
    https://doi.org/10.1021/acs.jmedchem.5b00586

    OPTIONS:
        path    PDB file path, glob or server result id.
        group   optional group name to put objects in.
        max_cs  the maximum number of consensus sites to consider.

    EXAMPLES:
        load_ftmap fftmap.1234.pdb
        load_ftmap fftmap.1234.pdb, GRP, 4
        load_ftmap 79781
    """
    if all(ch in "1234567890" for ch in path):
        fp, temp = tempfile.mkstemp(suffix=".pdb")
        open(fp).close()
        session = requests.Session()
        session.get("https://ftmap.bu.edu/nousername.php")
        with session.get(
            f"https://ftmap.bu.edu/file.php"
            f"?jobid={path}"
            f"&coeffi=0&model=0&filetype=model_file",
            stream=True,
        ) as ret:
            with open(temp, "wb") as fp:
                shutil.copyfileobj(ret.raw, fp)
            if not ret.content:
                raise CmdException(f"Invalid response id {path}")
        return process_session(
            Kozakov2015Ensemble.collect_ftmap,
            temp,
            group,
            int(max_cs),
            root="fftmap" + path,
        )
    return process_session(Kozakov2015Ensemble.collect_ftmap, path, group, int(max_cs))


def load_atlas(path, group=None, max_size=3):
    """
    Load an Atlas PDB file. See `help calculate_ftmap_hotspots`.
    """
    return process_session(
        Kozakov2015Ensemble.collect_atlas, path, group, int(max_size)
    )


def get_fractional_overlap(sel1, sel2, radius=2, state1=1, state2=1, verbose=1):
    """
    Compute the fractional overlap of ligand (A) respective to a hotspot
    ensemble (B).
        FO = Nc/Nt
    Nc is the number of atoms of A in contact with B. Nt is the number of atoms
    of A.

    Hydrogen atoms are ignored.

    If the contact radius is 0 then the VdW radii will be used.

    The states are for select a single state from a multi-state objects.

    OPTIONS:
        sel1    ligand object.
        sel2    hotspot object.
        radius  the radius so two atoms are in contact (default: 2).
        state1  state of sel1.
        state2  state of sel2.

    EXAMPLES:
        get_fractional_overlap ref_lig, ftmap1234.D.003
        get_fractional_overlap ref_lig, ftmap1234.CS.000_016
    """
    atoms1 = get_atoms(sel1, ["coords", "elem", "vdw"], state1)
    atoms2 = get_atoms(sel2, ["coords", "elem", "vdw"], state2)
    fo = fractional_overlap(atoms1, atoms2, int(radius) or None)
    if bool(verbose):
        print(f"  Fractional Overlap = {fo}")
    return fo


def get_fractional_overlap2(sel1, sel2, radius=2, state1=1, state2=1):
    """
    Compute the fractional overlap.

    SEE:
        get_fractional_overlap
    """
    print(sel1 + ' / ' + sel2)
    get_fractional_overlap(sel1, sel2, radius, state1, state2, 1)
    print(sel2 + ' / ' + sel1)
    get_fractional_overlap(sel2, sel1, radius, state2, state1, 1)


def calculate_kozakov2015(*args, **kwargs):
    """
Calculate a hotspot following Kozakov et al (2015).

SIGNATURE:
    calculate_kozakov2015 sel1 ...

EXAMPLES:
    calculate_kozakov2015 *CS.000_*, *CS.002_*
    calculate_kozakov2015 *.000_*, *.001_*

    """
    clusters = []
    for sel in args:
        cluster = Cluster('', sel, pm.get_coords(sel))
        clusters.append(cluster)

    ensemble = Ensemble(clusters)
    print(f"""
{ensemble}
S {ensemble.strength}
S0 {ensemble.strength0}
CD {ensemble.max_center_to_center}
MD {ensemble.max_dist}
    """)


def nearby_aminoacids_similarity(sel1, sel2, radius=4, state1=1, state2=1):
    """
    Compute the similarity of the aminoacids nearby two selections.

    SIGNATURE:
        nearby_aminoacids_similarity sel1, sel2, [radius=4,] [state1=1,] [state2=1]

    OPTIONS:
        sel1    Selection of object 1.
        sel2    Selection of object 2.
        radius  Radius to look for nearby aminoacids.
        state1  State of object 1.
        state2  State of object 2.

    EXAMPLES:
        nearby_aminoacids_similarity *CS.000_*, *CS.002_*
        nearby_aminoacids_similarity *.000_*, *.001_*
    """
    resis1 = get_atoms(sel1, ["resi"], state1)
    resis2 = get_atoms(sel2, ["resi"], state2)




def init_plugin_cli():
    pm.extend(load_ftmap)
    pm.extend(load_atlas)
    pm.extend(get_fractional_overlap)
    pm.extend(get_fractional_overlap2)
    pm.extend(calculate_kozakov2015)
    pm.extend(nearby_aminoacids_similarity)
