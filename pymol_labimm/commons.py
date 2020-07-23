import hashlib
import random
from contextlib import contextmanager

import numpy as np
import pandas as pd
import scipy.spatial
from pymol import cmd as pm
from pymol import stored

#
# Generic
#


def pairwise(iter1, iter2):
    """Pairwise element iteration."""
    for item1 in iter1:
        for item2 in iter2:
            yield (item1, item2)


#
# PyMOL specific
#


def new_object_name(prefix="", suffix=""):
    """Return an unused object name."""
    while True:
        seed = str(random.random())
        obj = hashlib.sha224(seed.encode()).hexdigest()
        obj = prefix + obj + suffix
        if obj not in pm.get_object_list():
            return obj


@contextmanager
def disable_feedback(what, level):
    """Disable feedback."""
    pm.feedback("disable", what, level)
    yield
    pm.feedback("enable", what, level)


@contextmanager
def settings(**kwargs):
    """Set options temporarily."""
    orig = {}
    for key, value in kwargs.items():
        orig[key] = pm.get(key)
        pm.set(key, value)
    yield
    for key, value in orig.items():
        pm.set(key, value)


#
# Atomic helpers
#


def contact_matrix(a, b, radius=None):
    """
    Compute the contact matrix between A and B.
    When radius is None the contact is defined by:
        Distance(a, b) < VDW(a) + VDW(b)
    When radius is a number the contact is defined by:
        Distance(a, b) < radius
    """
    d = scipy.spatial.distance_matrix(a[["x", "y", "z"]], b[["x", "y", "z"]])
    if radius is None:
        return ((d - b["vdw"].values).T - a["vdw"].values).T <= 0
    else:
        return (d - radius) <= 0


def fractional_overlap(a, b, radius=None):
    """Compute the fractional overlap of a respective to b."""
    a = a[a["elem"] != "H"]
    b = b[b["elem"] != "H"]
    contacts = np.any(contact_matrix(a, b, radius), axis=1)
    num_contacts = np.sum(contacts)
    total_atoms = len(a)
    fo = num_contacts / total_atoms
    return fo


def get_atoms(sel, attrs, state=1):
    """Get the atoms and attributes of a selection."""
    coords = None
    if "coords" in attrs:
        coords = pm.get_coords(sel, state)
        attrs.remove("coords")
    atoms = pd.DataFrame(coords, columns=["x", "y", "z"])

    if attrs:
        fields_str = ", ".join(attrs)
        stored.atoms = []
        pm.iterate_state(state, sel, f"stored.atoms.append(({fields_str}))")
        atoms = pd.concat([atoms, pd.DataFrame(stored.atoms, columns=attrs)], axis=1)
        del stored.atoms
    return atoms


def count_molecules(selection="all"):
    """
    By Thomas Holder.
    """
    tmpsele = pm.get_unused_name("_tmp")
    count = 0
    if pm.select(tmpsele, selection):
        count += 1
        while pm.select(tmpsele, f"{tmpsele} &! bm. first {tmpsele}"):
            count += 1
    pm.delete(tmpsele)
    return count
