from pymol import cmd as pm

from scipy.spatial import distance_matrix
import numpy as np
from matplotlib import pyplot as plt


@pm.extend
def dc(sel1, sel2, state1=1, state2=1, dist=1.25, verbose=1):
    """Compute the Density Correlation according to:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3264775/

USAGE:
    dc(sel1, sel2, state1=1, state2=1, dist=1.25, verbose=1)

sel1 and sel2 are respectively the molecule and hotspot; state1 and state2 are
the optional corresponding states (default to first state both). The threshold
distance can be changed with dist (default to 1.25).

verbose is the standard boolean API option to define verbosity.
    """
    xyz1 = pm.get_coords(f'({sel1}) and not elem H', state1)
    xyz2 = pm.get_coords(f'({sel2}) and not elem H', state2)
    dc_ = (distance_matrix(xyz1, xyz2) < float(dist)).any(1).sum()
    if bool(verbose):
        print(f'DC = {dc_:.2f}')
    return dc_

@pm.extend
def dce(sel1, sel2, state1=1, state2=1, dist=1.25, verbose=1):
    """Compute the Density Correlation Efficiency according to:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3264775/

USAGE:
    dce(sel1, sel2, state1=1, state2=1, dist=1.25, verbose=1)

All parameters are the same for dc. The reference selection wich atoms are
counted is sel1.
    """
    dc_ = dc(sel1, sel2, state1, state2, dist, verbose=False)
    dce_ = dc_/pm.count_atoms(f'({sel1}) and not elem H')
    if bool(verbose):
        print(f'DCE = {dce_:.2f}')
    return dce_


@pm.extend
def plot_dce_le(conf):
    with open(conf) as conf_file:
        hs = conf_file.readline().strip()
        ligs = []
        for line in conf_file:
            obj, pki = line.split()
            pki = float(pki)
            ha = pm.count_atoms(f'({obj}) and not elem H')
            le = pki / ha
            dce_ = dce(hs, obj, verbose=0)
            ligs.append((obj, pki, ha, le, dce_))
        ligs.sort(key=lambda l: l[2])
        _, _, _, ref_le, ref_dce = ligs.pop(0)
        
        x = []
        y = []
        labels = []

        for obj, _, _, le, dce_ in ligs:
            x.append(le/ref_le)
            y.append(dce_/ref_dce)
            labels.append(obj)
        
        m, b = np.polyfit(x, y, 1)
        fig, ax = plt.subplots()
        plt.scatter(x, y)
        plt.plot(x, m*np.array(x)+b)

        for index in range(len(x)):
            ax.text(x[index], y[index], labels[index], size=12)
        plt.show()
