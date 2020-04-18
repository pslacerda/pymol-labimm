import platform

import pymol.plugins


if platform.platform() == "MacOS":
    PREFS_GUESSES = {}

elif platform.platform() == "Windows":
    PREFS_GUESSES = {}

else:  # platform == 'Linux'
    PREFS_GUESSES = {
        "LABIMM_VINA": "/usr/bin/vina",
        "LABIMM_OBABEL": "/usr/bin/obabel",
        "LABIMM_ADT_PYTHON": "/usr/bin/python2.7",
        "LABIMM_PREPARE_RECEPTOR": "/usr/lib/python2.7/dist-packages/AutoDockTools/Utilities24/prepare_receptor4.py",
        "LABIMM_PREPARE_FLEXRECEPTOR": "/usr/lib/python2.7/dist-packages/AutoDockTools/Utilities24/prepare_flexreceptor4.py",
    }


def guess_prefs():
    for pref in PREFS_GUESSES:
        if not pymol.plugins.pref_get(pref):
            pymol.plugins.pref_set(pref, PREFS_GUESSES[pref])
