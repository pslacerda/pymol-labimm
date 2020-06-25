import os.path
import platform

import pymol.plugins

if platform.platform() == "MacOS":
    PREF_GUESSES = {
        "LABIMM_VINA": "",
        "LABIMM_OBABEL": "",
        "LABIMM_ADT_PYTHON": "",
        "LABIMM_PREPARE_RECEPTOR": "",
        "LABIMM_PREPARE_FLEXRECEPTOR": "",
    }

elif platform.platform() == "Windows":
    PREF_GUESSES = {
        "LABIMM_VINA": "",
        "LABIMM_OBABEL": "",
        "LABIMM_ADT_PYTHON": "",
        "LABIMM_PREPARE_RECEPTOR": "",
        "LABIMM_PREPARE_FLEXRECEPTOR": "",
    }

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


PLUGIN_DATA_DIR = os.path.expanduser("~/.pymol/labimm")
os.makedirs(PLUGIN_DATA_DIR, exist_ok=True)
