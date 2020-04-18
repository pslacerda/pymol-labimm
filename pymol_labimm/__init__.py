__version__ = "0.1.0"

import pymol.gui

from .ftmap.core import init_plugin_cli as ftmap_init_plugin_cli
from .ftmap.gui import init_plugin_gui as ftmap_init_plugin_gui
from .vina import init_plugin as vina_init_plugin
from .prefs import guess_prefs


def __init_plugin__(app=None):
    guess_prefs()

    # FTMap
    ftmap_init_plugin_cli()

    window = pymol.gui.get_qtwindow()
    menu_bar = window.menuBar()
    labimm_menu = menu_bar.addMenu("LaBiMM")

    ftmap_init_plugin_gui(labimm_menu)

    # Vina
    vina_init_plugin(labimm_menu)
