import subprocess
import sys


def pip(args):
    process = subprocess.Popen(
        [sys.executable, '-m', 'pip', '--disable-pip-version-check'] + args,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = process.communicate()
    return process, out, err


try:
    from pymol_labimm import __init_plugin__
except ImportError:
    proc, out, err = pip(['install', 'pymol-labimm'])
    if out:
        print(out)
    if err:
        print(err)
    from pymol_labimm import __init_plugin__


__init_plugin__()
