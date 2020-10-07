import subprocess
import sys
try:
    import conda.cli.python_api as conda
except:
    pass

def pip(args):
    process = subprocess.Popen(
        [sys.executable, "-m", "pip", "--disable-pip-version-check"] + args,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = process.communicate()
    return out, err, process.returncode


try:
    from pymol_labimm import init_plugin
except:
    try:
        out, err, rv = conda.run_command("install", "-c", "conda-forge", "-y", "-q", "openbabel", "plip")
        if out:
            print(out)
        if err:
            print(err)
    except:
        pass
    out, err, rv = pip(["install", "--upgrade", "pymol-labimm"])
    if out:
        print(out)
    if err:
        print(err)

    from pymol_labimm import init_plugin


init_plugin()
