import subprocess, sys


def pip(args):
    process = subprocess.Popen(
        [sys.executable, '-m', 'pip', '--disable-pip-version-check'] + args,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = process.communicate()
    return process, out, err


proc, out, err = pip(['install', 'git+https://github.com/pslacerda/labimm-pymol'])

if out:
    print(out)
if err:
    print(err)
