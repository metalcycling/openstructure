from ost import settings
from ost import bindings
import subprocess
import re

def GetHHblitsVersionString():
    version_string = None

    # raises if hhblits binary is not in path
    try:
        hhblits_bin = settings.Locate('hhblits')
    except:
        raise RuntimeError('Tried to determine HHblits version string which '\
                           'requires the hhblits binary to be in your path. '\
                           'Couldnt find it...')

    # run hhblits to determine version
    proc = subprocess.run([hhblits_bin, '-h'], stdout=subprocess.PIPE)

    # regular expression in the form
    # HHblits, whatever except newline, x.y.z
    # where x.y.z are the version numerals
    version_line = re.search(r'HHblits[^\n]+\d+\.\d+\.\d+', proc.stdout.decode())
    if version_line is not None:
        version = re.search(r'\d+\.\d+\.\d+', version_line.group())
        if version is not None:
            version_string = version.group()

    return version_string

hhblits_version_string = GetHHblitsVersionString()

if hhblits_version_string is None:
    raise RuntimeError('Could not determine HHblits version. Please '\
                       'import the hhblits2 or hhblits3 binding explicitely.')

hhblits_version = int(hhblits_version_string.split('.')[0])

if hhblits_version == 2:
    from ost.bindings.hhblits2 import *
elif hhblits_version == 3:
    from ost.bindings.hhblits3 import *
else:
    raise RuntimeError('Determined HHblits version to be %i. OpenStructure '\
                       'only supports 2 or 3. If you think the detected '\
                       'version is wrong and you have version 2 or 3, '\
                       'import the hhblits2 or hhblits3 binding explicitely.'\
                       %(hhblits_version))

