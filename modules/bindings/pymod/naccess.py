'''
Naccess module

Author: Florian Kiefer, Gerardo Tauriello (cleanup/speedup)

This module is for calculating surface areas
from OpenStructure using the external program naccess

How To Use This Module:
 1. Import it (e.g. as "from ost.bindings import naccess")
 2. Use it (e.g. as "sasa = naccess.CalculateSurfaceArea(entity)")

Requirement:
 - naccess installed
'''

import tempfile
import subprocess
import os
import re
from ost import io, mol, settings, LogWarning

def _GetExecutable(naccess_exe):
  """
  Method to check if naccess executable is present

  :param naccess:   Explicit path to naccess executable
  :returns:         Path to the executable
  :exception:       FileNotFound if executable is not found
  """
  return settings.Locate('naccess', explicit_file_name=naccess_exe)

def _CheckNaccessRoot(naccess_root):
  """
  :return: True, if given directory contains "accall" binary and files
           "vdw.radii" and "standard.data".
  :param naccess_root: Path to naccess folder to check.
  """
  accall_exe = os.path.join(naccess_root, "accall")
  check = (os.path.exists(accall_exe) and os.access(accall_exe, os.X_OK) \
           and os.path.exists(os.path.join(naccess_root, "vdw.radii")) \
           and os.path.exists(os.path.join(naccess_root, "standard.data")))
  if not check:
    LogWarning("NACCESS: Could not find required files to launch accall " \
               "directly in %s." % naccess_root)
  return check

def _SetupFiles(entity, selection, scratch_dir, max_number_of_atoms):
  """
  Setup files for naccess calculation in temporary directory

  :param entity:              OST entity to calculate surface
  :param selection:           Calculate surface for subset of entity
  :param scratch_dir:         Scratch directory. A subfolder for temporary files
                              is created in there. If not specified, a default
                              directory is used (see :func:`tempfile.mkdtemp`).
  :param max_number_of_atoms: Max Number of atoms in the entity (i.e. is limited
                              in the default NACCESS version to 50 000)
  :returns: Tuple containing temporary directory, input filename for naccess and
            directory of the input file
  :exception: RuntimeError if selection is not valid or too many atoms
  """

  # create temporary directory
  tmp_dir_name = ""
  if scratch_dir != None:
    tmp_dir_name = tempfile.mkdtemp(dir=scratch_dir)
  else:
    tmp_dir_name = tempfile.mkdtemp()

  # select as specified by user
  if selection != "":
    entity_view = entity.Select(selection)
  else:
    entity_view = entity
  if len(entity_view.atoms) > max_number_of_atoms:
    raise RuntimeError("Too much atoms for NACCESS (> %s)" % max_number_of_atoms)
  if not entity_view.IsValid():
    raise RuntimeError("Could not create view for selection (%s)"%(selection))
  
  # write entity to tmp file
  tmp_file_name = "entity.pdb"
  tmp_file_base = os.path.join(tmp_dir_name, "entity")
  io.SavePDB(entity_view, os.path.join(tmp_dir_name, tmp_file_name))
  return (tmp_dir_name, tmp_file_name, tmp_file_base)


def _ParseAsaFile(entity, file, asa_atom):
  """
  Reads Area file (.asa) and attach asa per atom to an entitiy

  :param entity:   EntityHandle or EntityView for attaching sasa on atom level
  :param file:     Filename of area file
  :param asa_atom: Name of the float property for SASA
  """
  
  asa_fh = open(file)
  asa_lines = asa_fh.readlines()
  asa_fh.close()
  
  for l in  asa_lines:
    if l.startswith("ATOM"):
      # get res_number, chain_id and atom name
      atom_name = l[12:16]
      chain_id = l[21]
      res_number = l[22:27]
      asa = l[54:63]
      atom_name = atom_name.strip()
      chain_id = chain_id
      res_number = res_number.strip()
      asa = asa.strip()
      m = re.match(r'(?P<num>-?\d+)(?P<ins>\w)?', res_number)
      di = m.groupdict()
      
      if di["ins"] == None:
        resNum = mol.ResNum(int(di["num"]))
      else:
        resNum = mol.ResNum(int(di["num"]), di["ins"])

      a = entity.FindAtom(chain_id, resNum, atom_name)
      if(a.IsValid()):
        a.SetFloatProp(asa_atom, float(asa))
      else:
        LogWarning("NACCESS: invalid asa entry %s %s %s" \
                   % (chain_id, resNum, atom_name))
      
def _ParseRsaFile(entity, file, asa_abs, asa_rel):
  """
  Reads Area file (.rsa) and attach asa (absolute + relative) per residue to an entitiy

  :param entity:     EntityHandle or EntityView for attaching sasa on atom level
  :param file:       Filename of .rsa file
  :param asa_atom:   Name of the float property for absolute SASA
  :param asa_atom:   Name of the float property for relative SASA
  :exception:        RuntimeError if residue names are not the same
  """
  area_fh = open(file)
  area_lines = area_fh.readlines()
  area_fh.close()
  # shift first line
  area_lines = area_lines[4:]
  # parse lines
  for l in area_lines:
    if l.startswith("RES"):
      # extract data
      p = re.compile(r'\s+')
      res_name = l[3:8]
      res_name = res_name.strip()
      chain_id = l[8:9]
      res_number = l[9:14]
      res_number = res_number.strip()
      abs_all, rel_all = l[15:28].strip().split()
      m = re.match(r'(?P<num>-?\d+)(?P<ins>\w)?', res_number)
      di = m.groupdict()
      if di["ins"] == None:
        resNum = mol.ResNum(int(di["num"]))
      else:
        resNum = mol.ResNum(int(di["num"]), di["ins"])
      # set res. props
      res = entity.FindResidue(chain_id, resNum)
      if res_name == res.name:
        res.SetFloatProp(asa_rel, float(rel_all) )
        res.SetFloatProp(asa_abs, float(abs_all) )
      else:
        raise RuntimeError("Residue Names are not the same for ResNumb: %s (%s vs %s)" % (res_number, res.name, res_name))
      

def __CleanupFiles(dir_name):
  """
  Method which recursively deletes a directory

  :warning: This method removes also non-empty directories without asking, so
            be careful!
  """
  import shutil
  shutil.rmtree(dir_name)


def _RunACCALL(command, temp_dir, query):
  """
  Fast method to run the Naccess surface calculation.

  This method starts the accall binary directly and pipes in the input provided
  in *query*. This is faster than calling the "naccess" script since the script
  has a constant overhead of roughly 1.3s in each call.

  :param command:  Command to execute
  :param temp_dir: Command is executed with this working directory
  :param query:    User input to pipe into *command*
  :returns:        stdout of command
  :exception:      CalledProcessError for non-zero return value
  """

  proc = subprocess.Popen(command, cwd=temp_dir, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, stdin=subprocess.PIPE)
  stdout_value, stderr_value = proc.communicate(query.encode())

  # check for successful completion of naccess
  if proc.returncode != 0:
    LogWarning("WARNING: naccess error\n%s\n%s" % (stdout_value.decode(), 
                                                   stderr_value.decode()))
    raise subprocess.CalledProcessError(proc.returncode, command)

  return stdout_value.decode()


def _RunNACCESS(command, temp_dir):
  """
  Method to run the Naccess surface calculation.

  This method starts the external Naccess executable and returns the stdout.

  :param command:  Command to execute
  :param temp_dir: Command is executed with this working directory
  :returns:        stdout of command
  :exception:      CalledProcessError for non-zero return value
  """
  proc = subprocess.Popen(command, cwd=temp_dir, shell=True, 
                          stdout=subprocess.PIPE)
  stdout_value, stderr_value = proc.communicate()

  # check for successful completion of naccess
  if proc.returncode != 0:
    LogWarning("WARNING: naccess error\n%s" % stdout_value.decode())
    raise subprocess.CalledProcessError(proc.returncode, command)

  return stdout_value.decode()


def CalculateSurfaceArea(entity,  radius=1.4,  
                         include_hydrogens=False, include_hetatm = False, 
                         include_water=False, selection="", naccess_exe=None,
                         naccess_root=None, keep_files=False, asa_abs= "asaAbs",
                         asa_rel="asaRel", asa_atom="asaAtom", scratch_dir=None,
                         max_number_of_atoms=50000):
  """
  Calculates analytical the solvent accessible surface area by using the
  external naccess program

  This method calculates the molecular surface areas by invoking the external
  program naccess. First, it is checked if the naccess executable is present, then,
  the necessary files are prepared in a temporary directory and naccess is
  executed. The last step is to remove the temporary directory.


  :param entity:              OST entity to calculate surface
  :param radius:              Surface probe radius
  :param include_hydrogens:   Calculate surface including hydrogens
  :param include_hetatm:      Calculate surface including hetatms
  :param include_water:       Calculate surface including water
  :param selection:           Calculate surface for subset of entity
  :param naccess_exe:         naccess executable (full path to executable)
  :param naccess_root:        Path to folder containing "accall" binary and
                              files "vdw.radii" and "standard.data". This is the
                              fastest way to call naccess!
  :param keep_files:          If True, do not delete temporary files
  :param asa_abs:             Attaches per residue absolute SASA to specified
                              FloatProp on residue level
  :param asa_rel:             Attaches per residue relative SASA to specified
                              FloatProp on residue level
  :param asa_atom:            Attaches per atom SASA to specified FloatProp at
                              atom level
  :param scratch_dir:         Scratch directory. A subfolder for temporary files
                              is created in there. If not specified, a default
                              directory is used (see :func:`tempfile.mkdtemp`).
  :param max_number_of_atoms: Max Number of atoms in the entity (i.e. is limited
                              in the default NACCESS version to 50 000)

  :returns:                   absolute SASA calculated using asa_atom
  """
  
  # check if naccess executable is specified
  if naccess_root and _CheckNaccessRoot(naccess_root):
    # use faster, direct call to accall binary
    fast_mode = True
  else:
    # get naccess executable
    naccess_executable = _GetExecutable(naccess_exe)
    # see if we can extract naccess_root from there (fallback to old mode)
    naccess_root = os.path.dirname(naccess_executable)
    fast_mode = _CheckNaccessRoot(naccess_root)
  
  # setup files for naccess
  (naccess_data_dir, naccess_data_file, naccess_data_base) \
    = _SetupFiles(entity, selection, scratch_dir, max_number_of_atoms)

  try:
    # call naccess
    if fast_mode:
      # cook up stdin query (same logic as naccess script)
      query = "PDBFILE %s\n" \
              "VDWFILE %s\n" \
              "STDFILE %s\n" \
              "PROBE %f\n" \
              "ZSLICE 0.05\n" \
              % (naccess_data_file, os.path.join(naccess_root, "vdw.radii"),
                 os.path.join(naccess_root, "standard.data"), radius)
      if include_hydrogens:
        query += "HYDROGENS\n"
      if include_water:
        query += "WATERS\n"
      if include_hetatm:
        query += "HETATOMS\n"
      # call it
      command = os.path.join(naccess_root, "accall")
      _RunACCALL(command, naccess_data_dir, query)
    else:
      LogWarning("NACCESS: Falling back to slower call to %s." \
                 % naccess_executable)
      # set command line
      command = "%s %s -p %f " % \
                (naccess_executable, naccess_data_file, radius)
      if include_hydrogens:
        command = "%s -y" % command
      if include_water:
        command = "%s -w" % command
      if include_hetatm:
        command = "%s -h" % command
      # execute naccess
      _RunNACCESS(command, naccess_data_dir)
    
    # parse outout
    new_asa = os.path.join(naccess_data_dir, "%s.asa" % naccess_data_base) 
    _ParseAsaFile(entity, new_asa, asa_atom)
      
    new_rsa = os.path.join(naccess_data_dir, "%s.rsa" % naccess_data_base) 
    _ParseRsaFile(entity, new_rsa, asa_abs, asa_rel)

  finally:
    # clean up
    if not keep_files:
      __CleanupFiles(naccess_data_dir)
  
  # sum up Asa for all atoms
  sasa = 0.0
  for a in entity.atoms:
    sasa += a.GetFloatProp(asa_atom, 0.0)

  return sasa
