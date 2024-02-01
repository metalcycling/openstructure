'''
MSMS module

Author: Tobias Schmidt

This module is for calculating MSMS surfaces as well as surface areas
(SESA, SASA) from OpenStructure using the external program MSMS.

How To Use This Module:
 1. Import it (e.g. as "from ost.bindings import msms")
 2. Use it (e.g. as "surfaces_list = msms.CalculateSurface(entity)"
                    "(sesa,sasa) = msms.CalculateSurfaceArea(entity)")

Requirement:
 - MSMS installed
'''

import tempfile
import subprocess
import os

from ost import io
from ost import mol
from ost import settings
from ost import geom


class MsmsProcessError(Exception):
  """
  Python 2.4 and older do not include the CalledProcessError exception. This
  class substitutes it.
  """
  def __init__(self, returncode,command):
    self.returncode = returncode
    self.command = command
  def __str__(self):
    return repr(self.returncode)


def GetVersion(msms_exe=None, msms_env=None):
  """
  Get version of MSMS executable
  """
  msms_executable = _GetExecutable(msms_exe, msms_env)
  command = "%s" % (msms_executable)
  proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  stdout_value, stderr_value = proc.communicate()

  version = ""
  for l in stdout_value.decode().splitlines():
    if l[0:4]=='MSMS':
      version = l.split(' ')[1]
      return version
  if version=="":
    LogWarning('Could not parse MSMS version string')
    return


def _GetExecutable(msms_exe, msms_env):
  """
  Function to check if MSMS executable is present

  :param msms_exe: Explicit path to msms executable
  :param msms_env: Environment variable pointing to msms executable
  :returns: Path to the executable
  :raises:  :class:`~ost.FileNotFound` if executable is not found
  """
  return settings.Locate('msms', explicit_file_name=msms_exe,
                         env_name=msms_env)


def _SetupFiles(entity, selection):
  """
  Setup files for MSMS calculation in temporary directory

  :param entity: The entity for which the surface is to be calculated
  :type entity: :class:`~ost.mol.EntityHandle` or :class:`~ost.mol.EntityHandle`
  :param selection:  Calculate surface for subset of entity
  :type selection: :class:`str`
  :returns: tuple containing temporary directory and msms input file
  :raises:         :class:`RuntimeError` if selection is not valid
  """
  #   create temporary directory
  tmp_dir_name=tempfile.mkdtemp()

  # select only heavy atoms if no_hydrogens is true
  entity_view=entity.Select(selection)
  if not entity_view.IsValid():
    raise RuntimeError("Could not create view for selection (%s)"%(selection))

  # write entity to tmp file
  tmp_file_name=os.path.join(tmp_dir_name,"entity")
  tmp_file_handle=open(tmp_file_name, 'w')
  for a in entity_view.GetAtomList():
    position=a.GetPos()
    tmp_file_handle.write('%8.3f %8.3f %8.3f %4.2f\n' % (position[0],
                          position[1], position[2], a.radius))
  tmp_file_handle.close()
  
  return (tmp_dir_name, tmp_file_name)


def _ParseAreaFile(entity, selection, file, asa_prop, esa_prop):
  """
   Reads Area file (-af) and attach sasa and sesa per atom to an entitiy

  :param entity:   :class:`~ost.mol.EntityHandle` or :class:`~ost.mol.EntityView`
                   for attaching sasa and sesa on atom level
  :param file:     Filename of area file
  :param asa_prop: Name of the float property for SASA
  :param esa_prop: Name of the float property for SESA
  :raises: :class:`RuntimeError` if number of atoms in file != number of atoms in entity
  """
  view=entity.Select(selection)
  area_fh = open(file)
  area_lines = area_fh.readlines()
  area_fh.close()
  # shift first line
  area_lines = area_lines[1:]
  if view.GetAtomCount() != len(area_lines):
      raise RuntimeError("Atom count (%d) unequeal to number of atoms in area file (%d)" % (view.GetAtomCount(), len(area_lines)))
  for l in area_lines:
      atom_no, sesa, sasa = l.split()
      a = view.atoms[int(atom_no)]
      if asa_prop:
        a.SetFloatProp(asa_prop, float(sasa))
      if esa_prop:
        a.SetFloatProp(esa_prop, float(sesa))
    
    

def _CleanupFiles(dir_name):
  """
  Function which recursively deletes a directory and all the files contained
  in it. *Warning*: This method removes also non-empty directories without
  asking, so be careful!
  """
  import shutil
  shutil.rmtree(dir_name)

def _RunMSMS(command):
  """
  Run the MSMS surface calculation

  This functions starts the external MSMS executable and returns the stdout of
  MSMS.

  :param command:          Command to execute
  :returns:                 stdout of MSMS
  :raises:              :class:`CalledProcessError` for non-zero return value
  """
  proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  stdout_value, stderr_value = proc.communicate()

  #check for successful completion of msms
  if proc.returncode!=0:
    print("WARNING: msms error\n", stdout_value.decode())
    raise MsmsProcessError(proc.returncode, command)

  return stdout_value.decode()
  


def CalculateSurfaceArea(entity, density=1.0, radius=1.5,  all_surf=False,
                         no_hydrogens=False, no_hetatoms=False, no_waters=False,
                         selection='',
                         msms_exe=None, msms_env=None, keep_files=False, 
                         attach_asa=None, attach_esa=None):
  """
  Calculates analytical solvent excluded and solvent accessible surface
  area by using the external MSMS program.

  This method calculates the molecular surface areas by invoking the external
  program MSMS. First, it is checked if the MSMS executable is present, then,
  the necessary files are prepared in a temporary directory and MSMS is
  executed. The last step is to remove the temporary directory.


  :param entity:        OST entity to calculate surface
  :param density:       Surface point density
  :param radius:       Surface probe radius
  :param all_surf:      Calculate surface area for all cavities (returns multiple
      surfaces areas as a list)
  :param no_hydrogens:  Calculate surface only for hevy atoms
  :param selection:     Calculate surface for subset of entity
  :param msms_exe:      msms executable (full path to executable)
  :param msms_env:      msms environment variable
  :param keep_files:    Do not delete temporary files
  :param attach_asa:    Attaches per atom SASA to specified FloatProp at atom level
  :param attach_esa:    Attaches per atom SESA to specified FloatProp at atom level
  :returns:             Tuple of lists for (SES, SAS)
  """
  import re 

  # check if msms executable is specified
  msms_executable=_GetExecutable(msms_exe, msms_env)

  # parse selection
  if no_hydrogens:
    if selection!='':
      selection+=" and "
    selection+="ele!=H"
  
  if no_hetatoms:
    if selection!='':
      selection+=" and "
    selection+="ishetatm=False"
  
  if no_waters:
    if selection!='':
      selection+=" and "
    selection+="rname!=HOH"

  # setup files for msms
  (msms_data_dir, msms_data_file)=_SetupFiles(entity, selection)

  # set command line
  command="%s -if %s -of %s -density %s -probe_radius %s -surface ases" % \
          (msms_executable, msms_data_file, msms_data_file, density, radius)
  if all_surf:
    command+=" -all"
  if attach_asa != None or attach_esa != None:
    command+=" -af %s" % os.path.join(msms_data_dir, "asa_atom")
  # run msms
  stdout_value=_RunMSMS(command)
  
  # add sesa and asa to entity if attach_asa is specified
  if attach_asa != None or attach_esa != None:
      _ParseAreaFile(entity, selection, os.path.join(msms_data_dir, "asa_atom.area"),
                     attach_asa, attach_esa)
  
  # parse MSMS output
  msms_ases=[]
  msms_asas=[]
  data_paragraph=False
  for line in stdout_value.splitlines():
    if re.match('MSMS terminated normally', line):
      data_paragraph=False
    if data_paragraph:
      (ses_,sas_)=line.split()[5:7]
      msms_ases.append(float(ses_))
      msms_asas.append(float(sas_))
    if re.match('    Comp. probe_radius,   reent,    toric,   contact    SES       SAS', line):
      data_paragraph=True

  # clean up
  if not keep_files:
    _CleanupFiles(msms_data_dir)

  return (msms_ases, msms_asas)
  
def CalculateSurfaceVolume(entity, density=1.0, radius=1.5,  all_surf=False,
                         no_hydrogens=False, no_hetatoms=False, no_waters=False,
                         selection='',
                         msms_exe=None, msms_env=None, keep_files=False, 
                         attach_asa=None, attach_esa=None):
  """
  Calculates the volume of the solvent excluded surface by using the external MSMS program.

  This method calculates the volume of the molecular surface by invoking the external
  program MSMS. First, it is checked if the MSMS executable is present, then,
  the necessary files are prepared in a temporary directory and MSMS is
  executed. The last step is to remove the temporary directory.


  :param entity:        OST entity to calculate surface
  :param density:       Surface point density
  :param radius:       Surface probe radius
  :param all_surf:      Calculate surface area for all cavities (returns multiple
      surfaces areas as a list)
  :param no_hydrogens:  Calculate surface only for hevy atoms
  :param selection:     Calculate surface for subset of entity
  :param msms_exe:      msms executable (full path to executable)
  :param msms_env:      msms environment variable
  :param keep_files:    Do not delete temporary files
  :param attach_asa:    Attaches per atom SASA to specified FloatProp at atom level
  :param attach_esa:    Attaches per atom SESA to specified FloatProp at atom level
  :returns:             Tuple of lists for (SES, SAS)
  """
  import re 

  # check if msms executable is specified
  msms_executable=_GetExecutable(msms_exe, msms_env)

  # parse selection
  if no_hydrogens:
    if selection!='':
      selection+=" and "
    selection+="ele!=H"
  
  if no_hetatoms:
    if selection!='':
      selection+=" and "
    selection+="ishetatm=False"
  
  if no_waters:
    if selection!='':
      selection+=" and "
    selection+="rname!=HOH"

  # setup files for msms
  (msms_data_dir, msms_data_file)=_SetupFiles(entity, selection)

  # set command line
  command="%s -if %s -of %s -density %s -probe_radius %s " % \
          (msms_executable, msms_data_file, msms_data_file, density, radius)
  if all_surf:
    command+=" -all"
  if attach_asa != None or attach_esa != None:
    command+=" -af %s" % os.path.join(msms_data_dir, "asa_atom")
  # run msms
  stdout_value=_RunMSMS(command)
  
  # add sesa and asa to entity if attach_asa is specified
  if attach_asa != None or attach_esa != None:
      _ParseAreaFile(entity, selection, os.path.join(msms_data_dir, "asa_atom.area"),
                     attach_asa, attach_esa)
  
  # parse MSMS output
  ses_volume=0
  for line in stdout_value.splitlines():
    if re.match('    Total ses_volume:', line):
      ses_volume=float(line.split(':')[1])

  # clean up
  if not keep_files:
    _CleanupFiles(msms_data_dir)

  return ses_volume


def CalculateSurface(entity, density=1.0, radius=1.5, all_surf=False,
                     no_hydrogens=False, no_hetatoms=False, no_waters=False,
                     selection='',
                     msms_exe=None, msms_env=None, keep_files=False,
                     attach_asa=None, attach_esa=None):
  
  """
  Calculates molecular surface by using the external MSMS program

  This method calculates a molecular surface by invoking the external program
  MSMS. First, it is checked if the MSMS executable is present, then, the
  necessary files are prepared in a temporary directory and MSMS is executed.
  The last step is to remove the temporary directory.


  :param entity:        Entity for which the surface is to be calculated
  :param density:       Surface point density
  :param radius:        Surface probe radius
  :param all_surf:      Calculate surface for all cavities (returns multiple
                        surfaces as a list)
  :param no_hydrogens:  Calculate surface only for heavy atoms
  :param selection:     Calculate surface for subset of entity
  :param msms_exe:      msms executable (full path to executable)
  :param msms_env:      msms environment variable
  :param keep_files:    Do not delete temporary files
  :param attach_asa:    Attaches per atom SASA to specified FloatProp at atom level
  :param attach_esa:    Attaches per atom SESA to specified FloatProp at atom level
  :returns:             list of :class:`~ost.mol.SurfaceHandle` objects
  """
  import os
  import re

  # check if msms executable is specified
  msms_executable=_GetExecutable(msms_exe, msms_env)

  # parse selection
  if no_hydrogens:
    if selection!='':
      selection+=" and "
    selection+="ele!=H"
      
  if no_hetatoms:
    if selection!='':
      selection+=" and "
    selection+="ishetatm=False"
      
  if no_waters:
    if selection!='':
      selection+=" and "
    selection+="rname!=HOH"

  # setup files for msms
  (msms_data_dir, msms_data_file)=_SetupFiles(entity, selection)

  # set command line
  command="%s -if %s -of %s -density %s -probe_radius %s" % (msms_executable,
          msms_data_file, msms_data_file, density, radius)
  if all_surf:
    command+=" -all"
  if attach_asa != None or attach_esa != None:
    command+=" -af %s" % os.path.join(msms_data_dir, "asa_atom")

  # run msms
  stdout_value=_RunMSMS(command)

  # add sesa and asa to entity if attach_asa is specified
  if attach_asa != None or attach_esa != None:
      _ParseAreaFile(entity, selection, os.path.join(msms_data_dir, "asa_atom.area"),
                     attach_asa, attach_esa)

  # parse msms output
  num_surf=0
  for line in stdout_value.splitlines():
    if re.search('RS component [0-9]+ identified',line):
      num_surf=int(line.split()[2])

  # get surfaces
  entity_sel = entity.Select(selection)
  msms_surfaces=[]
  s = io.LoadSurface(msms_data_file, "msms")
  s.Attach(entity_sel, 3+radius)
  msms_surfaces.append(s)
  for n in range(1,num_surf+1):
    filename=msms_data_file+'_'+str(n)
    s = io.LoadSurface(filename, "msms")
    s.Attach(entity_sel, 3+radius)
    msms_surfaces.append(s)

  # clean up
  if not keep_files:
    _CleanupFiles(msms_data_dir)

  return msms_surfaces

