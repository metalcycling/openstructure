#------------------------------------------------------------------------------
# This file is part of the OpenStructure project <www.openstructure.org>
#
# Copyright (C) 2008-2020 by the OpenStructure authors
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------
"""
Wrapper for the CAD score.

References:

Olechnovic K, Kulberkyte E, Venclovas C., CAD-score: A new contact area 
difference-based function for evaluation of protein structural models
Proteins. 2012 Aug 30. [Epub ahead of print]

Authors: Valerio Mariani, Alessandro Barbato
"""

import subprocess, os, tempfile, platform, re
from ost import settings, io, mol

def _SetupFiles(model,reference):
  # create temporary directory
  tmp_dir_name=tempfile.mkdtemp()
  dia = 'PDB'
  for chain in model.chains:
    if chain.name==" ":
      raise RuntimeError("One of the chains in the model has no name. Cannot "
                         "calculate CAD score")
    if len(chain.name) > 1:
      dia = 'CHARMM'
      break;
    for res in chain.residues:
      if len(res.name) > 3:
        dia = 'CHARMM'
        break;
  io.SavePDB(model, os.path.join(tmp_dir_name, 'model.pdb'), dialect=dia)
  dia = 'PDB'
  for chain in reference.chains:
    if chain.name==" ":
      raise RuntimeError("One of the chains in the reference has no name. "
                         "Cannot calculate CAD score")
    if len(chain.name) > 1:
      dia = 'CHARMM'
      break;
    for res in chain.residues:
      if len(res.name) > 3:
        dia = 'CHARMM'
        break;
  io.SavePDB(reference, os.path.join(tmp_dir_name, 'reference.pdb'),dialect=dia)
  return tmp_dir_name

def _CleanupFiles(dir_name):
  import shutil
  shutil.rmtree(dir_name)

class CADResult:
  """
  Holds the result of running CAD
  
  .. attribute:: globalAA
    
    The global CAD's atom-atom (AA) score
  
  .. attribute:: localAA
  
    Dictionary containing local CAD's atom-atom (AA) scores. 
    
    :type: dictionary (key: tuple(chain, resnum) (e.g.: 
           ("A", ost.mol.ResNum(24)), value: CAD local AA score 
           (see CAD Documentation online)
  """
  def __init__(self, globalAA, localAA):    
    self.globalAA=globalAA
    self.localAA=localAA    

def _ParseCADGlobal(lines):
  header = lines[0].split()
  aa_idx = header.index("AA")
  aa_score=float(lines[1].split()[aa_idx])
  return aa_score

def _ParseCADLocal(lines):
  local_scores_idx = None
  for line_idx in range(len(lines)):
    if "local_scores" in lines[line_idx]:
      local_scores_idx = line_idx
      break
  if local_scores_idx == None:
    raise RuntimeError("Failed to parse local cadscores")
  local_aa_dict={}
  for line_idx in range(local_scores_idx+2, len(lines)):
    items=lines[line_idx].split()
    local_aa = float(items[2])
    if local_aa < 0.0:
      continue # invalid CAD score
    key = (items[0], mol.ResNum(int(items[1])))
    local_aa_dict[key] = local_aa
    
  return local_aa_dict

def _ParseVoronotaGlobal(lines):
  return float(lines[0].split()[4])

def _ParseVoronotaLocal(lines):
  local_aa_dict={}
  chain_name_regex = r'c\<\D+\>'
  resnum_regex = r'r\<\d+\>'
  insertion_code_regex = r'i\<\D\>'
  for line in lines:
    local_aa = float(line.split()[-1])
    if local_aa < 0.0:
      continue # invalid CAD score
    chain_data = re.findall(chain_name_regex, line)
    resnum_data = re.findall(resnum_regex, line)
    insertion_code_data = re.findall(insertion_code_regex, line)
    resnum = None
    if len(insertion_code_data) == 0:
      resnum = mol.ResNum(int(resnum_data[0][1:].strip('><')))
    else:
      resnum = mol.ResNum(int(resnum_data[0][1:].strip('><')), 
                          insertion_code_data[0][1:].strip('><'))
    key = (chain_data[0][1:].strip('><'), resnum)
    local_aa_dict[key] = local_aa
  return local_aa_dict

def _RunCAD(tmp_dir, mode, cad_bin_path, old_regime):

  model_filename=os.path.join(tmp_dir, 'model.pdb')
  reference_filename=os.path.join(tmp_dir, 'reference.pdb')
  globalAA = None
  localAA = None

  if platform.system() == "Windows":
    raise RuntimeError('CAD score not available on Windows')

  if mode == "classic":
    cad_calc_path = None
    cad_read_g_path = None
    cad_read_l_path = None
    if cad_bin_path:
      cad_calc_path = settings.Locate('CADscore_calc.bash', 
                                      search_paths=[cad_bin_path],
                                      search_system_paths=False)  
      cad_read_g_path = settings.Locate('CADscore_read_global_scores.bash', 
                                        search_paths=[cad_bin_path],
                                        search_system_paths=False)  
      cad_read_l_path=settings.Locate('CADscore_read_local_scores.bash', 
                                      search_paths=[cad_bin_path],
                                      search_system_paths=False)
      # also try to locate the actual executable that is called from the 
      # bash scripts 
      executable_path = settings.Locate('voroprot2', 
                                        search_paths=[cad_bin_path],
                                        search_system_paths=False) 
    else:
      cad_calc_path = settings.Locate('CADscore_calc.bash')  
      cad_read_g_path = settings.Locate('CADscore_read_global_scores.bash')  
      cad_read_l_path = settings.Locate('CADscore_read_local_scores.bash')
      # also try to locate the actual executable that is called from the 
      # bash scripts
      executable_path = settings.Locate('voroprot2')  
    command1="\"%s\" -m \"%s\" -t \"%s\" -D \"%s\"" %(cad_calc_path, 
                                                      model_filename, 
                                                      reference_filename, 
                                                      os.path.join(tmp_dir,
                                                                   "cadtemp"))
    command2="\"%s\" -D \"%s\"" %(cad_read_g_path, os.path.join(tmp_dir,
                                                                "cadtemp"))
    command3="\"%s\" -m \"%s\" -t \"%s\" -D \"%s\" -c AA" %(cad_read_l_path, 
                                                            model_filename, 
                                                            reference_filename,
                                                            os.path.join(tmp_dir,
                                                            "cadtemp"))

    ps1=subprocess.Popen(command1, shell=True, stdout=subprocess.PIPE)
    ps1.communicate()
    ps2=subprocess.Popen(command2, shell=True, stdout=subprocess.PIPE)
    stdout,_ = ps2.communicate()
    lines=stdout.decode().splitlines()
    try:
      globalAA=_ParseCADGlobal(lines)
    except:
      raise RuntimeError("CAD calculation failed")
    ps3=subprocess.Popen(command3, shell=True, stdout=subprocess.PIPE)
    stdout,_ = ps3.communicate()
    lines=stdout.decode().splitlines()
    try:
      localAA=_ParseCADLocal(lines)
    except:
      raise RuntimeError("CAD calculation failed")

  elif mode == "voronota":
    local_score_filename = os.path.join(tmp_dir, "local_scores.txt")
    voronota_cadscore_path = None
    if cad_bin_path:
      voronota_cadscore_path = settings.Locate("voronota-cadscore", 
                                               search_paths=[cad_bin_path],
                                               search_system_paths=False)
      # also try to locate the actual executable that is called from the 
      # bash script
      executable_path = settings.Locate("voronota", 
                                        search_paths=[cad_bin_path],
                                        search_system_paths=False)      
    else:
      voronota_cadscore_path = settings.Locate("voronota-cadscore")
      # also try to locate the actual executable that is called from the 
      # bash script
      executable_path = settings.Locate("voronota")  
    cmd = [voronota_cadscore_path, '-m', model_filename, '-t', 
           reference_filename, '--contacts-query-by-code AA', 
           '--output-residue-scores', local_score_filename]
    if old_regime:
      cmd.append("--old-regime")
    cmd = ' '.join(cmd)
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    stdout, _ = ps.communicate()
    try:
      globalAA = _ParseVoronotaGlobal(stdout.decode().splitlines())
    except:
      raise RuntimeError("CAD calculation failed")
    try:
      with open(local_score_filename) as f: 
        localAA = _ParseVoronotaLocal(f.readlines())
    except:
      raise RuntimeError("CAD calculation failed")

  else:
    raise RuntimeError("Invalid CAD mode! Allowed are: "
                       "[\"classic\", \"voronota\"]")


  return CADResult(globalAA,localAA)

def _HasInsertionCodes(model, reference):
  for r in model.residues:
    if r.GetNumber().GetInsCode() != "\0":
      return True
  for r in reference.residues:
    if r.GetNumber().GetInsCode() != "\0":
      return True
  return False

def _MapLabels(model, cad_results, label):
  for k,v in cad_results.localAA.items():
    r = model.FindResidue(k[0], k[1])
    if not r.IsValid():
      raise RuntimeError("Failed to map cadscore on residues: " +
                         "CAD score estimated for residue in chain \"" +
                         k[0] + "\" with ResNum " + str(k[1]) + ". Residue " +
                         "could not be found in model.")
    r.SetFloatProp(label, v)

def CADScore(model, reference, mode = "voronota", label = "localcad",
             old_regime = False, cad_bin_path = None):
  """
  Calculates global and local atom-atom (AA) CAD Scores. 

  You can either access the original implementation available from
  https://bitbucket.org/kliment/cadscore/downloads/
  or the new implementation which is part of the Voronota package 
  available from https://bitbucket.org/kliment/voronota/downloads/.

  The scores of the two implementations differ but strongly correlate
  as the contacts between atoms are estimated differently. When using
  the "voronota" *mode* you can minimize those discrepancies by
  setting the *old_regime* flag to True.

  Furthermore, the "voronota" *mode* generates per-residue scores that 
  are inverted when compared to the classical implementation 
  (0.0: bad, 1.0 good). 

  :param model: The model structure. 
  :type model: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param reference: The reference structure
  :type reference: :class:`~ost.mol.EntityView` or 
                   :class:`~ost.mol.EntityHandle`
  :param mode:  What CAD score implementation to use, must be one in 
                ["classic", "voronota"]
  :param label: Local CAD scores will be mapped on residues of model as 
                float property with this name
  :type label: :class:`str`
  :param old_regime: Only has an effect if *mode* is "voronota". If set to true,
                     the discrepancies between the two modes is minimized but
                     the behaviour of inverted scores persists.
  :type old_regime: :class:`bool`
  :param cad_bin_path: Path to search for the required executables 
                       (["CADscore_calc.bash", 
                       "CADscore_read_global_scores.bash",
                       "CADscore_read_local_scores.bash"] for "classic" *mode* 
                       or ["voronota-cadscore"] for "voronota" *mode*). If not
                       set, the env path is searched.
  :type cad_bin_path: :class:`str`
  :returns: The result of the CAD score calculation
  :rtype: :class:`CADResult`
  
  :raises: :class:`~ost.settings.FileNotFound` if any of the CAD score 
           executables could not be located.
  :raises: :class:`RuntimeError` if the calculation failed
  """
  if mode == "classic" and _HasInsertionCodes(model, reference):
    raise RuntimeError("The classic CAD score implementation does not support "
                       "insertion codes in residues")
  tmp_dir_name=_SetupFiles(model, reference)
  result=_RunCAD(tmp_dir_name, mode, cad_bin_path, old_regime)
  _CleanupFiles(tmp_dir_name)
  _MapLabels(model, result, label)
  return result
