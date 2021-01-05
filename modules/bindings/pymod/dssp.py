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
Utility functions to load secondary structure information from DSSP files
and assign them to entities.

Authors: Pascal Benkert, Marco Biasini
"""

import os
import tempfile,subprocess

from ost import io,mol
from ost import settings


def _SkipHeader(stream):
    line=stream.readline()
    while line:
      if line.strip().find('#')==0:
        return True
      line=stream.readline()
    return False

def _Cleanup(temp_dssp_path, pdb_path):
  if os.path.exists(temp_dssp_path):
    os.remove(temp_dssp_path)
  if os.path.exists(pdb_path):
    os.remove(pdb_path)

def _CalcRelativeSA(residue_type, absolute_sa):
  solvent_max_list=[118,317,238,243,183,262,286,154,258,228,
                    243,278,260,271,204,234,206,300,303,216] #TODO: source?
  residue_indices = "ARNDCQEGHILKMFPSTWYV"
  # cysteine bridges are marked with lower-case letters by DSSP. We don't 
  # really care which cysteines are forming covalent bonds, so let's set the 
  # one-letter-code to "C".
  if residue_type.islower():
    residue_type='C'
  if residue_indices.find(residue_type)==-1:
    raise RuntimeError('residue %s is a non-standard residue' %(residue_type))	
  else:
    rel=float(absolute_sa)/(solvent_max_list[residue_indices.find(residue_type)])
  return rel
  

def AssignDSSP(ent, pdb_path="", extract_burial_status=False, tmp_dir=None, 
               dssp_bin=None):
  """
  Assign secondary structure states to peptide residues in the structure. This
  function uses the DSSP command line program.

  If you already have a DSSP output file and would like to assign the secondary 
  structure states to an entity, use :func:`LoadDSSP`.
  
  :param ent: The entity for which the secondary structure should be calculated
  :type ent: :class:`~ost.mol.EntityHandle` or :class:`~ost.mol.EntityView`
  :param extract_burial_status: If true, also extract burial status and store
                                as float-property
                                ``relative_solvent_accessibility`` at residue
                                level
  :param tmp_dir: If set, overrides the default tmp directory of the
                  operating system
  :param dssp_bin: The path to the DSSP executable
  :raises: :class:`~ost.settings.FileNotFound` if the dssp executable is not 
      in the path.
  :raises: :class:`RuntimeError` when dssp is executed with errors
  """

  if not ent.IsValid():
    raise ValueError('model entity is not valid')
  if ent.atom_count==0:
    raise ValueError('model entity does not contain any atoms')

  dssp_abs_path=settings.Locate(['dsspcmbi','dssp','mkdssp'], env_name='DSSP_EXECUTABLE', 
                                explicit_file_name=dssp_bin)
  if os.path.isdir(dssp_abs_path):
    raise RuntimeError('"%s" is a directory. Specify path to DSSP binary' % dssp_abs_path)
  if not os.access(dssp_abs_path, os.X_OK):
    raise RuntimeError('"%s" is not executable' % dssp_abs_path)
    
  pdb_path=tempfile.mktemp(suffix=".pdb",prefix="temp_entity", dir=tmp_dir)
  io.SavePDB(ent, pdb_path)
  temp_dssp_path=tempfile.mktemp(suffix=".out",prefix="dssp", dir=tmp_dir)  

  cmd = [dssp_abs_path, pdb_path, temp_dssp_path]
  s = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
  if s.returncode == 0:
    try:
      LoadDSSP(temp_dssp_path, ent, extract_burial_status)
      _Cleanup(temp_dssp_path, pdb_path)
      return ent
    except:
      pass # continue with dssp 4 fallback
  
  # either dssp execution failed or output file cannot be loaded
  # both can be the result of using dssp 4 (https://github.com/PDB-REDO/dssp) 
  # => try fallback

  # dssp 4 desperately wants a CRYST1 record (REMOVE IF NOT NEEDED ANYMORE!)
  with open(pdb_path, 'r') as fh:
    file_content = fh.read()
  with open(pdb_path, 'w') as fh:
    fh.write('CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          ' +
             os.linesep + file_content)

  # explicitely request dssp output format
  cmd = [dssp_abs_path, '--output-format', 'dssp', pdb_path, temp_dssp_path]
  s_fallback = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
  
  if s_fallback.returncode == 0:
    try:
      LoadDSSP(temp_dssp_path, ent, extract_burial_status)
      _Cleanup(temp_dssp_path, pdb_path)
      return ent
    except:
      pass # continue with cleanup and raise error
      
  _Cleanup(temp_dssp_path, pdb_path)
  raise RuntimeError('Failed to assign DSSP')


def LoadDSSP(file_name, model, extract_burial_status=False,
             entity_saved=False, calculate_relative_sa=True):
    """
    Loads DSSP output and assigns secondary structure states to the peptidic
    residues.

    If you would like to run dssp *and* assign the secondary structure,
    use :func:`AssignDSSP` instead.

    :param file_name: The filename of the DSSP output file
    :param model: The entity to which the secondary structure states should be
                  assigned
    :param extract_burial_status: If true also calculates burial status of
        residues and assigns it to the burial_status string property.
    :param calculate_relative_sa: If true also relative solvent accessibility and
        and assigns it to the relative_solvent_accessibility float property of
        the residue.
    :param entity_save: Whether the entity was saved.
    """
    if not model.IsValid():
      raise ValueError('model entity is not valid')
    if model.atom_count==0:
      raise ValueError('model entity does not contain any atoms')
    stream=open(file_name)
    if not _SkipHeader(stream):
      stream.close()
      raise RuntimeError('Ill-formatted DSSP file')

    for line in stream:
      num=line[6:10].strip()
      ins_code=line[10].strip()
      chain_name=line[11]
      solvent_accessibility=float(line[34:39].strip())
      #solvent_accessibility=line[34:39].strip()
      amino_acid=line[13]
      #print line


      if isinstance(model,mol.ChainView):
        chain=model
      else:
        chain=model.FindChain(chain_name)
        
      if not chain.IsValid():
        continue
      if num=='':
        continue
      residue=None
      try:
        if ins_code == "":
          residue=chain.FindResidue(mol.ResNum(int(num)))
        else:
          residue=chain.FindResidue(mol.ResNum(int(num),ins_code))

        # set property "burial status:
        if extract_burial_status:
          #set default (dummy) burial status for incomplete residues:
          residue.SetStringProp("burial_status", 'X')

          #handle seleno-methionine appearing as amino acid 'X' in DSSP:
          if residue.name=="MSE" and amino_acid=='X':
            amino_acid='M'

          residue.SetFloatProp("solvent_accessibility",
                                          solvent_accessibility)
          if calculate_relative_sa:
            relative_sa=_CalcRelativeSA(amino_acid,solvent_accessibility)
            residue.SetFloatProp("relative_solvent_accessibility",
                                            relative_sa)
            if relative_sa < 0.25:
              residue.SetStringProp("burial_status", 'b')
            else:
              residue.SetStringProp("burial_status", 'e')
      except Exception as e:
        print("ERROR:",e)
        continue

      rtype=line[16:17]
      rt=mol.SecStructure.COIL
      if rtype=='H':
        rt=mol.SecStructure.ALPHA_HELIX
      elif rtype=='E':
        rt=mol.SecStructure.EXTENDED
      elif rtype=='B':
        rt=mol.SecStructure.BETA_BRIDGE
      elif rtype=='S':
        rt=mol.SecStructure.BEND
      elif rtype=='T':
        rt=mol.SecStructure.TURN
      elif rtype=='I':
        rt=mol.SecStructure.PI_HELIX
      elif rtype=='G':
        rt=mol.SecStructure.THREE_TEN_HELIX
      # for corrupted DSSP files. Catch in calling routine:
      if not residue.IsValid():
        #Todo: if residues with occupancy 0 have been removed before
        #using a selection statement, they are missed here
        #IMPORTANT: asign DSSP before any selections
        stream.close()
        raise RuntimeError('Ill-formatted DSSP file: invalid residue')	
      residue.SetSecStructure(mol.SecStructure(rt))
    stream.close()
