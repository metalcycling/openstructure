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
Wrappers for the tmalign and tmscore utilities.

References:

tmscore: Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710 
tmalign: Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9 


Authors: Pascal Benkert, Marco Biasini
"""

import subprocess, os, tempfile, platform
import ost
from ost import settings, io, geom, seq

def _SetupFiles(models):
  # create temporary directory
  tmp_dir_name=tempfile.mkdtemp()
  dia = 'PDB'
  for index, model in enumerate(models):
    for chain in model.chains:
      if len(chain.name) > 1:
        dia = 'CHARMM'
        break;
      for res in chain.residues:
        if len(res.name) > 3:
          dia = 'CHARMM'
          break;
    io.SavePDB(model, os.path.join(tmp_dir_name, 'model%02d.pdb' % (index+1)), dialect=dia)
  return tmp_dir_name

def _CleanupFiles(dir_name):
  import shutil
  shutil.rmtree(dir_name)

def _ParseTmAlign(lines,lines_matrix):
  info_line=lines[12].split(',')
  aln_length=int(info_line[0].split('=')[1].strip())
  rmsd=float(info_line[1].split('=')[1].strip())
  tm_score_swapped=float(lines[13].split('=')[1].split('(')[0].strip())
  tm_score=float(lines[14].split('=')[1].split('(')[0].strip())
  tf1=[float(i.strip()) for i in lines_matrix[2].split()]
  tf2=[float(i.strip()) for i in lines_matrix[3].split()]
  tf3=[float(i.strip()) for i in lines_matrix[4].split()]
  rot=geom.Mat3(tf1[2], tf1[3], tf1[4], tf2[2], tf2[3],
                tf2[4], tf3[2], tf3[3], tf3[4])
  tf=geom.Mat4(rot)
  tf.PasteTranslation(geom.Vec3(tf1[1], tf2[1], tf3[1]))
  seq1 = seq.CreateSequence("1",lines[18].strip())
  seq2 = seq.CreateSequence("2",lines[20].strip())
  alignment = seq.CreateAlignment()
  alignment.AddSequence(seq2)
  alignment.AddSequence(seq1)
  return ost.bindings.TMAlignResult(rmsd, tm_score, tm_score_swapped,
                                    aln_length, tf, alignment)

def _ParseUSAlign(lines,lines_matrix):

  # stuff that is immediately parsed
  rmsd = None
  tm_score = None
  tm_score_swapped = None
  aligned_length = None

  # first goes into intermediate data structures
  aln_data = list()
  mapping_data1 = list()
  mapping_data2 = list()
  in_aln = False

  for line in lines:
    if in_aln:
      if len(line.strip()) == 0:
        in_aln = False
      else:
        aln_data.append(line.strip('*'))
    elif line.startswith("Name of Structure_1:"):
      tmp = [item.strip() for item in line.split()[3].split(':')[1:]]
      for item in tmp:
        if len(item) > 0:
          mapping_data1.append(item.split(',')[1])
        else:
          mapping_data1.append("")
    elif line.startswith("Name of Structure_2:"):
      tmp = [item.strip() for item in line.split()[3].split(':')[1:]]
      for item in tmp:
        if len(item) > 0:
          mapping_data2.append(item.split(',')[1])
        else:
          mapping_data2.append("")
    elif line.startswith("Aligned length="):
      data = [item.strip() for item in line.split(',')]
      for item in data:
        if item.startswith("Aligned length="):
          aligned_length = int(item.split("=")[1])
        elif item.startswith("RMSD="):
          rmsd = float(item.split("=")[1])
    elif line.startswith("TM-score="):
      if "(normalized by length of Structure_1" in line:
        tm_score_swapped = float(line.split('(')[0].split('=')[1].strip())
      elif "(normalized by length of Structure_2" in line:
        tm_score = float(line.split('(')[0].split('=')[1].strip())
    elif line.startswith("(\":\" denotes residue pairs of"):
      in_aln = True

  assert(len(aln_data)==3)
  aln_sequences1 = aln_data[0].split('*')
  aln_sequences2 = aln_data[2].split('*')

  # do mapping/aln data
  alns = ost.seq.AlignmentList()
  ent1_mapped_chains = ost.StringList()
  ent2_mapped_chains = ost.StringList()
  assert(len(mapping_data1) == len(mapping_data2))
  assert(len(aln_sequences1) == len(aln_sequences2))
  assert(len(mapping_data1) == len(aln_sequences1))
  for a, b, c, d in zip(mapping_data1, mapping_data2,
                        aln_sequences1, aln_sequences2):
    if len(a) > 0 and len(b) > 0:
      ent1_mapped_chains.append(a)
      ent2_mapped_chains.append(b)
      assert(len(c) == len(d))
      aln = seq.CreateAlignment()
      aln.AddSequence(seq.CreateSequence(a, c))
      aln.AddSequence(seq.CreateSequence(b, d))
      alns.append(aln)

  # parse transformation matrix
  tf1=[float(i.strip()) for i in lines_matrix[2].split()[1:]]
  tf2=[float(i.strip()) for i in lines_matrix[3].split()[1:]]
  tf3=[float(i.strip()) for i in lines_matrix[4].split()[1:]]
  mat = geom.Mat4(tf1[0], tf1[1], tf1[2], tf1[3],
                  tf2[0], tf2[1], tf2[2], tf2[3],
                  tf3[0], tf3[1], tf3[2], tf3[3],
                  0.0,    0.0,    0.0,    1.0)

  return ost.bindings.MMAlignResult(rmsd, tm_score, tm_score_swapped,
                                    aligned_length, mat, alns,
                                    ent1_mapped_chains, ent2_mapped_chains)

def _RunTmAlign(tmalign, tmp_dir):
  model1_filename=os.path.join(tmp_dir, 'model01.pdb')
  model2_filename=os.path.join(tmp_dir, 'model02.pdb')
  if platform.system() == "Windows":
    tmalign_path=settings.Locate('tmalign.exe', explicit_file_name=tmalign)
    command="\"%s\" %s %s -m %s" %(os.path.normpath(tmalign_path), model1_filename, model2_filename, os.path.join(tmp_dir,'matrix.txt'))
  else:
    tmalign_path=settings.Locate('tmalign', explicit_file_name=tmalign)  
    command="\"%s\" \"%s\" \"%s\" -m \"%s\"" %(tmalign_path, model1_filename, model2_filename, os.path.join(tmp_dir,'matrix.txt'))
  ps=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  stdout,_=ps.communicate()
  lines=stdout.decode().splitlines()
  if (len(lines))<22:
    _CleanupFiles(tmp_dir)
    raise RuntimeError("tmalign superposition failed")
  matrix_file=open(os.path.join(tmp_dir,'matrix.txt'))
  lines_matrix=matrix_file.readlines()
  matrix_file.close() 
  return _ParseTmAlign(lines,lines_matrix)

def _RunUSAlign(usalign, tmp_dir):
  model1_filename=os.path.join(tmp_dir, 'model01.pdb')
  model2_filename=os.path.join(tmp_dir, 'model02.pdb')
  mat_filename = os.path.join(tmp_dir, "mat.txt")
  usalign_path=settings.Locate('USalign', explicit_file_name=usalign)  
  command = f"{usalign_path} {model1_filename} {model2_filename}  -mm 1 -ter 0 -m {mat_filename}"
  ps=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  stdout,_=ps.communicate()
  lines=stdout.decode().splitlines()
  if (len(lines))<22:
    _CleanupFiles(tmp_dir)
    raise RuntimeError("USalign superposition failed")
  with open(mat_filename) as fh:
    lines_matrix = fh.readlines()
  return _ParseUSAlign(lines,lines_matrix)

class TMScoreResult:
  """
  Holds the result of running TMscore
  
  .. attribute:: rmsd_common
    
    The RMSD of the common Calpha atoms of both structures

  .. attribute:: rmsd_below_five

    The RMSD of all Calpha atoms that can be superposed below five Angstroem
    
  .. attribute:: tm_score
  
    The TM-score of the structural superposition
  
  .. attribute:: transform
  
    The transform that superposes the model onto the reference structure.
    
    :type: :class:`~ost.geom.Mat4`
  
  .. attribute:: gdt_ha
  
    The GDT_HA of the model to the reference structure.

  .. attribute:: gdt_ts

    The GDT_TS of the model to the reference structure.

  """
  def __init__(self, rmsd_common, tm_score, max_sub, 
               gdt_ts, gdt_ha, rmsd_below_five, transform):
    self.rmsd_common=rmsd_common
    self.tm_score=tm_score    
    self.max_sub=max_sub
    self.gdt_ts=gdt_ts
    self.gdt_ha=gdt_ha
    self.rmsd_below_five=rmsd_below_five
    self.transform=transform
    
def _ParseTmScore(lines):
  tf1=[float(i.strip()) for i in lines[23].split()]
  tf2=[float(i.strip()) for i in lines[24].split()]
  tf3=[float(i.strip()) for i in lines[25].split()]
  rot=geom.Mat3(tf1[2], tf1[3], tf1[4], tf2[2], tf2[3],
                  tf2[4], tf3[2], tf3[3], tf3[4])
  tf=geom.Mat4(rot)
  tf.PasteTranslation(geom.Vec3(tf1[1], tf2[1], tf3[1]))
  result=TMScoreResult(float(lines[14].split()[-1].strip()),
                       float(lines[16].split()[2].strip()),
                       float(lines[17].split()[1].strip()),
                       float(lines[18].split()[1].strip()),
                       float(lines[19].split()[1].strip()),
                       float(lines[27].split()[-1].strip()),
                       tf)
  return result

def _RunTmScore(tmscore, tmp_dir):
  model1_filename=os.path.join(tmp_dir, 'model01.pdb')
  model2_filename=os.path.join(tmp_dir, 'model02.pdb')  
  if platform.system() == "Windows":
    tmscore_path=settings.Locate('tmscore.exe', explicit_file_name=tmscore)
    command="\"%s\" %s %s" %(os.path.normpath(tmscore_path), model1_filename, 
                             model2_filename)
  else:
    tmscore_path=settings.Locate('tmscore', explicit_file_name=tmscore)
    command="\"%s\" \"%s\" \"%s\"" % (tmscore_path, model1_filename, 
                                      model2_filename)
  ps=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  stdout,_=ps.communicate()
  lines=stdout.decode().splitlines()
  if (len(lines))<22:
    _CleanupFiles(tmp_dir)
    raise RuntimeError("tmscore superposition failed")
  return _ParseTmScore(lines)


def TMAlign(model1, model2, tmalign=None):
  """
  Performs a sequence independent superposition of model1 onto model2, the 
  reference.
  

  :param model1: The model structure. If the superposition is successful, will 
                 be superposed onto the reference structure
  :type model1: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param model2: The reference structure
  :type model2: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param tmalign: If not None, the path to the tmalign executable.
  :returns: The result of the tmscore superposition
  :rtype: :class:`ost.bindings.TMAlignResult`
  
  :raises: :class:`~ost.settings.FileNotFound` if tmalign could not be located.
  :raises: :class:`RuntimeError` if the superposition failed
  """
  tmp_dir_name=_SetupFiles((model1, model2))
  result=_RunTmAlign(tmalign, tmp_dir_name)
  model1.handle.EditXCS().ApplyTransform(result.transform)
  _CleanupFiles(tmp_dir_name)
  return result

def TMScore(model1, model2, tmscore=None):
  """
  Performs a sequence dependent superposition of model1 onto model2, 
  the reference.

  :param model1: The model structure. If the superposition is successful, will 
                 be superposed onto the reference structure
  :type model1: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param model2: The reference structure
  :type model2: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param tmscore: If not None, the path to the tmscore executable.
  :returns: The result of the tmscore superposition
  :rtype: :class:`TMScoreResult`
  
  :raises: :class:`~ost.settings.FileNotFound` if tmalign could not be located.
  :raises: :class:`RuntimeError` if the superposition failed
  """
  tmp_dir_name=_SetupFiles((model1, model2))
  result=_RunTmScore(tmscore, tmp_dir_name)
  model1.handle.EditXCS().ApplyTransform(result.transform)  
  _CleanupFiles(tmp_dir_name)
  return result


def USAlign(model1, model2, usalign=None):
  """
  Performs a sequence independent superposition of model1 onto model2, the 
  reference. Can deal with multimeric complexes and RNA.

  Creates temporary model files on disk and runs USalign with:
  ``USalign model1.pdb model2.pdb -mm 1 -ter 0 -m rotmat.txt``
  
  :param model1: The model structure. If the superposition is successful, will 
                 be superposed onto the reference structure
  :type model1: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param model2: The reference structure
  :type model2: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param usalign: If not None, the path to the USalign executable. Searches
                  for executable with name ``USalign`` in PATH if not given.
  :returns: The result of the superposition
  :rtype: :class:`ost.bindings.MMAlignResult`
  
  :raises: :class:`~ost.settings.FileNotFound` if executable could not be located.
  :raises: :class:`RuntimeError` if the superposition failed
  """
  tmp_dir_name=_SetupFiles((model1, model2))
  result=_RunUSAlign(usalign, tmp_dir_name)
  model1.handle.EditXCS().ApplyTransform(result.transform)
  _CleanupFiles(tmp_dir_name)
  return result
