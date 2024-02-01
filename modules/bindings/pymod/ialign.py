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
Wrappers for the iAlign utility.

References:

Mu Gao and Jeffrey Skolnick, 2010. iAlign: a method for the structural comparison
of protein-protein interfaces. Bioinformatics, 26(18):2259-65.
Mu Gao and Jeffrey Skolnick, 2010. Structural space of protein-protein interfaces
is degenerate, close to complete, and highly connected. PNAS 107(52):22517-22.

Authors: Pascal Benkert, Marco Biasini, Martino Bertoni
"""

import subprocess, os, tempfile, platform
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

class iAlignResult:
  """
  Holds the result of running iAlign
  
  .. attribute:: rmsd
    
    The RMSD of the common Calpha atoms of both structures
  
  .. attribute:: transform
  
    The transform that superposes the model onto the reference structure.
    
    :type: :class:`~ost.geom.Mat4`
  
  .. attribute:: alignment
  
    The alignment of the structures, that is the pairing of Calphas of both 
    structures. Since the programs only read ATOM records, residues consisting 
    of HETATMs (MSE) are not included in the alignment.
    
    :type: :class:`~ost.seq.AlignmentHandle`

  .. attribute:: IS_score

    The IS-score of the structural superposition

  .. attribute:: aligned_residues

    The total number of aligned residues

  .. attribute:: aligned_contacts

    The total number of aligned contacts

  """
  def __init__(self, rmsd, transform, alignment, is_score, 
               aligned_residues, aligned_contacts):
    
    self.rmsd=rmsd
    self.transform=transform
    self.alignment=alignment
    self.is_score=is_score    
    self.aligned_residues=aligned_residues
    self.aligned_contacts=aligned_contacts

def _ParseiAlign(lines):
  info_line=lines[18].split(',')
  is_score=float(info_line[0].split('=')[1].strip())
  aln_residues=int(lines[19].split('=')[1].strip())
  aln_contacts=int(lines[20].split('=')[1].strip())
  info_line=lines[21].split(',')
  rmsd=float(info_line[0].split('=')[1].strip())

  tf1=[float(i.strip()) for i in lines[25][1:].split()]
  tf2=[float(i.strip()) for i in lines[26][1:].split()]
  tf3=[float(i.strip()) for i in lines[27][1:].split()]
  rot=geom.Mat3(tf1[2], tf1[3], tf1[4], tf2[2], tf2[3],
                tf2[4], tf3[2], tf3[3], tf3[4])
  tf=geom.Mat4(rot)
  tf.PasteTranslation(geom.Vec3(tf1[1], tf2[1], tf3[1]))
  seq1 = seq.CreateSequence("1",lines[32].strip())
  seq2 = seq.CreateSequence("2",lines[34].strip())
  alignment = seq.CreateAlignment()
  alignment.AddSequence(seq2)
  alignment.AddSequence(seq1)
  return iAlignResult(rmsd, tf, alignment, is_score, aln_residues, aln_contacts)

def _RuniAlign(ialign, tmp_dir, options={}):
  opts = {'a' : 1,  # concise output
          'w' : tmp_dir
         } 
  opts.update(options)
  cmd_opts = []
  for k, v in opts.items():
    if type(v) == type(True):
      if v == True:
        cmd_opts.append('-%s' % str(k))
    else:
      cmd_opts.append('-%s %s' % (str(k), str(v)))
  cmd_opts = ' '.join(cmd_opts)
  model1_filename=os.path.join(tmp_dir, 'model01.pdb')
  model2_filename=os.path.join(tmp_dir, 'model02.pdb')
  if platform.system() == "Windows":
    ialign_path=settings.Locate('ialign.pl', explicit_file_name=ialign)
    command="\"%s\" %s %s %s" % (os.path.normpath(ialign_path), model1_filename, model2_filename, cmd_opts)
  else:
    ialign_path=settings.Locate('ialign.pl', explicit_file_name=ialign)  
    command="\"%s\" \"%s\" \"%s\" %s" % (ialign_path, model1_filename, model2_filename, cmd_opts)

  ps=subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
  stdout,_=ps.communicate()
  lines=stdout.decode().splitlines()

  if (len(lines))<22:
    _CleanupFiles(tmp_dir)
    #for l in lines:
    #  print l
    raise RuntimeError("iAlign superposition failed")
  return _ParseiAlign(lines)


def iAlign(model1, model2, ialign=None):
  """
  Compare protein-protein interfaces of the structures of two pairs of
  protein complexes and suporpose them.
  

  :param model1: The model structure. If the superposition is successful, will 
                 be superposed onto the reference structure
  :type model1: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param model2: The reference structure
  :type model2: :class:`~ost.mol.EntityView` or :class:`~ost.mol.EntityHandle`
  :param ialign: If not None, the path to the ialign executable.
  :returns: The result of the tmscore superposition
  :rtype: :class:`iAlignResult`
  
  :raises: :class:`~ost.settings.FileNotFound` if tmalign could not be located.
  :raises: :class:`RuntimeError` if the superposition failed
  """
  tmp_dir_name=_SetupFiles((model1, model2))
  result=_RuniAlign(ialign, tmp_dir_name)
  model1.handle.EditXCS().ApplyTransform(result.transform)
  _CleanupFiles(tmp_dir_name)
  return result
