import os
from ost import seq, mol, io
import tempfile
def SaveToTempDir(objects, seq_format='fasta', structure_format='pdb'):
  """
  Take all objects and saves them to a temporary directory. The paths of the 
  saved files will be returned as a tuple. This works for alignments, and 
  structure files. The output format for sequence files and structure files can
  be changed by setting to seq_format and structure_format parameters 
  appropriately.
  """
  # create temporary directory
  tmp_dir_name=tempfile.mkdtemp()
  file_names=[]
  for index, obj in enumerate(objects):
    if isinstance(obj, seq.AlignmentHandle):
      name=os.path.join(tmp_dir_name, 'aln%02d.fasta' % (index+1))
      io.SaveAlignment(obj, name, seq_format)
      file_names.append(name)
      continue
    if isinstance(obj, seq.SequenceHandle):
      name=os.path.join(tmp_dir_name, 'seq%02d.fasta' % (index+1))
      io.SaveSequence(obj, name, seq_format)
      file_names.append(name)
      continue
    if isinstance(obj, seq.ConstSequenceList) or isinstance(obj, seq.SequenceList):
      name=os.path.join(tmp_dir_name, 'sql%02d.fasta' % (index+1))
      io.SaveSequenceList(obj, name, seq_format)
      file_names.append(name)
      continue
    if isinstance(obj, mol.EntityView) or isinstance(obj, mol.EntityHandle):
      name=os.path.join(tmp_dir_name, 'mol%02d.pdb' % (index+1))
      io.SaveEntity(obj, name, structure_format)
      file_names.append(name)
      continue
  return file_names
  
class TempDirWithFiles:
  def __init__(self, objects, seq_format='fasta', structure_format='pdb'):
    self.files=SaveToTempDir(objects, seq_format=seq_format, 
                             structure_format=structure_format)
    self.dirname=os.path.dirname(self.files[0])
  def Cleanup(self):
    import shutil
    shutil.rmtree(self.dirname)


