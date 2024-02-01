import os.path
import string
import os

from ost import io

class ModelRepository:
  """
  Model repository. A model repository abstracts the way that PDB files are 
  loaded. Instead of explicitly specifying the PDB filename, only the PDB 
  id (and optionally a chain) needs to be specified. The actual files are then 
  resolved by the repository.
  
  Usage
  -----
  The usage pattern of the model repository is simple. After construction, 
  models may be loaded by passing in a model id and optionally a number of 
  chain names (see documentation for io.LoadPDB).
  
  Example:
  import string
  repos=repository.ModelRepository('path_to_pdbs', 
                                   file_pattern='pdb%(id)s.ent.gz', 
                                   transform=string.lower)
  # load 1ake (note that the name is transformed by string.lower)
  m=repos.Load('1AKE')
  """
  def __init__(self, directory=None, 
               file_pattern='%(id)s.pdb',transform=str):
    """
    Construct new model repository
    """
    if directory==None:
      self.directory_=os.getenv('PDB_PATH', '')
    else:
      self.directory_=directory;
    self.file_pattern_=file_pattern
    self.transform_=transform or string.__init__
  def FilenameForModel(self, pdb_id, chain):
    pdb_id=self.transform_(pdb_id)
    basename=self.file_pattern_ % {'id' : pdb_id, 'chain' :chain, 'dir' : pdb_id[1:3]}
    return os.path.join(self.directory_, basename)
    
  def Load(self, pdb_id, chains='', calpha_only=False, fault_tolerant=False):
    return io.LoadPDB(self.FilenameForModel(pdb_id, chains),
                      chains, calpha_only=calpha_only, 
                      fault_tolerant=fault_tolerant)
                      
  def LoadMulti(self, pdb_id, chains=""):
    return io.LoadMultiPDB(self.FilenameForModel(pdb_id, chains))
    
