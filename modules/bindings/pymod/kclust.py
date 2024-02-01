import subprocess, os, tempfile
from ost import io, seq, settings

"""
Wrapper for kClust a protein sequence clustering algorithm

unpublished, but mentioned in:

Michael Remmert, Andreas Biegert, Andreas Hauser & Johannes Soeding
HHblits: lighting-fast iterative protein sequence searching by
HMM-HMM alignment
Nature Methods 9, 173-175 (2012)

source code can be downloaded from:

ftp://toolkit.genzentrum.lmu.de/pub/
"""

class cluster:

  """
  Holds the information of one cluster
  
  .. attribute:: sequences
    
    SequenceList containing all sequences of the cluster
  
  .. attribute:: representative_id
  
    a string, that contains the name of the representative sequence of the cluster
    as estimated by kClust. 
    
  .. attribute:: alignment

    alignment handle containing a multiple sequence alignment of all sequences of
    the cluster. Gets only calculated when binding is called with appropriate flag.
  """

  def __init__(self, sequences, representative_id):

    self.sequences=sequences
    self.representative_id=representative_id
    self.alignment=None

def _SetupFiles(sequences):

  # create temporary directory
  tmp_dir_name=tempfile.mkdtemp()
  io.SaveSequenceList(sequences, os.path.join(tmp_dir_name,'fastadb.fasta'))
  return tmp_dir_name

def _CleanupFiles(tmp_dir_name):

  import shutil
  shutil.rmtree(tmp_dir_name)

def _ParseOutput(tmp_dir_name):

  with open(os.path.join(tmp_dir_name,'headers.dmp'),'r') as f:
    header_data=f.readlines()
  with open(os.path.join(tmp_dir_name,'clusters.dmp'),'r') as f:
    cluster_data=f.readlines()
  sequences=io.LoadSequenceList(os.path.join(tmp_dir_name,'fastadb.fasta'))

  clusters=dict()
  header_mapper=dict()
  for line in header_data:
    header_mapper[int(line.split()[0])]=line.split()[1].strip().strip('>')

  #find numeric ids of the representatives of the clusters
  unique_representatives=list()
  for line in cluster_data[1:]:
    actual_cluster=int(line.split()[1])
    try:
      unique_representatives.index(actual_cluster)
    except:
      unique_representatives.append(actual_cluster)

  #assign every header to its corresponding cluster, where the
  #cluster id is given by the id of the representative of the cluster
  for idx in unique_representatives:
    clusters[idx]=seq.CreateSequenceList()
  for line in cluster_data[1:]:
    clusters[int(line.split()[1])].AddSequence(sequences.FindSequence(header_mapper[int(line.split()[0])]))

  #translate into final output

  res=list()
  for k, v in clusters.items():
    res.append(cluster(v, header_mapper[k]))

  return res
  
def _RunkClust(tmp_dir_name, clustering_thresh, create_alignments):

  bitscore=clustering_thresh*0.060269-0.68498

  executable=settings.Locate('kClust') 

  cmd=[]
  cmd.append(executable)
  cmd.append('-i')
  cmd.append(os.path.join(tmp_dir_name,'fastadb.fasta'))
  cmd.append('-d')
  cmd.append(tmp_dir_name)
  cmd.append('-s')
  cmd.append(str(bitscore))

  cmd=' '.join(cmd)
  ps=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  stdout, stderr = ps.communicate()

  result=_ParseOutput(tmp_dir_name)

  if(create_alignments):
    from ost.bindings import clustalw
    for c in result:
      if len(c.sequences)>1:
        c.alignment=clustalw.ClustalW(c.sequences)
      else:
        aln=seq.CreateAlignment()
        aln.AddSequence(c.sequences[0])
        c.alignment=aln

  return result

def kClust(sequences, clustering_thresh=30, create_alignments=False):

  """
  Uses kClust to generate clusters of amino acid sequences.
  
  :param sequences: All sequences you want to cluster. 
  :type sequences: :class:`ost.seq.SequenceList`

  :param clustering_thres: Sequence identity threshold to build clusters. Note,
                           that clustering_thresh is more a rule of thumb, since
                           compositional bias in the sequence can also play a role.
                           The value gets transformed in a bitscore, that is used
                           as an input parameter of kClust 

  :param create_alignments: Flag, wether the alignments of the clusters get calculated.
                            Requires clustalw in the path.

  :returns: A list of cluster instances
  
  :raises: :class:`~ost.settings.FileNotFound` if kClust could not be located.
  """

  tmp_dir_name=_SetupFiles(sequences)
  result=_RunkClust(tmp_dir_name, clustering_thresh, create_alignments)
  _CleanupFiles(tmp_dir_name)
  return result
