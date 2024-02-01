import sys
import unittest
from ost import *
from ost import settings
from ost.bindings import kclust

class TestkClustBindings(unittest.TestCase):
  
  def testkClustExec(self):

    seq_list=io.LoadSequenceList('testfiles/multiple.fasta')

    clusters40=kclust.kClust(seq_list,clustering_thresh=40)
    clusters90=kclust.kClust(seq_list,clustering_thresh=90)

    self.assertEqual(len(clusters40),16)
    self.assertEqual(len(clusters90),27)

  def testkClustParseOutput(self):

    representative_ids=['Q5KTS5.1','P46331.2','P66776.1','P42317.2','P50199.1',
                        'Q59787.1','Q9WYG0.1','P50197.1','P50198.1','Q7Z4W1.2',
                        'Q04520.1']

    cluster_sizes=[3,3,4,1,7,3,2,1,1,3,1]

    result=kclust._ParseOutput('testfiles')

    self.assertEqual(len(result),11)

    for c, r_id in zip(result, representative_ids):
      self.assertEqual(c.representative_id, r_id)
    for c, s in zip(result, cluster_sizes):
      self.assertEqual(len(c.sequences),s)

if __name__ == "__main__":
  # test if kClust is available on system, otherwise ignore tests
  try:
    blastpath=settings.Locate(('kClust'))
  except(settings.FileNotFound):
    print("Could not find kClust executable: ignoring unit tests")
    sys.exit(0)
  from ost import testutils
  testutils.RunTests()
