import unittest
from ost import *
from ost import settings
from ost import testutils
from ost.seq.alg import SequenceIdentity
from ost.bindings import ialign

class TestIAlign(unittest.TestCase):
  
  def setUp(self):
    self.protein = io.LoadEntity("testfiles/testprotein.pdb")


  def testIAlign(self):

    try:
      ialign_exec = settings.Locate('ialign.pl')  
    except:
      print("Could not find ialign master perl script: ignoring unit tests")
      return

    ialign_result = ialign.iAlign(self.protein, self.protein)

    # model and reference are the same, we expect pretty good results
    self.assertEqual(ialign_result.rmsd, 0.0)
    self.assertEqual(ialign_result.is_score, 1.0)
    self.assertEqual(SequenceIdentity(ialign_result.alignment), 100.0)
    self.assertEqual(ialign_result.aligned_residues, 78)
    self.assertEqual(ialign_result.aligned_contacts, 91)

    # transformation should be identity matrix (no transformation at all...)
    identity = geom.Mat4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1)
    self.assertEqual(ialign_result.transform, identity)


if __name__ == "__main__":
  testutils.RunTests()
