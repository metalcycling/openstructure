import unittest
from ost import *
from ost import settings
from ost import testutils
from ost.seq.alg import SequenceIdentity
from ost.bindings import tmtools
from ost.bindings import WrappedTMAlign

class TestTMBindings(unittest.TestCase):
  
  def setUp(self):
    self.protein = io.LoadEntity("testfiles/testprotein.pdb")


  def testTMAlign(self):

    try:
      cad_calc_path = settings.Locate('tmalign')  
    except:
      print("Could not find tmalign executable: ignoring unit tests")
      return

    tm_result = tmtools.TMAlign(self.protein, self.protein)

    # model and reference are the same, we expect pretty good results
    self.assertEqual(tm_result.rmsd, 0.0)
    self.assertEqual(tm_result.tm_score, 1.0)
    self.assertEqual(tm_result.aligned_length, len(self.protein.chains[0].residues))
    self.assertEqual(SequenceIdentity(tm_result.alignment), 100.0)

    # transformation should be identity matrix (no transformation at all...)
    identity = geom.Mat4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1)
    self.assertEqual(tm_result.transform, identity)


  def testTMScore(self):

    try:
      cad_calc_path = settings.Locate('tmscore')  
    except:
      print("Could not find tmalign executable: ignoring unit tests")
      return

    tm_result = tmtools.TMScore(self.protein, self.protein)

    # model and reference are the same, we expect pretty good results
    self.assertEqual(tm_result.rmsd_common, 0.0)
    self.assertEqual(tm_result.tm_score, 1.0)
    self.assertEqual(tm_result.max_sub, 1.0)
    self.assertEqual(tm_result.gdt_ts, 1.0)
    self.assertEqual(tm_result.gdt_ha, 1.0)
    self.assertEqual(tm_result.rmsd_below_five, 0.0)

    # transformation should be identity matrix (no transformation at all...)
    identity = geom.Mat4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1)
    self.assertEqual(tm_result.transform, identity)

  def testWrappedTMAlign(self):

    tm_result = WrappedTMAlign(self.protein.CreateFullView().chains[0],
                               self.protein.CreateFullView().chains[0])

    # model and reference are the same, we expect pretty good results
    self.assertAlmostEqual(tm_result.rmsd, 0.0, places=4)
    self.assertAlmostEqual(tm_result.tm_score, 1.0, places=4)
    self.assertEqual(tm_result.aligned_length, len(self.protein.chains[0].residues))
    self.assertEqual(SequenceIdentity(tm_result.alignment), 100.0)

    # transformation should be identity matrix (no transformation at all...)
    identity = geom.Mat4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1)
    for i in range(4):
      for j in range(4):
        self.assertAlmostEqual(tm_result.transform[i,j], identity[i,j])



if __name__ == "__main__":
  testutils.RunTests()
