import unittest
from ost import *
from ost import settings
from ost.bindings import lga
from ost import testutils

class TestLGABindings(unittest.TestCase):
  
  def setUp(self):
    self.protein = io.LoadEntity("testfiles/testprotein.pdb")
    self.chain_a = self.protein.Select("cname=A")

  def testLGA(self):

    try:
      lga_path = settings.Locate('lga')  
    except:
      print("Could not find lga executable: ignoring unit tests")
      return

    lga_result = lga.GDT(self.chain_a, self.chain_a, 
                         reference_length=len(self.chain_a.residues))

    expected_transform = geom.Mat4(1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1)
    self.assertEqual(lga_result.gdt_ts, 100.0)
    self.assertEqual(lga_result.gdt_ha, 100.0)
    self.assertEqual(lga_result.GetTransform(), expected_transform)

if __name__ == "__main__":
  testutils.RunTests()
