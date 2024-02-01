import unittest
import sys
from ost import *
from ost import settings
from ost.bindings import hbplus

class TestHBPlusBinding(unittest.TestCase):

  def setUp(self):
    self.protein = io.LoadEntity("testfiles/testprotein.pdb")

  def testHBondList(self):
    hbond_list = hbplus.HBondList(self.protein)
    self.assertEqual(len(hbond_list), 499)

  def testHBondScore(self):
    self.assertEqual(hbplus.HBondScore(self.protein, self.protein), 1.0)

if __name__ == "__main__":
  try:
    settings.Locate("hbplus")
  except:
    print("Could not find hbplus, could not test binding...")
    sys.exit(0)
  from ost import testutils
  testutils.RunTests()
