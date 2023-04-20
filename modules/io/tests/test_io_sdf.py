import unittest
from ost import *
import subprocess

class TestSDF(unittest.TestCase):
  def setUp(self):
    pass

  def test_LoadEntity(self):
    ent = io.LoadSDF('testfiles/sdf/compound.sdf')
    self.assertEqual(len(ent.chains), 4)
    self.assertEqual(len(ent.atoms), 180)
    self.assertEqual(len(ent.bonds), 188)

  def test_LoadEntity_crlf(self):
    ent = io.LoadSDF('testfiles/sdf/6d5w_rank1_crlf.sdf.gz')
    self.assertEqual(len(ent.atoms), 21)
    self.assertEqual(len(ent.bonds), 24)
    
if __name__== '__main__':
  from ost import testutils
  testutils.RunTests()


 
