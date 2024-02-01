import unittest
from ost import *
from ost import settings
from ost import seq

class TestWeightMatrix(unittest.TestCase):
  
  def test_GetWeight(self):
    mat = seq.alg.BLOSUM62
    self.assertEqual(mat.GetWeight('A', 'A'), 4)
    self.assertEqual(mat.GetWeight('A', 'B'), -2)
    self.assertEqual(mat.GetWeight('a', 'A'), 4)
    self.assertEqual(mat.GetWeight('A', 'b'), -2)
    self.assertEqual(mat.GetWeight('A', '?'), 0)
    self.assertEqual(mat.GetWeight('E', '?'), 0)
    self.assertEqual(mat.GetWeight('Y', '?'), 0)
    self.assertEqual(mat.GetWeight('?', 'A'), 0)
    self.assertEqual(mat.GetWeight('?', 'E'), 0)
    self.assertEqual(mat.GetWeight('?', 'Y'), 0)
    self.assertEqual(mat.GetWeight('?', 'y'), 0)
    
    
  def test_SetWeight(self):
    mat = seq.alg.BLOSUM62
    self.assertEqual(mat.GetWeight('A', 'A'), 4)
    mat.SetWeight('A', 'A', -1)
    self.assertEqual(mat.GetWeight('A', 'A'), -1)
    
    self.assertEqual(mat.GetWeight('A', 'B'), -2)
    mat.SetWeight('A', 'B', 10)
    self.assertEqual(mat.GetWeight('A', 'B'), 10)
    
    self.assertEqual(mat.GetWeight('E', '?'), 0)
    mat.SetWeight('E', '?', 10)
    self.assertEqual(mat.GetWeight('E', '?'), 0)
    
    self.assertEqual(mat.GetWeight('?', 'E'), 0)
    mat.SetWeight('?', 'E', 10)
    self.assertEqual(mat.GetWeight('?', 'E'), 0)

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
