import sys
import unittest
from ost.seq.alg import aaindex

class TestAAIndex(unittest.TestCase):
  
    
  def testStuff(self):
      index = aaindex.AAIndex()

      # test some random aaindex annotations from all three files
      # the files themselves are lazy loaded until the required
      # annotation is found

      # entry from aaindex1
      self.assertEqual(index["BHAR880101"].GetScore("C"), 0.346)

      # entries from aaindex2
      # symmetric one:
      self.assertEqual(index["HENS920102"].GetPairScore("C", "A"), -1.0)
      self.assertEqual(index["HENS920102"].GetPairScore("A", "C"), -1.0)
      # non-symmetric one:
      self.assertEqual(index["LINK010101"].GetPairScore("C", "A"), 0.035)
      self.assertEqual(index["LINK010101"].GetPairScore("A", "C"), 0.000)


      # entries from aaindex3
      # symmetric one
      self.assertEqual(index["TANS760102"].GetPairScore("H", "R"), 17)
      self.assertEqual(index["TANS760102"].GetPairScore("R", "H"), 17)
      # non symmetric one
      self.assertEqual(index["ZHAC000102"].GetPairScore("H", "R"), 2.19)
      self.assertEqual(index["ZHAC000102"].GetPairScore("R", "H"), 0.94)

      # depending on annotation type (single amino acids or pairs), we need
      # to call the right functions
      with self.assertRaises(RuntimeError):
          index["TANS760102"].GetScore("H")

      with self.assertRaises(RuntimeError):
          index["BHAR880101"].GetPairScore("H", "R")


if __name__ == "__main__":
  from ost import testutils
  # the function below indirectly enables GetSharedDataPath when
  # calling stuff from python which is the case in unit tests
  testutils.SetDefaultCompoundLib()
  testutils.RunTests()
