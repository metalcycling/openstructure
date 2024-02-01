import unittest
from ost import *
from ost import settings
from ost import seq
from ost.bindings.clustalw import *

class TestSemiGlobalAlign(unittest.TestCase):
  def testSemiGlobalAlignment(self):
    seq_a=seq.CreateSequence('A', 'abcdefghijklmnok')
    seq_b=seq.CreateSequence('B', 'cdehijk')
    alns=seq.alg.SemiGlobalAlign(seq_a, seq_b, seq.alg.BLOSUM62)
    self.assertEqual(len(alns), 1)
    self.assertEqual(alns[0].sequences[0].name, 'A')
    self.assertEqual(alns[0].sequences[1].name, 'B')
    self.assertEqual(str(alns[0].sequences[0]), 'abcdefghijklmnok')
    self.assertEqual(str(alns[0].sequences[1]), '--cde--hijk-----')
    self.assertEqual(alns[0].sequences[0].offset, 0)
    self.assertEqual(alns[0].sequences[1].offset, 0)


if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
