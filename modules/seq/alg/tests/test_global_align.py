import unittest
from ost import *
from ost import settings
from ost import seq
from ost.bindings.clustalw import *

class TestGlobalAlign(unittest.TestCase):
  def testDeletionInSeqB(self):
    seq_a=seq.CreateSequence('A', 'aacdefghiklmn')
    seq_b=seq.CreateSequence('B', 'acdhiklmn')
    alns=seq.alg.GlobalAlign(seq_a, seq_b, seq.alg.BLOSUM62)
    self.assertEqual(len(alns), 1)
    self.assertEqual(alns[0].sequences[0].name, 'A')
    self.assertEqual(alns[0].sequences[1].name, 'B')
    self.assertEqual(str(alns[0].sequences[0]), 'aacdefghiklmn')
    self.assertEqual(str(alns[0].sequences[1]), '-acd---hiklmn')
    self.assertEqual(alns[0].sequences[0].offset, 0)
    self.assertEqual(alns[0].sequences[1].offset, 0)

  def testDeletionInSeqA(self):
    seq_a=seq.CreateSequence('A', 'acdhiklmn')
    seq_b=seq.CreateSequence('B', 'acdefghiklmn')
    alns=seq.alg.GlobalAlign(seq_a, seq_b, seq.alg.BLOSUM62)
    self.assertEqual(len(alns), 1)
    self.assertEqual(alns[0].sequences[0].name, 'A')
    self.assertEqual(alns[0].sequences[1].name, 'B')

    self.assertEqual(str(alns[0].sequences[0]), 'acd---hiklmn')
    self.assertEqual(str(alns[0].sequences[1]), 'acdefghiklmn')
    self.assertEqual(alns[0].sequences[0].offset, 0)
    self.assertEqual(alns[0].sequences[1].offset, 0)

  def testOffset(self):
    seq_a=seq.CreateSequence('A', 'acdhiklmn')
    seq_b=seq.CreateSequence('B', 'ggiklmn')
    alns=seq.alg.GlobalAlign(seq_a, seq_b, seq.alg.BLOSUM62)
    self.assertEqual(len(alns), 1)
    self.assertEqual(alns[0].sequences[0].name, 'A')
    self.assertEqual(alns[0].sequences[1].name, 'B')

    self.assertEqual(str(alns[0].sequences[0]), 'acdhiklmn')
    self.assertEqual(str(alns[0].sequences[1]), 'g--giklmn')
    self.assertEqual(alns[0].sequences[0].offset, 0)
    self.assertEqual(alns[0].sequences[1].offset, 0)

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
