import sys
import unittest
from ost import *
from ost import settings
from ost.bindings import clustalw

class TestClustalWBindings(unittest.TestCase):
  
  def setUp(self):
    self.protein = io.LoadEntity("testfiles/testprotein.pdb").Select("cname=A")

    self.targetseq = io.LoadSequence("testfiles/test.fasta")
    self.targetseq.AttachView(self.protein)
    self.templseq = io.LoadSequence("testfiles/similar.fasta")
    self.multseq = io.LoadSequenceList("testfiles/multiple.fasta")
    self.pw_alignment = io.LoadAlignment("testfiles/pairwise_aln.fasta")
    self.nopgap_pw_alignment = io.LoadAlignment("testfiles/nopgap_pairwise_aln.fasta")
    self.mult_alignment = io.LoadAlignment("testfiles/multiple_aln.fasta")

    self.strseq1 = self.targetseq.GetGaplessString()
    self.strseq2 = self.templseq.GetGaplessString()

    self.seq1 = io.LoadSequence("testfiles/seq1.fasta")
    self.seq2 = io.LoadSequence("testfiles/seq2.fasta")
    self.seq1_seq2_alignment = io.LoadAlignment("testfiles/seq1_seq2_aln.fasta")
    self.seq1_seq2_alignment_options_changed = io.LoadAlignment("testfiles/seq1_seq2_aln_options_changed.fasta")

  def testPairwiseClustalW(self):
    aln=clustalw.ClustalW(self.targetseq, self.templseq)
    assert self.pw_alignment.ToString(80) == aln.ToString(80), \
           "Pairwise alignment differs from precomputed one"
  
  def testNoPGapPariwiseClustalW(self):
    aln=clustalw.ClustalW(self.targetseq, self.templseq, nopgap=True)
    assert self.nopgap_pw_alignment.ToString(80) == aln.ToString(80), \
           "NoPGap pairwise alignment differs from precomputed one"

  def testAttachedViewClustalW(self):
    aln=clustalw.ClustalW(self.targetseq, self.templseq)
    assert aln.FindSequence("testseq").HasAttachedView(), \
           "Aligned sequence doesn't have an attached view"
    
  def testMultipleClustalW(self):
    aln=clustalw.ClustalW(self.multseq)
    # order of sequences in aln may differ (reported by Andrius Merkys who
    # found architecture dependent behaviour)
    # We therefore check for matching sequence names and respective sequences
    # but not for order.
    ref_sequences = {s.GetName(): str(s) for s in self.mult_alignment.sequences}
    sequences = {s.GetName(): str(s) for s in aln.sequences}
    self.assertEqual(sorted(ref_sequences.keys()), sorted(sequences.keys()))
    for sname, s in ref_sequences.items():
      self.assertEqual(s, sequences[sname])


  def testStringClustalW(self):
    aln=clustalw.ClustalW(self.strseq1, self.strseq2)
    aln.SetSequenceName(0,self.targetseq.GetName())
    aln.SetSequenceName(1,self.templseq.GetName())
    assert self.pw_alignment.ToString(80) == aln.ToString(80), \
           "Pairwise alignment using two strings differs from precomputed one \n%s \n%s" \
           %(self.pw_alignment.ToString(80),aln.ToString(80))

  def testPairwiseClustalWChangedOptions(self):
    # five residues removed two positions before the end of seq2
    aln=clustalw.ClustalW(self.seq1,self.seq2)
    assert self.seq1_seq2_alignment.ToString(80) == aln.ToString(80), \
           "Pairwise alignment with default gap penalties differs from precomputed one"
    aln=clustalw.ClustalW(self.seq1,self.seq2,clustalw_option_string="-GAPOPEN=2 -GAPEXT=0")
    assert self.seq1_seq2_alignment_options_changed.ToString(80) == aln.ToString(80), \
           "Pairwise alignment with modified gap penalties differs from precomputed one"

  def testUniqueIdentifier(self):
    # common case
    seq1 = seq.CreateSequence('heelloo', 'AWESOME')
    seq2 = seq.CreateSequence('heelloo', 'AWESOME')
    self.assertRaises(ValueError, clustalw.ClustalW, seq1, seq2)
    # nasty case with spaces
    seq2 = seq.CreateSequence('heelloo dear', 'AWESOME')
    self.assertRaises(ValueError, clustalw.ClustalW, seq1, seq2)

if __name__ == "__main__":
  # test if clustalw package is available on system, otherwise ignore tests
  try:
    clustalw_path=settings.Locate(('clustalw', 'clustalw2'))
  except(settings.FileNotFound):
    print("Could not find clustalw executable: ignoring unit tests")
    sys.exit(0)
  from ost import testutils
  testutils.RunTests()
