import sys
import unittest
from ost import *
from ost import settings
from ost.bindings import blast
import re
import tempfile
import os
import shutil

class TestBlastBindings(unittest.TestCase):
  def setUp(self):
    self.query=seq.CreateSequence('A', 'MKPHPWFFGKIPRAKAEEMLSKQRHDGAFLIRESESAPG'+
                                       'DFSLSVKFGNDVQHFKVLRDGAGKYFLWVVKFNSLNELV'+
                                       'DYHRSTSVSRNQQIFLRDIEQVP')
  def testAllowedMatricesAndPenalties(self):
    blast.Blast(self.query, 'testfiles/seqdb', matrix='BLOSUM45', 
                gap_open=13, gap_ext=3)
    blast.Blast(self.query, 'testfiles/seqdb', matrix='BLOSUM80')
    blast.Blast(self.query, 'testfiles/seqdb', matrix='PAM70')
    blast.Blast(self.query, 'testfiles/seqdb', matrix='PAM30', gap_open=7,
                gap_ext=2)
    blast.Blast(self.query, 'testfiles/seqdb', matrix='BLOSUM62')
    self.assertRaises(ValueError, blast.Blast, self.query, 'testfiles/seqdb', 
                      matrix='MUAHA')
    self.assertRaises(RuntimeError, blast.Blast, self.query, 'testfiles/seqdb',
                      blast_location='/I/do/not/exist')

  def testMissingDB(self):
    self.assertRaises(IOError, blast.Blast, self.query, 
                      'testfiles/nonexistentdb')
  
  def testBlastExec(self):

    hits = blast.Blast(self.query, 'testfiles/seqdb')

    expected_perfect_hit_ids = ["1fhsA", "3n84E", "1griA"]
    perfect_hit_ids = []
    
    for h in hits:
      perfect_match = False
      for aln_patch in h.aligned_patches:
        if aln_patch.seqid == 1.0:
          perfect_match = True
      if perfect_match:
        perfect_hit_ids.append(h.identifier)

    expected_perfect_hit_ids.sort()    
    perfect_hit_ids.sort()
    
    self.assertEqual(perfect_hit_ids, expected_perfect_hit_ids)

  def testBlastParseOutput(self):

    with open('testfiles/raw_blastout.txt','r') as f:
      raw_out=f.read()

    parsed_out=blast.ParseBlastOutput(raw_out)

    # numbers and ids can be specific to a certain blast version
    expected_ids=['1griA','1fhsA','3n84E','1mw9X']
    expected_evalues=[4.80893e-59,4.80893e-59,9.06925e-58,2.96523]
    expected_scores=[534,534,523,27]
    expected_bitscores=[210.305,210.305,206.068,15.0086]
    expected_seqid=[1.0,1.0,1.0,0.30769]

    for i in range(4):
      patch=parsed_out[i].aligned_patches[0]
      self.assertEqual(expected_ids[i],parsed_out[i].identifier)
      self.assertAlmostEqual(patch.evalue,expected_evalues[i])
      self.assertAlmostEqual(patch.score,expected_scores[i])
      self.assertAlmostEqual(patch.bit_score,expected_bitscores[i])
      self.assertAlmostEqual(patch.seqid,expected_seqid[i],places=5)

  def testBlastVersion(self):

    self.assertRaises(RuntimeError, blast.Blast, self.query, 'testfiles/seqdb',
                      blast_location='/I/do/not/exist')
    version=blast.BlastVersion()
    self.assertTrue(isinstance(version,str))
    re_v = re.compile('\d+\.\d+\.\d+')
    match=False
    if re_v.match(version):
      match=True
    self.assertEqual(match,True)

  def testBuildDatabase(self):

    tmp_dir_name=tempfile.mkdtemp()
    db_name=os.path.join(tmp_dir_name,'testdb')
    blast.CreateDB('testfiles/multiple.fasta',db_name)
    self.assertRaises(IOError,blast.CreateDB, 'testfiles/multiple.fasta',db_name,
                      mkdb_cmd='I/am/not/a/damn/executable')

    test_seq=seq.CreateSequence('B','ALRLVKDGFAVAIADYNDATATAVAAEINQAGGRAVAIKVDV'+
                                    'SRRDQVFAAVEQARKALGGFNVIVNNAGIAPSTPIESIT')
    blastout=blast.Blast(test_seq,db_name)
    self.assertTrue(len(blastout)>10)
    shutil.rmtree(tmp_dir_name)

if __name__ == "__main__":
  # test if blast package is available on system, otherwise ignore tests
  try:
    blastpath=settings.Locate(('blastp','blastall'))
  except(settings.FileNotFound):
    print("Could not find blast executable: ignoring unit tests")
    sys.exit(0)
  from ost import testutils
  testutils.RunTests()
