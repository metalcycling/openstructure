import unittest
import ost
from ost import seq, io, mol

class TestAlignToSeqRes(unittest.TestCase):
  def testAlignWorking(self):
    ent, seqres = io.LoadMMCIF("testfiles/align_to_seqres.mmcif", seqres = True)
    chain =  ent.FindChain("A")
    sequence = seqres.FindSequence(chain.GetName());
    seqres_aln = seq.alg.AlignToSEQRES(chain, sequence, try_resnum_first=False)
    self.assertEqual(str(sequence), "MYTNSDFVVIKALEDGVNVIGLTRGADTRFHHSEKLDKGEV"+
                     "LIAQFTEHTSAIKVRGKAYIQTRHGVIESEGKK")
    self.assertEqual(str(seqres_aln.sequences[1]), "----SDFVVIKALEDGVNVIGLTR--"+
                     "-TRFHHSEKLDKGEVLIAQFTEHTSAIKVRGKAYIQTRHGVIESEGK-")

    seqres_aln = seq.alg.AlignToSEQRES(chain, sequence, try_resnum_first=True)
    self.assertEqual(str(seqres_aln.sequences[1]), "----SDFVVIKALEDGVNVIGLTR--"+
                     "-TRFHHSEKLDKGEVLIAQFTEHTSAIKVRGKAYIQTRHGVIESEGK-")

  def testAlignFail(self):
    ent, seqres = io.LoadMMCIF("testfiles/align_to_seqres_valueerror.mmcif",
                               seqres = True)
    chain =  ent.FindChain("A")
    sequence = seqres.FindSequence(chain.GetName());
    ost.PushVerbosityLevel(0)
    self.assertRaises(ValueError, seq.alg.AlignToSEQRES, chain, sequence, True)
    ost.PopVerbosityLevel()

  def testValidateWorking(self):
    ent, seqres = io.LoadMMCIF("testfiles/align_to_seqres.mmcif", seqres = True)
    chain =  ent.FindChain("A")
    sequence = seqres.FindSequence(chain.GetName());
    seqres_aln = seq.alg.AlignToSEQRES(chain, sequence, try_resnum_first=False)
    self.assertEqual(seq.alg.ValidateSEQRESAlignment(seqres_aln, chain), True)

  def testValidateWorkingOnAttachedView(self):
    ent, seqres = io.LoadMMCIF("testfiles/align_to_seqres.mmcif", seqres = True)
    chain =  ent.FindChain("A")
    sequence = seqres.FindSequence(chain.GetName());
    seqres_aln = seq.alg.AlignToSEQRES(chain, sequence, try_resnum_first=False)
    seqres_aln.AttachView(1, chain.Select(''))
    self.assertEqual(seq.alg.ValidateSEQRESAlignment(seqres_aln), True)

  def testValidateEmptySequenceWorking(self):
    alignment = seq.CreateAlignment(seq.CreateSequence('SEQRES', ''), 
                                    seq.CreateSequence('atoms', ''))
    chain = mol.ChainHandle()
    self.assertEqual(seq.alg.ValidateSEQRESAlignment(alignment, chain), True)

  def testValidateStrandBreakageFail(self):
    ent, seqres = io.LoadMMCIF("testfiles/validate_segres_aln_breakage.mmcif",
                               seqres = True)
    chain = ent.FindChain("A")
    sequence = seqres.FindSequence(chain.GetName());
    seqres_aln = seq.alg.AlignToSEQRES(chain, sequence, try_resnum_first = True,
                                       validate = False)
    self.assertEqual(seq.alg.ValidateSEQRESAlignment(seqres_aln, chain), False)

  def testValidateGapConnectedFail(self):
    ent, seqres = io.LoadMMCIF("testfiles/validate_seqres_aln_connected.mmcif",
                               seqres = True,)
    chain = ent.FindChain("A")
    sequence = seqres.FindSequence(chain.GetName());
    seqres_aln = seq.alg.AlignToSEQRES(chain, sequence, try_resnum_first = True,
                                       validate = False)
    self.assertEqual(seq.alg.ValidateSEQRESAlignment(seqres_aln, chain), False)


if __name__ == "__main__":
  from ost import testutils
  if testutils.DefaultCompoundLibIsSet():
    testutils.RunTests()
  else:
    print('No compound library available. Ignoring test_aligntoseqres.py tests.')
