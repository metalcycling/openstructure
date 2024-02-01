import unittest
import os
import ost
from ost import seq, io, mol

class TestAlignmentFromChainView(unittest.TestCase):

  def testAlignment(self):

    prot = io.LoadPDB(os.path.join("testfiles", "dummy.pdb"))
    view_resnums = [2,3,4,7,8]
    query = "rnum=" + ','.join([str(item) for item in view_resnums])
    prot_view = prot.Select(query)

    exp_handle_seq = ''.join([r.one_letter_code for r in prot.residues])

    exp_view_seq = ['-'] * len(prot.residues)
    for rnum in view_resnums:
      exp_view_seq[rnum-1] = prot.residues[rnum-1].one_letter_code
    exp_view_seq = ''.join(exp_view_seq)

    aln = seq.alg.AlignmentFromChainView(prot_view.chains[0], 
                                         handle_seq_name='custom_handle_name',
                                         view_seq_name='custom_view_name')

    handle_seq = aln.GetSequence(0)
    view_seq = aln.GetSequence(1)

    self.assertEqual(handle_seq.GetName(), 'custom_handle_name')
    self.assertEqual(view_seq.GetName(), 'custom_view_name')
    self.assertEqual(exp_handle_seq, handle_seq.GetString())
    self.assertEqual(exp_view_seq, view_seq.GetString())

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
