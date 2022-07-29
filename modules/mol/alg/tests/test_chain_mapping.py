import unittest, os, sys
import ost
from ost import conop
from ost import io, mol, seq, settings
# check if we can import: fails if numpy or scipy not available
try:
  from ost.mol.alg.chain_mapping import *
except ImportError:
  print("Failed to import chain_mapping.py. Happens when numpy or scipy "\
        "missing. Ignoring test_chain_mapping.py tests.")
  sys.exit(0)

def _LoadFile(file_name):
  """Helper to avoid repeating input path over and over."""
  return io.LoadPDB(os.path.join('testfiles', file_name))

class TestChainMapper(unittest.TestCase):

  def test_chem_grouping(self):

    ent = _LoadFile("3l1p.1.pdb")
    mapper = ChainMapper(ent)

    # manually extract polypeptide and nucleotide views/sequences
    pep_view_one = ent.Select("cname=A")
    pep_view_two = ent.Select("cname=B")
    nuc_view_one = ent.Select("cname=C")
    nuc_view_two = ent.Select("cname=D")

    pep_s_one = ''.join([r.one_letter_code for r in pep_view_one.residues])
    pep_s_one = seq.CreateSequence("A", pep_s_one)
    pep_s_two = ''.join([r.one_letter_code for r in pep_view_two.residues])
    pep_s_two = seq.CreateSequence("B", pep_s_two)

    nuc_s_one = ''.join([r.one_letter_code for r in nuc_view_one.residues])
    nuc_s_one = seq.CreateSequence("C", nuc_s_one)
    nuc_s_two = ''.join([r.one_letter_code for r in nuc_view_two.residues])
    nuc_s_two = seq.CreateSequence("D", nuc_s_two)

    self.assertEqual(len(mapper.polypep_seqs), 2)
    self.assertEqual(len(mapper.polynuc_seqs), 2)
    self.assertEqual(mapper.polypep_seqs[0].GetName(), pep_s_one.GetName())
    self.assertEqual(mapper.polypep_seqs[1].GetName(), pep_s_two.GetName())
    self.assertEqual(mapper.polynuc_seqs[0].GetName(), nuc_s_one.GetName())
    self.assertEqual(mapper.polynuc_seqs[1].GetName(), nuc_s_two.GetName())
    self.assertEqual(str(mapper.polypep_seqs[0]), str(pep_s_one))
    self.assertEqual(str(mapper.polypep_seqs[1]), str(pep_s_two))
    self.assertEqual(str(mapper.polynuc_seqs[0]), str(nuc_s_one))
    self.assertEqual(str(mapper.polynuc_seqs[1]), str(nuc_s_two))

    # peptide sequences should be in the same group, the nucleotides not
    self.assertEqual(len(mapper.chem_group_alignments), 3)
    self.assertEqual(len(mapper.chem_groups), 3)
    self.assertEqual(mapper.chem_groups[0], ["A", "B"])    
    self.assertEqual(mapper.chem_groups[1], ["C"])    
    self.assertEqual(mapper.chem_groups[2], ["D"])

    # check chem_group_types attribute
    self.assertEqual(mapper.chem_group_types, [ost.mol.ChemType.AMINOACIDS,
                                               ost.mol.ChemType.NUCLEOTIDES,
                                               ost.mol.ChemType.NUCLEOTIDES])

    # check chem_group_ref_seqs attribute
    self.assertEqual(len(mapper.chem_group_ref_seqs), 3)
    self.assertEqual(str(mapper.chem_group_ref_seqs[0]), str(pep_s_one))
    self.assertEqual(str(mapper.chem_group_ref_seqs[1]), str(nuc_s_one))
    self.assertEqual(str(mapper.chem_group_ref_seqs[2]), str(nuc_s_two))

    # check chem_group_alignments attribute
    self.assertEqual(len(mapper.chem_group_alignments), 3)
    self.assertEqual(mapper.chem_group_alignments[0].GetCount(), 2)
    self.assertEqual(mapper.chem_group_alignments[1].GetCount(), 1)
    self.assertEqual(mapper.chem_group_alignments[2].GetCount(), 1)
    s0 = mapper.chem_group_alignments[0].GetSequence(0)
    s1 = mapper.chem_group_alignments[0].GetSequence(1)
    self.assertEqual(s0.GetGaplessString(), str(pep_s_one))
    self.assertEqual(s1.GetGaplessString(), str(pep_s_two))
    s0 = mapper.chem_group_alignments[1].GetSequence(0)
    self.assertEqual(s0.GetGaplessString(), str(nuc_s_one))
    s0 = mapper.chem_group_alignments[2].GetSequence(0)
    self.assertEqual(s0.GetGaplessString(), str(nuc_s_two))

    # ensure that error is triggered if there are insertion codes
    tmp_ent = ent.Copy()
    ed = tmp_ent.EditXCS()
    r = tmp_ent.residues[0]
    ed.SetResidueNumber(r, mol.ResNum(r.GetNumber().GetNum(), 'A'))
    self.assertRaises(Exception, ChainMapper, tmp_ent)

    # ensure that error is triggered if resnums are not strictly increasing
    tmp_ent = ent.Copy()
    ed = tmp_ent.EditXCS()
    r = tmp_ent.residues[0]
    ed.SetResidueNumber(r, mol.ResNum(r.GetNumber().GetNum() + 42))
    self.assertRaises(Exception, ChainMapper, tmp_ent)

    # chain B has a missing Valine... set pep_gap_thr to 0.0 should give an
    # additional chem group
    mapper = ChainMapper(ent, pep_gap_thr=0.0)
    self.assertEqual(len(mapper.chem_groups), 4)

    # introduce single point mutation in ent... increasing pep_seqid_thr to 100
    # should give an additional chem group too
    mapper = ChainMapper(ent, pep_seqid_thr=100.)
    self.assertEqual(len(mapper.chem_groups), 3)
    tmp_ent = ent.Copy()
    tmp_ent.residues[42].SetOneLetterCode('X')
    mapper = ChainMapper(tmp_ent) # default values should tolerate mutation
    self.assertEqual(len(mapper.chem_groups), 3)
    mapper = ChainMapper(tmp_ent, pep_seqid_thr=100.) # strict seqid thresh
    self.assertEqual(len(mapper.chem_groups), 4)

    # Introduce a copy of one of the nucleotide chains so we can play the same
    # game, i.e. delete one nucleotide and mutate one nucleotide in the newly
    # added chain
    tmp_ent = ent.Copy()
    ed = tmp_ent.EditXCS()
    ed.InsertChain("X", tmp_ent.chains[3], deep=True)
    mapper = ChainMapper(tmp_ent)
    self.assertEqual(len(mapper.polynuc_seqs), 3)
    self.assertEqual(len(mapper.chem_groups), 3)
    ed.DeleteResidue(tmp_ent.chains[-1].residues[4])
    tmp_ent.chains[-1].residues[3].SetOneLetterCode('X')
    mapper = ChainMapper(tmp_ent)
    self.assertEqual(len(mapper.polynuc_seqs), 3)
    self.assertEqual(len(mapper.chem_groups), 3)
    mapper = ChainMapper(tmp_ent, nuc_seqid_thr=100.)
    self.assertEqual(len(mapper.polynuc_seqs), 3)
    self.assertEqual(len(mapper.chem_groups), 4)
    mapper = ChainMapper(tmp_ent, nuc_gap_thr=0.0)
    self.assertEqual(len(mapper.polynuc_seqs), 3)
    self.assertEqual(len(mapper.chem_groups), 4)


if __name__ == "__main__":
  from ost import testutils
  if testutils.SetDefaultCompoundLib():
    testutils.RunTests()
  else:
    print('No compound lib available. Ignoring test_chain_mapping.py tests.')
