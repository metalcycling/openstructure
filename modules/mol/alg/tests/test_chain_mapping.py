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

def _CompareViews(v1, v2):
  """Iterates over atoms of the two views and compares qualified names / pos
  """
  if len(v1.atoms) != len(v2.atoms):
    return False
  for a,b in zip(v1.atoms, v2.atoms):
    if a.GetQualifiedName() != b.GetQualifiedName():
      return False
    a_p = a.GetPos()
    b_p = b.GetPos()
    d = a_p - b_p
    if max([abs(d[0]), abs(d[1]), abs(d[2])]) > 0.001:
      return False
  return True

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

    for s in mapper.polypep_seqs:
      self.assertTrue(s.HasAttachedView())
    for s in mapper.polynuc_seqs:
      self.assertTrue(s.HasAttachedView())
    self.assertTrue(_CompareViews(mapper.polypep_seqs[0].GetAttachedView(), pep_view_one))
    self.assertTrue(_CompareViews(mapper.polypep_seqs[1].GetAttachedView(), pep_view_two))
    self.assertTrue(_CompareViews(mapper.polynuc_seqs[0].GetAttachedView(), nuc_view_one))
    self.assertTrue(_CompareViews(mapper.polynuc_seqs[1].GetAttachedView(), nuc_view_two))

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
    for s in mapper.chem_group_ref_seqs:
      self.assertTrue(s.HasAttachedView()) 
    self.assertTrue(_CompareViews(mapper.chem_group_ref_seqs[0].GetAttachedView(), pep_view_one))
    self.assertTrue(_CompareViews(mapper.chem_group_ref_seqs[1].GetAttachedView(), nuc_view_one))
    self.assertTrue(_CompareViews(mapper.chem_group_ref_seqs[2].GetAttachedView(), nuc_view_two))

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
    self.assertTrue(mapper.chem_group_alignments[0].GetSequence(0).HasAttachedView())
    self.assertTrue(mapper.chem_group_alignments[0].GetSequence(1).HasAttachedView())
    self.assertTrue(mapper.chem_group_alignments[1].GetSequence(0).HasAttachedView())
    self.assertTrue(mapper.chem_group_alignments[2].GetSequence(0).HasAttachedView())
    self.assertTrue(_CompareViews(mapper.chem_group_alignments[0].GetSequence(0).GetAttachedView(), pep_view_one))
    self.assertTrue(_CompareViews(mapper.chem_group_alignments[0].GetSequence(1).GetAttachedView(), pep_view_two))
    self.assertTrue(_CompareViews(mapper.chem_group_alignments[1].GetSequence(0).GetAttachedView(), nuc_view_one))
    self.assertTrue(_CompareViews(mapper.chem_group_alignments[2].GetSequence(0).GetAttachedView(), nuc_view_two))

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

  def test_chem_mapping(self):

    ref = _LoadFile("3l1p.1.pdb")
    mdl = _LoadFile("3l1p.1_model.pdb")

    # manually extract polypeptide and nucleotide views/sequences
    ref_pep_view_one = ref.Select("cname=A")
    ref_pep_view_two = ref.Select("cname=B")
    ref_nuc_view_one = ref.Select("cname=C")
    ref_nuc_view_two = ref.Select("cname=D")

    ref_pep_s_one = ''.join([r.one_letter_code for r in ref_pep_view_one.residues])
    ref_pep_s_one = seq.CreateSequence("A", ref_pep_s_one)
    ref_pep_s_two = ''.join([r.one_letter_code for r in ref_pep_view_two.residues])
    ref_pep_s_two = seq.CreateSequence("B", ref_pep_s_two)

    ref_nuc_s_one = ''.join([r.one_letter_code for r in ref_nuc_view_one.residues])
    ref_nuc_s_one = seq.CreateSequence("C", ref_nuc_s_one)
    ref_nuc_s_two = ''.join([r.one_letter_code for r in ref_nuc_view_two.residues])
    ref_nuc_s_two = seq.CreateSequence("D", ref_nuc_s_two)


    mdl_pep_view_one = mdl.Select("cname=X")
    mdl_pep_view_two = mdl.Select("cname=Y")
    mdl_nuc_view_one = mdl.Select("cname=Z")

    mdl_pep_s_one = ''.join([r.one_letter_code for r in mdl_pep_view_one.residues])
    mdl_pep_s_one = seq.CreateSequence("X", mdl_pep_s_one)
    mdl_pep_s_two = ''.join([r.one_letter_code for r in mdl_pep_view_two.residues])
    mdl_pep_s_two = seq.CreateSequence("Y", mdl_pep_s_two)

    mdl_nuc_s_one = ''.join([r.one_letter_code for r in mdl_nuc_view_one.residues])
    mdl_nuc_s_one = seq.CreateSequence("Z", mdl_nuc_s_one)

    mapper = ChainMapper(ref)

    chem_mapping, alns, mdl_view = mapper.GetChemMapping(mdl)

    self.assertEqual(len(mapper.chem_groups), 3)
    self.assertEqual(len(chem_mapping), len(mapper.chem_groups))
    self.assertEqual(chem_mapping[0], ['X', 'Y'])
    self.assertEqual(chem_mapping[1], [])
    self.assertEqual(chem_mapping[2], ['Z'])

    self.assertEqual(len(alns), 3)
    self.assertEqual(len(alns[0]), 2)
    self.assertEqual(len(alns[1]), 0)
    self.assertEqual(len(alns[2]), 1)

    self.assertEqual(alns[0][0].GetSequence(0).GetGaplessString(), str(ref_pep_s_one))
    self.assertEqual(alns[0][0].GetSequence(1).GetGaplessString(), str(mdl_pep_s_one))
    self.assertEqual(alns[0][1].GetSequence(0).GetGaplessString(), str(ref_pep_s_one))
    self.assertEqual(alns[0][1].GetSequence(1).GetGaplessString(), str(mdl_pep_s_two))
    self.assertEqual(alns[2][0].GetSequence(0).GetGaplessString(), str(ref_nuc_s_two))
    self.assertEqual(alns[2][0].GetSequence(1).GetGaplessString(), str(mdl_nuc_s_one))

    self.assertTrue(alns[0][0].GetSequence(0).HasAttachedView())
    self.assertTrue(alns[0][0].GetSequence(1).HasAttachedView())
    self.assertTrue(alns[0][1].GetSequence(0).HasAttachedView())
    self.assertTrue(alns[0][1].GetSequence(1).HasAttachedView())
    self.assertTrue(alns[2][0].GetSequence(0).HasAttachedView())
    self.assertTrue(alns[2][0].GetSequence(1).HasAttachedView())
    self.assertTrue(_CompareViews(alns[0][0].GetSequence(0).GetAttachedView(),ref_pep_view_one))
    self.assertTrue(_CompareViews(alns[0][0].GetSequence(1).GetAttachedView(),mdl_pep_view_one))
    self.assertTrue(_CompareViews(alns[0][1].GetSequence(0).GetAttachedView(),ref_pep_view_one))
    self.assertTrue(_CompareViews(alns[0][1].GetSequence(1).GetAttachedView(),mdl_pep_view_two))
    self.assertTrue(_CompareViews(alns[2][0].GetSequence(0).GetAttachedView(),ref_nuc_view_two))
    self.assertTrue(_CompareViews(alns[2][0].GetSequence(1).GetAttachedView(),mdl_nuc_view_one))

  def test_chain_mapping(self):
    ref = _LoadFile("3l1p.1.pdb")
    mdl = _LoadFile("3l1p.1_model.pdb")
    mapper = ChainMapper(ref)

    # This is not supposed to be in depth algorithm testing, we just check
    # whether the various algorithms return sensible chain mappings

    # lDDT based chain mappings
    naive_lddt_res = mapper.GetlDDTMapping(mdl, strategy="naive")
    self.assertEqual(naive_lddt_res.mapping, [['X', 'Y'],[None],['Z']])

    # the "fast" strategy produces actually a suboptimal mapping in this case...
    greedy_lddt_res = mapper.GetlDDTMapping(mdl, strategy="greedy_fast")
    self.assertEqual(greedy_lddt_res.mapping, [['Y', 'X'],[None],['Z']])

    greedy_lddt_res = mapper.GetlDDTMapping(mdl, strategy="greedy_full")
    self.assertEqual(greedy_lddt_res.mapping, [['X', 'Y'],[None],['Z']])

    greedy_lddt_res = mapper.GetlDDTMapping(mdl, strategy="greedy_block")
    self.assertEqual(greedy_lddt_res.mapping, [['X', 'Y'],[None],['Z']])


    # QS score based chain mappings
    naive_qsscore_res = mapper.GetQSScoreMapping(mdl, strategy="naive")
    self.assertEqual(naive_qsscore_res.mapping, [['X', 'Y'],[None],['Z']])

    greedy_qsscore_res = mapper.GetQSScoreMapping(mdl, strategy="greedy_fast")
    self.assertEqual(naive_qsscore_res.mapping, [['X', 'Y'],[None],['Z']])

    greedy_qsscore_res = mapper.GetQSScoreMapping(mdl, strategy="greedy_full")
    self.assertEqual(naive_qsscore_res.mapping, [['X', 'Y'],[None],['Z']])

    greedy_qsscore_res = mapper.GetQSScoreMapping(mdl, strategy="greedy_block")
    self.assertEqual(naive_qsscore_res.mapping, [['X', 'Y'],[None],['Z']])


    # rigid chain mappings
    greedy_rigid_res = mapper.GetRigidMapping(mdl, strategy="greedy_single_gdtts")
    self.assertEqual(greedy_rigid_res.mapping, [['X', 'Y'],[None],['Z']])

    greedy_rigid_res = mapper.GetRigidMapping(mdl, strategy="greedy_iterative_gdtts")
    self.assertEqual(greedy_rigid_res.mapping, [['X', 'Y'],[None],['Z']])

    greedy_rigid_res = mapper.GetRigidMapping(mdl, strategy="greedy_iterative_rmsd")
    self.assertEqual(greedy_rigid_res.mapping, [['X', 'Y'],[None],['Z']])

    # test flat mapping functionality of MappingResult
    flat_map = greedy_rigid_res.GetFlatMapping()
    self.assertEqual(len(flat_map), 3)
    self.assertEqual(flat_map[greedy_rigid_res.chem_groups[0][0]], 'X')
    self.assertEqual(flat_map[greedy_rigid_res.chem_groups[0][1]], 'Y')
    self.assertEqual(flat_map[greedy_rigid_res.chem_groups[2][0]], 'Z')
    flat_map = greedy_rigid_res.GetFlatMapping(mdl_as_key=True)
    self.assertEqual(len(flat_map), 3)
    self.assertEqual(greedy_rigid_res.chem_groups[0][0], flat_map['X'])
    self.assertEqual(greedy_rigid_res.chem_groups[0][1], flat_map['Y'])
    self.assertEqual(greedy_rigid_res.chem_groups[2][0], flat_map['Z'])


if __name__ == "__main__":
  from ost import testutils
  if testutils.SetDefaultCompoundLib():
    testutils.RunTests()
  else:
    print('No compound lib available. Ignoring test_chain_mapping.py tests.')
