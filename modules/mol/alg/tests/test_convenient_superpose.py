from ost import io,mol,geom
import unittest
import os
import random

class TestConvenientSuperpose(unittest.TestCase):
  
  def setUp(self):
    pass
  
  def assertEqualAtomOrder(self, view1, view2):
    self.assertEqual(len(view1.atoms),len(view2.atoms))
    for a1, a2 in zip(view1.atoms, view2.atoms):
      self.assertEqual(a1.element, a2.element)
      self.assertTrue(geom.Equal(a1.pos, a2.pos))

  def testAssertion(self):
    ent1_ent = io.LoadEntity(os.path.join("testfiles","1aho.pdb"))
    ent1_view = ent1_ent.Select("")
    self.assertEqualAtomOrder(ent1_ent, ent1_ent)
    self.assertEqualAtomOrder(ent1_ent, ent1_view)
    self.assertEqualAtomOrder(ent1_view, ent1_ent)
    self.assertEqualAtomOrder(ent1_view, ent1_view)
    self.assertRaises(AssertionError, self.assertEqualAtomOrder, ent1_view.Select("ele!=H"), ent1_view)
    
  def testCorrectlyOrdered(self):
    ent1_ent = io.LoadEntity(os.path.join("testfiles","1aho.pdb"))
    ent1_view = ent1_ent.Select("")
    ## test MatchResidueByNum
    view1, view2 = mol.alg.MatchResidueByNum(ent1_ent, ent1_ent)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent1_view, ent1_ent)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent1_ent, ent1_view)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent1_view, ent1_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByIdx
    view1, view2 = mol.alg.MatchResidueByIdx(ent1_ent, ent1_ent)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByIdx(ent1_view, ent1_ent)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByIdx(ent1_ent, ent1_view)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByIdx(ent1_view, ent1_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByGlobalAln
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent1_ent, ent1_ent)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent1_view, ent1_ent)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent1_ent, ent1_view)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent1_view, ent1_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByLocalAln
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent1_ent, ent1_ent)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent1_view, ent1_ent)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent1_ent, ent1_view)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent1_view, ent1_view)
    self.assertEqualAtomOrder(view1, view2)

  def testMissingFirstAtom(self):
    ent_view = io.LoadEntity(os.path.join("testfiles","1aho.pdb")).Select("")
    ent_view_missing = ent_view.Select("not (cname=A and rnum=1 and aname=N)")
    ## test MatchResidueByNum
    view1, view2 = mol.alg.MatchResidueByNum(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByIdx
    view1, view2 = mol.alg.MatchResidueByIdx(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByIdx(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByGlobalAln
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByLocalAln
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)

  def testMissingManyAtoms(self):
    ent_view = io.LoadEntity(os.path.join("testfiles","1aho.pdb")).Select("")
    ent_view_missing = ent_view.Select("not (cname=A and rnum=3,19,32 and aname=CB,CA,CD)")
    ## test MatchResidueByNum
    view1, view2 = mol.alg.MatchResidueByNum(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByIdx
    view1, view2 = mol.alg.MatchResidueByIdx(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByIdx(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByGlobalAln
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByLocalAln
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)

  def testMissingFirstResidue(self):
    ent_view = io.LoadEntity(os.path.join("testfiles","1aho.pdb")).Select("")
    ent_view_missing = ent_view.Select("not (cname=A and rnum=1)")
    ## test MatchResidueByNum
    view1, view2 = mol.alg.MatchResidueByNum(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByGlobalAln
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByLocalAln
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)

  def testMissingHydrogens(self):
    ent_view = io.LoadEntity(os.path.join("testfiles","1aho.pdb")).Select("")
    ent_view_missing = ent_view.Select("ele!=H")
    ## test MatchResidueByNum
    view1, view2 = mol.alg.MatchResidueByNum(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByIdx
    view1, view2 = mol.alg.MatchResidueByIdx(ent_view, ent_view_missing)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByIdx(ent_view_missing, ent_view)
    self.assertEqualAtomOrder(view1, view2)
  
  def testWrongAtomOrder(self):
    ent_view = io.LoadEntity(os.path.join("testfiles","1aho.pdb")).Select("")
    ent_view_wrong = ent_view.CreateEmptyView()
    for c in ent_view.chains:
      ent_view_wrong.AddChain(c)
      for r in c.residues:
        ent_view_wrong.AddResidue(r)
        atoms = list(r.atoms)
        random.shuffle(atoms)
        for a in atoms:
          ent_view_wrong.AddAtom(a)
    ## test MatchResidueByNum
    view1, view2 = mol.alg.MatchResidueByNum(ent_view, ent_view_wrong)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent_view_wrong, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByIdx
    view1, view2 = mol.alg.MatchResidueByIdx(ent_view, ent_view_wrong)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByIdx(ent_view_wrong, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByGlobalAln
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_view, ent_view_wrong)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_view_wrong, ent_view)
    self.assertEqualAtomOrder(view1, view2)
    ## test MatchResidueByLocalAln
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_view, ent_view_wrong)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_view_wrong, ent_view)
    self.assertEqualAtomOrder(view1, view2)

  def testWrongResidueOrder(self):
    ent_view = io.LoadEntity(os.path.join("testfiles","1aho.pdb")).Select("")
    ent_view_wrong = ent_view.CreateEmptyView()
    for c in ent_view.chains:
      ent_view_wrong.AddChain(c)
      residues = list(c.residues)
      random.shuffle(residues)
      for r in residues:
        ent_view_wrong.AddResidue(r)
        for a in r.atoms:
          av=ent_view_wrong.AddAtom(a)
    ## test MatchResidueByNum
    view1, view2 = mol.alg.MatchResidueByNum(ent_view, ent_view_wrong)
    self.assertEqualAtomOrder(view1, view2)
    view1, view2 = mol.alg.MatchResidueByNum(ent_view_wrong, ent_view)
    self.assertEqualAtomOrder(view1, view2)

  def testPeptideFilter(self):
    # check that CA only chains are ok and ligand chains are skipped
    ent_1_full = io.LoadPDB(os.path.join("testfiles", "5tglA.pdb"))
    ent_2_full = io.LoadPDB(os.path.join("testfiles", "5tglA_modified.pdb"))
    exp_atom_count = 253
    # check CA-only chain
    ent_1 = ent_1_full.Select("cname=A")
    ent_2 = ent_2_full.Select("cname=A")
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_1, ent_2)
    self.assertEqual(view1.atom_count, exp_atom_count)
    self.assertEqual(view2.atom_count, exp_atom_count)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_1, ent_2)
    self.assertEqual(view1.atom_count, exp_atom_count)
    self.assertEqual(view2.atom_count, exp_atom_count)
    # check ligand chain
    ent_1_ligand = ent_1_full.Select("cname='_'")
    ent_2_ligand = ent_2_full.Select("cname='_'")
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_1_ligand, ent_2_ligand)
    self.assertEqual(view1.atom_count, 0)
    self.assertEqual(view2.atom_count, 0)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_1_ligand, ent_2_ligand)
    self.assertEqual(view1.atom_count, 0)
    self.assertEqual(view2.atom_count, 0)
    # check both together
    view1, view2 = mol.alg.MatchResidueByLocalAln(ent_1_full, ent_2_full)
    self.assertEqual(view1.atom_count, exp_atom_count)
    self.assertEqual(view2.atom_count, exp_atom_count)
    view1, view2 = mol.alg.MatchResidueByGlobalAln(ent_1_full, ent_2_full)
    self.assertEqual(view1.atom_count, exp_atom_count)
    self.assertEqual(view2.atom_count, exp_atom_count)
    # try case where local alignment fails
    ev1 = ent_1.Select('rindex<2')   # seq = GI
    ev2 = ent_2.Select('rindex=2:3') # seq = RA
    view1, view2 = mol.alg.MatchResidueByLocalAln(ev1, ev2)
    self.assertEqual(view1.atom_count, 0)
    self.assertEqual(view2.atom_count, 0)

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
