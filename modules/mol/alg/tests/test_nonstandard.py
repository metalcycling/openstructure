import unittest
from ost import io, mol


class TestNonStandard(unittest.TestCase):
 
  def test_fastModified(self):
    # phoshoserine: test if we correctly strip off modifications
    tpl=io.LoadPDB('testfiles/sep.pdb')
    new_hdl=mol.CreateEntity();
    ed=new_hdl.EditXCS()
    c=ed.InsertChain('A')
    ed.AppendResidue(c, 'SER')

    err, has_cbeta=mol.alg.CopyConserved(tpl.residues[0], new_hdl.residues[0], ed)
    self.assertTrue(err)
    self.assertTrue(has_cbeta)
    residues=new_hdl.residues
    self.assertEqual(len(residues), 1)
    self.assertEqual(len(residues[0].atoms), 6)
    self.assertTrue(new_hdl.FindAtom("A", mol.ResNum(1), "N").IsValid())
    self.assertTrue(new_hdl.FindAtom("A", mol.ResNum(1), "CA").IsValid())
    self.assertTrue(new_hdl.FindAtom("A", mol.ResNum(1), "C").IsValid())
    self.assertTrue(new_hdl.FindAtom("A", mol.ResNum(1), "O").IsValid())
    self.assertTrue(new_hdl.FindAtom("A", mol.ResNum(1), "CB").IsValid())
    self.assertTrue(new_hdl.FindAtom("A", mol.ResNum(1), "OG").IsValid())
    
    
  def test_CBeta(self):
    # test if the dst residues contain cbeta, unless they are glycines
    tpl=io.LoadPDB('testfiles/cbeta.pdb')
    new_hdl=mol.CreateEntity();
    ed=new_hdl.EditXCS()
    c=ed.InsertChain('A')
    ed.AppendResidue(c, 'MET')
    ed.AppendResidue(c, 'GLY')
    ed.AppendResidue(c, 'GLY')
    ed.AppendResidue(c, 'HIS')
    err, has_cbeta=mol.alg.CopyConserved(tpl.residues[0], new_hdl.residues[0], ed)
    self.assertTrue(has_cbeta)
    self.assertTrue(err)
    err, has_cbeta=mol.alg.CopyConserved(tpl.residues[1], new_hdl.residues[1], ed)
    self.assertFalse(has_cbeta)
    self.assertTrue(err)
    err, has_cbeta=mol.alg.CopyConserved(tpl.residues[2], new_hdl.residues[2], ed)
    self.assertFalse(has_cbeta)
    self.assertTrue(err)
    err, has_cbeta=mol.alg.CopyConserved(tpl.residues[3], new_hdl.residues[3], ed)
    self.assertTrue(has_cbeta)
    self.assertTrue(err)
      
    residues=new_hdl.residues
    self.assertEqual(len(residues), 4)
    self.assertTrue(residues[0].FindAtom("CB").IsValid())
    self.assertFalse(residues[1].FindAtom("CB").IsValid())
    self.assertFalse(residues[2].FindAtom("CB").IsValid())
    self.assertTrue(residues[3].FindAtom("CB").IsValid())


  def test_CBetaFail(self):
    # make sure that we can handle cases where CB reconstruction fails

    # NOTES:
    # - SNN is (since June 2011) labeled as a modified ASN but has a weird
    #   backbone structure without any overlap in atom names with ASN
    #   -> we hence expect it to a) always fall back to CopyNonConserved
    #      and b) fail to copy any atoms (and hence also fail to)
    # - SNN also removes N of following residue which is expected to have an
    #   incomplete backbone which would make it impossible to create a CB pos.
    
    # source of file: residues A.198 and A.199 from PDB ID 2YHW
    tpl = io.LoadPDB('testfiles/cbeta_fail.pdb')
    new_hdl = mol.CreateEntity();
    ed = new_hdl.EditXCS()
    c = ed.InsertChain('A')
    ed.AppendResidue(c, 'GLY')
    ed.AppendResidue(c, 'ALA')
    ed.AppendResidue(c, 'ASN')
    ed.AppendResidue(c, 'GLY')
    ed.AppendResidue(c, 'ALA')
    # SNN to GLY
    err = mol.alg.CopyResidue(tpl.residues[0], new_hdl.residues[0], ed)
    self.assertFalse(err)
    self.assertEqual(new_hdl.residues[0].atom_count, 0)
    # SNN to ALA
    err = mol.alg.CopyResidue(tpl.residues[0], new_hdl.residues[1], ed)
    self.assertFalse(err)
    self.assertEqual(new_hdl.residues[1].atom_count, 0)
    # SNN to ASN
    err = mol.alg.CopyResidue(tpl.residues[0], new_hdl.residues[2], ed)
    self.assertFalse(err)
    self.assertEqual(new_hdl.residues[2].atom_count, 0)
    # GLY to GLY
    err = mol.alg.CopyResidue(tpl.residues[1], new_hdl.residues[3], ed)
    self.assertTrue(err)
    self.assertEqual(new_hdl.residues[3].atom_count, 3)
    self.assertFalse(new_hdl.residues[3].FindAtom("CB").IsValid())
    # GLY to ALA
    err = mol.alg.CopyResidue(tpl.residues[1], new_hdl.residues[4], ed)
    self.assertFalse(err)
    self.assertEqual(new_hdl.residues[4].atom_count, 3)
    self.assertFalse(new_hdl.residues[4].FindAtom("CB").IsValid())


  def test_CopyResidue(self):
    tpl = io.LoadPDB('testfiles/cbeta.pdb')
    new_hdl = mol.CreateEntity();
    ed = new_hdl.EditXCS()
    c = ed.InsertChain('A')
    ed.AppendResidue(c, 'MET')
    ed.AppendResidue(c, 'GLY')
    ed.AppendResidue(c, 'GLY')
    ed.AppendResidue(c, 'HIS')
    ed.AppendResidue(c, 'HIS')
    ed.AppendResidue(c, 'GLY')
    ed.AppendResidue(c, 'HIS')
    ed.AppendResidue(c, 'MET')
    
    # MET to MET
    err = mol.alg.CopyResidue(tpl.residues[0], new_hdl.residues[0], ed)
    self.assertTrue(err)
    #GLY to GLY
    err = mol.alg.CopyResidue(tpl.residues[1], new_hdl.residues[1], ed)
    self.assertTrue(err)
    # GLY to GLY
    err = mol.alg.CopyResidue(tpl.residues[2], new_hdl.residues[2], ed)
    self.assertTrue(err)
    #now we copy a HIS to a HIS
    err = mol.alg.CopyResidue(tpl.residues[3], new_hdl.residues[3], ed)
    self.assertTrue(err)
    # copy a GLY to a HIS
    err, has_cbeta = mol.alg.CopyNonConserved(tpl.residues[1], new_hdl.residues[4], ed)
    self.assertFalse(has_cbeta)
    # copy a MET to a GLY
    err = mol.alg.CopyResidue(tpl.residues[0], new_hdl.residues[5], ed)
    self.assertFalse(err)
    # copy a MET to a HIS
    err = mol.alg.CopyResidue(tpl.residues[0], new_hdl.residues[6], ed)
    self.assertFalse(err)
    # copy a GLY to a MET with adding CB
    err = mol.alg.CopyResidue(tpl.residues[1], new_hdl.residues[7], ed)
    self.assertFalse(err)
      
    residues = new_hdl.residues
    self.assertEqual(len(residues), 8)
    # MET to MET
    self.assertTrue(residues[0].FindAtom("CB").IsValid())
    #GLY to GLY
    self.assertFalse(residues[1].FindAtom("CB").IsValid())
    #now we copy a GLY to a GLY
    self.assertFalse(residues[2].FindAtom("CB").IsValid())
    #now we copy a HIS to a HIS
    self.assertTrue(residues[3].FindAtom("CB").IsValid())
    #now we copy a GLY to a HIS without adding CB
    self.assertFalse(residues[4].FindAtom("CB").IsValid())
    #now we copy a MET to a GLY
    self.assertFalse(residues[5].FindAtom("CB").IsValid())
    # copy a MET to a HIS
    self.assertTrue(residues[6].FindAtom("CB").IsValid())
    # copy a GLY to a MET with adding CB
    self.assertTrue(residues[7].FindAtom("CB").IsValid())
    

if __name__ == "__main__":
  from ost import testutils
  if testutils.DefaultCompoundLibIsSet():
    testutils.RunTests()
  else:
    print('No compound library available. Ignoring test_nonstandard.py tests.')
