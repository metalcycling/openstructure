import unittest
from ost import *
import subprocess

class TestPDB(unittest.TestCase):
  def setUp(self):
    pass

  def test_compnd_parser(self):
    e=io.LoadPDB('testfiles/pdb/compnd.pdb', restrict_chains="A")
    self.assertEqual(e.GetChainCount(), 1)
    ch = e.FindChain("A");
    self.assertEqual(ch.GetIntProp("mol_id"), 1)

  def test_properly_assigns_profile_properties(self):
    io.profiles['TEST'] = io.IOProfile()
    io.profiles['TEST'].quack_mode = False
    self.assertFalse(io.profiles['TEST'].quack_mode)
    self.assertFalse(io.profiles['TEST'].Copy().quack_mode)
    io.profiles['TEST'].quack_mode = True
    self.assertTrue(io.profiles['TEST'].quack_mode)
    self.assertTrue(io.profiles['TEST'].Copy().quack_mode)
  def test_no_bond_feasibility(self):
    io.profiles['FEAS_CHECK']=io.IOProfile(processor=conop.HeuristicProcessor(check_bond_feasibility=True))
    io.profiles['NO_FEAS_CHECK']=io.IOProfile(processor=conop.HeuristicProcessor(check_bond_feasibility=False))
        
    e1=io.LoadPDB('testfiles/pdb/simple_defective.pdb', restrict_chains="A",profile='FEAS_CHECK')
    e2=io.LoadPDB('testfiles/pdb/simple_defective.pdb', restrict_chains="A",profile='NO_FEAS_CHECK')
    ed=io.LoadPDB('testfiles/pdb/simple_defective.pdb', restrict_chains="A")
    
    res1=e1.FindResidue('A',3)
    res2=e2.FindResidue('A',3)
    resd=ed.FindResidue('A',3)

    self.assertFalse(mol.BondExists(res1.FindAtom("CA"),res1.FindAtom("CB")))
    self.assertTrue(mol.BondExists(res2.FindAtom("CA"),res2.FindAtom("CB")))
    self.assertTrue(mol.BondExists(resd.FindAtom("CA"),resd.FindAtom("CB")))

  def test_remote_loading(self):

    if subprocess.call(['ping google.com -c 1'], shell=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0:
      print("No internet connection (or wrong OS) to test remote loading in "
            "io.LoadPDB: ignoring unit test")
      return

    with self.assertRaises(IOError):
      io.LoadPDB('1ake', remote=True, remote_repo="cheeseisgoodforyou")

    # let's hope that crambin stays the same forever
    crambin_pdb = io.LoadPDB('1crn', remote=True, remote_repo='pdb')
    self.assertEqual(len(crambin_pdb.residues), 46)
    self.assertEqual(len(crambin_pdb.atoms), 327)

  def test_conect(self):
    """ See whether read_conect has an effect on reading CONECT
    """
    prot = io.LoadPDB("testfiles/pdb/conect.pdb")
    res = prot.FindResidue("A", mol.ResNum(3))
    a1 = res.FindAtom("N")
    a2 = res.FindAtom("CA")
    tmp = sorted([str(a1), str(a2)])
    bond = None
    for b in prot.bonds:
      if sorted([str(b.GetFirst()), str(b.GetSecond())]) == tmp:
        bond = b
        break
    self.assertTrue(bond is not None)
    self.assertEqual(bond.bond_order, 1)

    prot = io.LoadPDB("testfiles/pdb/conect.pdb", read_conect=True)
    res = prot.FindResidue("A", mol.ResNum(3))
    a1 = res.FindAtom("N")
    a2 = res.FindAtom("CA")
    tmp = sorted([str(a1), str(a2)])
    bond = None
    for b in prot.bonds:
      if sorted([str(b.GetFirst()), str(b.GetSecond())]) == tmp:
        bond = b
        break
    self.assertTrue(bond is not None)
    self.assertEqual(bond.bond_order, 2) # now it should be two

if __name__== '__main__':
  from ost import testutils
  testutils.RunTests()


 
