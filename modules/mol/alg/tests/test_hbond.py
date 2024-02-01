from ost import io,mol
import unittest
import os

class TestHBond(unittest.TestCase):
  
  def setUp(self):
    p=io.IOProfile(dialect='CHARMM')
    self.eh1=io.LoadPDB(os.path.join("testfiles","hbond1.pdb"),profile=p)
    self.eh2=io.LoadPDB(os.path.join("testfiles","hbond2.pdb"),profile=p)
    self.hb_da_dict=mol.alg.hbond.BuildCHARMMHBondDonorAcceptorDict()
  
  def testHBondList(self):
    hbl1=mol.alg.hbond.GetHbondListFromView(self.eh1.Select('cname=PROA'),self.hb_da_dict)
    for hb in hbl1:self.assertTrue(hb.IsFormed())
    hbl=mol.alg.hbond.GetHbondListBetweenViews(self.eh1,self.eh2,self.hb_da_dict)
    self.assertEqual(len(hbl),1)
    
  def testHBondScore(self):
    self.assertEqual(mol.alg.hbond.CalculateHBondScore(self.eh1,self.eh2,hbond_donor_acceptor_dict=self.hb_da_dict),0.8)
  

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
