import unittest, os, sys
import ost
from ost import io, geom, mol
from ost.mol.alg import BUInfo

class TestBioUnit(unittest.TestCase):

  def test_bu(self):
    ent, seqres, info = io.LoadMMCIF("testfiles/1out.cif.gz", 
                                     seqres=True,
                                     info=True)

    # Create BUInfo from MMCifInfoBioUnit
    biounits = info.GetBioUnits()
    self.assertEqual(len(biounits), 1)
    bu_info = BUInfo(biounits[0])

    # directly use the dump and load mechanism
    # IF ANY UNIT TEST FAILS, DISABLE THE FOLLOWING TWO LINES
    # TO EXCLUDE THE POSSIBILITY OF ISSUES IN BUINFO IO FUNCTIONALITY
    bytes_str = bu_info.ToBytes()
    bu_info = BUInfo.FromBytes(bytes_str)

    # check whether properties in BUInfo object are correctly set
    asu_chains = bu_info.GetAUChains()
    self.assertEqual(asu_chains, [["A", "B", "C", "D", "E", "F"]])
    transformations = bu_info.GetTransformations()
    self.assertEqual(len(transformations), 1)
    self.assertEqual(len(transformations[0]), 2)
    self.assertEqual(transformations[0][0], geom.Mat4()) # identity

    # reconstruct biounit
    bu = mol.alg.CreateBU(ent, bu_info)
    self.assertEqual([ch.GetName() for ch in bu.chains],
                     ["1.A", "1.B", "1.C", "1.D", "1.E", "1.F",
                      "2.A", "2.B", "2.C", "2.D", "2.E", "2.F"])

    # extract copies of original assymetric units from biounit
    for ch in bu.chains:
      ch.SetIntProp("asu_copy_idx", int(ch.GetName().split('.')[0]))
    asu_copy_one = bu.Select("gcasu_copy_idx=1")
    asu_copy_two = bu.Select("gcasu_copy_idx=2")

    # compare bonds - well, this tests whether bonds have been correctly
    # transferred. BUT: we can be pretty sure that all atoms/residues/chains
    # have been transferred too when explicitely checking the bonds.
    asu_bonds = ent.GetBondList()
    asu_bonds_desc = list()
    for b in asu_bonds:
      at_one = b.GetFirst()
      at_one_cname = at_one.GetChain().GetName()
      at_one_rname = at_one.GetResidue().GetName()
      at_one_rnum = at_one.GetResidue().GetNumber().GetNum()
      at_one_name = at_one.GetName()
      at_one_desc = f"{at_one_cname}.{at_one_rname}.{at_one_rnum}.{at_one_name}"
      at_two = b.GetSecond()
      at_two_cname = at_two.GetChain().GetName()
      at_two_rname = at_two.GetResidue().GetName()
      at_two_rnum = at_two.GetResidue().GetNumber().GetNum()
      at_two_name = at_two.GetName()
      at_two_desc = f"{at_two_cname}.{at_two_rname}.{at_two_rnum}.{at_two_name}"
      asu_bonds_desc.append((at_one_desc, at_two_desc))

    asu_copy_one_bonds = asu_copy_one.GetBondList()
    asu_copy_one_bonds_desc = list()
    for b in asu_copy_one_bonds:
      at_one = b.GetFirst()
      at_one_cname = at_one.GetChain().GetName().split('.')[1]
      at_one_rname = at_one.GetResidue().GetName()
      at_one_rnum = at_one.GetResidue().GetNumber().GetNum()
      at_one_name = at_one.GetName()
      at_one_desc = f"{at_one_cname}.{at_one_rname}.{at_one_rnum}.{at_one_name}"
      at_two = b.GetSecond()
      at_two_cname = at_two.GetChain().GetName().split('.')[1]
      at_two_rname = at_two.GetResidue().GetName()
      at_two_rnum = at_two.GetResidue().GetNumber().GetNum()
      at_two_name = at_two.GetName()
      at_two_desc = f"{at_two_cname}.{at_two_rname}.{at_two_rnum}.{at_two_name}"
      asu_copy_one_bonds_desc.append((at_one_desc, at_two_desc))
    
    asu_copy_two_bonds = asu_copy_two.GetBondList()
    asu_copy_two_bonds_desc = list()
    for b in asu_copy_two_bonds:
      at_one = b.GetFirst()
      at_one_cname = at_one.GetChain().GetName().split('.')[1]
      at_one_rname = at_one.GetResidue().GetName()
      at_one_rnum = at_one.GetResidue().GetNumber().GetNum()
      at_one_name = at_one.GetName()
      at_one_desc = f"{at_one_cname}.{at_one_rname}.{at_one_rnum}.{at_one_name}"
      at_two = b.GetSecond()
      at_two_cname = at_two.GetChain().GetName().split('.')[1]
      at_two_rname = at_two.GetResidue().GetName()
      at_two_rnum = at_two.GetResidue().GetNumber().GetNum()
      at_two_name = at_two.GetName()
      at_two_desc = f"{at_two_cname}.{at_two_rname}.{at_two_rnum}.{at_two_name}"
      asu_copy_two_bonds_desc.append((at_one_desc, at_two_desc))

    self.assertEqual(asu_bonds_desc, asu_copy_one_bonds_desc)
    self.assertEqual(asu_bonds_desc, asu_copy_two_bonds_desc)


if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
