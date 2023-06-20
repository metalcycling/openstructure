import sys
import unittest
from ost import geom, conop, io
from ost.conop import cleanup

class TestCleanUp(unittest.TestCase):

  def setUp(self):
    self.comp_lib = conop.GetDefaultLib()
    io.profiles['DEFAULT'].processor = conop.RuleBasedProcessor(self.comp_lib)
    self.ent = io.LoadPDB("sample_test_cleanup.pdb")
    self.ent_no_wat = io.LoadPDB("sample_nowater.pdb")
    self.ent_no_lig = io.LoadPDB("sample_noligands.pdb")

  def testStripWater(self):
    self.new_ent = cleanup.Cleanup(self.ent, strip_water=True, canonicalize=False, remove_ligands=False)
    self.assertEqual( self.new_ent.residue_count, self.ent_no_wat.residue_count )
    self.assertTrue( self.new_ent.residues[0].IsValid() )
    self.assertEqual( self.new_ent.residues[0].qualified_name, self.ent_no_wat.residues[0].qualified_name)
    self.assertTrue( self.new_ent.residues[1].IsValid() )
    self.assertEqual( self.new_ent.residues[1].qualified_name, self.ent_no_wat.residues[1].qualified_name)
    self.assertTrue( self.new_ent.residues[2].IsValid() )
    self.assertEqual( self.new_ent.residues[2].qualified_name, self.ent_no_wat.residues[2].qualified_name)
    self.assertTrue( self.new_ent.residues[3].IsValid() )
    self.assertEqual( self.new_ent.residues[3].qualified_name, self.ent_no_wat.residues[3].qualified_name)
    self.assertTrue( self.new_ent.residues[4].IsValid() )
    self.assertEqual( self.new_ent.residues[4].qualified_name, self.ent_no_wat.residues[4].qualified_name)
    self.assertTrue( self.new_ent.residues[5].IsValid() )
    self.assertEqual( self.new_ent.residues[5].qualified_name, self.ent_no_wat.residues[5].qualified_name)
    self.assertTrue( self.new_ent.residues[6].IsValid() )
    self.assertEqual( self.new_ent.residues[6].qualified_name, self.ent_no_wat.residues[6].qualified_name)
    self.assertTrue( self.new_ent.residues[7].IsValid() )
    self.assertEqual( self.new_ent.residues[7].qualified_name, self.ent_no_wat.residues[7].qualified_name)

  def testCanonicalize(self):
    self.new_ent = cleanup.Cleanup(self.ent, strip_water=False, canonicalize=True, remove_ligands=False)
    #standard residue must be the same
    self.gly = self.ent.residues[1]
    self.new_gly = self.new_ent.residues[1]
    self.assertTrue(self.new_gly.IsValid())
    self.assertTrue(self.new_gly.IsPeptideLinking())
    self.assertEqual(self.gly.atom_count, self.new_gly.atom_count)
    #TEMP del sidechain of incomplete residue and OXT if present
    self.new_cys = self.new_ent.residues[4]
    self.new_cys_atoms = set([atm.name for atm in self.new_cys.atoms])
    self.assertEqual( len(self.new_cys_atoms), 4, msg = repr(self.new_cys_atoms))
    self.assertTrue( "CB" in self.new_cys_atoms)
    self.assertTrue( "CA" in self.new_cys_atoms)
    self.assertTrue( "C" in self.new_cys_atoms)
    self.assertFalse( "OXT" in self.new_cys_atoms)
    self.assertTrue( "N" in self.new_cys_atoms)
    #test replacement of atoms
    self.mse = self.ent.residues[0]
#    self.assertTrue( self.mse.IsValid())
#    self.assertTrue( self.mse.IsPeptideLinking())
    self.sel = self.mse.FindAtom("SE")
#    self.assertTrue( self.sel.IsValid())
    self.met = self.new_ent.residues[0]
    self.assertTrue(self.met.IsValid())
    self.assertEqual(self.mse.atom_count, self.met.atom_count)
    self.assertEqual(self.met.name, "MET")
    self.assertEqual(self.met.one_letter_code, "M")
    self.assertTrue(self.met.IsPeptideLinking())
    self.sul = self.met.FindAtom("SD")
    self.assertTrue(self.sul.IsValid())
    self.assertTrue(geom.Equal(self.sul.pos,self.sel.pos), msg = "sul:%s sel:%s"%(str(self.sul.pos), str(self.sel.pos)) )
    self.assertEqual(self.sul.element, "S")
#    self.AssertTrue( sul.mass == conop.Conopology.Instance().GetDefaultAtomMass("S"))
#    self.AssertTrue( sul.radius == conop.Conopology.Instance().GetDefaultAtomRadius("S"))
    for atm in self.met.atoms:
      self.assertFalse( atm.is_hetatom)
    #test addition
    self.mly = self.ent.residues[2]
#    self.assertTrue( self.mly.IsValid())
#    self.assertTrue( self.mly.IsPeptideLinking())
    self.new_lys = self.new_ent.residues[2]
    self.assertTrue(self.new_lys.IsValid())
    self.assertTrue(self.new_lys.IsPeptideLinking())
    self.assertEqual(self.new_lys.name, "LYS")
    self.assertEqual(self.new_lys.one_letter_code, "K")
    self.new_lys_atoms = set([atm.name for atm in self.new_lys.atoms])
    self.canon_lys = self.comp_lib.FindCompound("LYS")
    self.canon_lys_atoms = set([atom.name for atom in self.canon_lys.atom_specs
                                             if atom.element != "H" and atom.element != "D"  and not atom.is_leaving ])
    self.assertEqual(self.canon_lys_atoms, self.new_lys_atoms)
    self.assertFalse(self.canon_lys_atoms - self.new_lys_atoms)
    self.assertFalse(self.new_lys_atoms - self.canon_lys_atoms) #test the reverse
    for atm in self.new_lys.atoms:
      self.assertFalse( atm.is_hetatom)
    #deletions
    self.dha = self.ent.residues[3]
#    self.assertTrue( self.dha.IsValid())
#    self.assertTrue( self.dha.IsPeptideLinking())
    self.new_ser = self.new_ent.residues[3]
    self.assertTrue(self.new_ser.IsValid())
    self.assertTrue(self.new_ser.IsPeptideLinking())
    self.assertEqual(self.new_ser.name, "SER")
    self.assertEqual(self.new_ser.one_letter_code, "S")
    self.new_ser_atoms = set([atm.name for atm in self.new_ser.atoms])
    self.canon_ser = self.comp_lib.FindCompound("SER")
    self.canon_ser_atoms = set([atom.name for atom in self.canon_ser.atom_specs
                                             if atom.element != "H" and atom.element != "D"  and not atom.is_leaving ])
    #TEMP
    self.assertEqual( len(self.new_ser_atoms), 5)
    self.assertTrue( "CB" in self.new_ser_atoms)
    self.assertTrue( "CA" in self.new_ser_atoms)
    self.assertTrue( "C" in self.new_ser_atoms)
    self.assertTrue( "O" in self.new_ser_atoms)
    self.assertTrue( "N" in self.new_ser_atoms)
    #AFTER TEMP
    #self.assertEqual( self.canon_ser_atoms, self.new_ser_atoms)
    #self.assertFalse(self.canon_ser_atoms - self.new_ser_atoms)
    #self.assertFalse(self.new_ser_atoms - self.canon_ser_atoms) #test the reverse
    for atm in self.new_ser.atoms:
      self.assertFalse( atm.is_hetatom)
    #test deletion of whole residue
    self.assertEqual(self.ent.residues[5].chem_class, "D_PEPTIDE_LINKING")
    self.assertNotEqual(self.new_ent.residues[5].name, "DAL")
    self.assertNotEqual(self.ent.residue_count, self.new_ent.residue_count)

  def testRemoveLigands(self):
    self.new_ent = cleanup.Cleanup(self.ent, strip_water=False, canonicalize=False, remove_ligands=True)
    self.assertEqual(self.new_ent.residue_count, self.ent_no_lig.residue_count )
    #MSE
    self.assertTrue(self.new_ent.residues[0].IsValid() )
    self.assertEqual(self.new_ent.residues[0].qualified_name, self.ent_no_lig.residues[0].qualified_name)
    self.assertTrue(self.new_ent.residues[0].IsPeptideLinking())
    self.assertTrue(self.new_ent.residues[0].atoms[0].is_hetatom)
    #GLY
    self.assertTrue(self.new_ent.residues[1].IsValid() )
    self.assertEqual(self.new_ent.residues[1].qualified_name, self.ent_no_lig.residues[1].qualified_name)
    self.assertTrue(self.new_ent.residues[1].IsPeptideLinking())
    self.assertFalse(self.new_ent.residues[1].atoms[0].is_hetatom)
    #MLY
    self.assertTrue(self.new_ent.residues[2].IsValid() )
    self.assertEqual(self.new_ent.residues[2].qualified_name, self.ent_no_lig.residues[2].qualified_name)
    self.assertTrue(self.new_ent.residues[2].IsPeptideLinking())
    self.assertTrue(self.new_ent.residues[2].atoms[0].is_hetatom)
    #DHA
    self.assertTrue(self.new_ent.residues[3].IsValid() )
    self.assertEqual(self.new_ent.residues[3].qualified_name, self.ent_no_lig.residues[3].qualified_name)
    self.assertTrue(self.new_ent.residues[3].IsPeptideLinking())
    self.assertTrue(self.new_ent.residues[3].atoms[0].is_hetatom)
    #CYS
    self.assertTrue(self.new_ent.residues[4].IsValid() )
    self.assertEqual(self.new_ent.residues[4].qualified_name, self.ent_no_lig.residues[4].qualified_name)
    self.assertTrue(self.new_ent.residues[4].IsPeptideLinking())
    self.assertFalse(self.new_ent.residues[4].atoms[0].is_hetatom)
    #DAL
    self.assertTrue(self.new_ent.residues[5].IsValid() )
    self.assertEqual(self.new_ent.residues[5].qualified_name, self.ent_no_lig.residues[5].qualified_name)
    self.assertTrue(self.new_ent.residues[5].IsPeptideLinking())
    self.assertTrue(self.new_ent.residues[5].atoms[0].is_hetatom)
    #HOH
    self.assertTrue(self.new_ent.residues[6].IsValid() )
    self.assertEqual(self.new_ent.residues[6].qualified_name, self.ent_no_lig.residues[6].qualified_name)
    self.assertFalse(self.new_ent.residues[6].IsPeptideLinking()) # here assertFalse instead of assertTrue
    self.assertTrue(self.new_ent.residues[6].atoms[0].is_hetatom)


if __name__ == "__main__":
  from ost import testutils
  if testutils.DefaultCompoundLibIsSet():
    testutils.RunTests()
  else:
    print('No compound library available. Ignoring test_cleanup.py tests.')
