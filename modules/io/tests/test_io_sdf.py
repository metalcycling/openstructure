import unittest
from ost import *
import subprocess

class TestSDF(unittest.TestCase):
  def setUp(self):
    pass

  def test_LoadEntity(self):
    ent = io.LoadSDF('testfiles/sdf/compound.sdf')
    self.assertEqual(len(ent.chains), 4)
    self.assertEqual(len(ent.atoms), 180)
    self.assertEqual(len(ent.bonds), 188)

  def test_LoadEntity_crlf(self):
    ent = io.LoadSDF('testfiles/sdf/6d5w_rank1_crlf.sdf.gz')
    self.assertEqual(len(ent.atoms), 21)
    self.assertEqual(len(ent.bonds), 24)

  def test_Charge(self):
    ent = io.LoadSDF('testfiles/sdf/simple.sdf')
    self.assertEqual(ent.FindAtom("00001_Simple Ligand", 1, "6").charge,  0)

    # Write and read charges properly
    for chg in range(-3, 4):
      ent.FindAtom("00001_Simple Ligand", 1, "6").charge = chg
      sdf_str = io.EntityToSDFStr(ent)
      ent = io.SDFStrToEntity(sdf_str)
      self.assertEqual(ent.FindAtom("00001_Simple Ligand", 1, "6").charge,  chg)

    # Only -3 to +3 is supported
    # If M CHG is implemented the following tests can be removed
    with self.assertRaises(Exception):
      ent.FindAtom("00001_Simple Ligand", 1, "6").charge = 4
      io.EntityToSDFStr(ent)

    with self.assertRaises(Exception):
      ent.FindAtom("00001_Simple Ligand", 1, "6").charge = -4
      io.EntityToSDFStr(ent)

  def test_MChg(self):
    ent = io.LoadSDF('testfiles/sdf/m_chg.sdf')
    cl_at = ent.FindAtom("00001_Simple Ligand", 1, "6")
    self.assertEqual(cl_at.charge, -1)
    # Charge from atom line is ignored
    n_at = ent.FindAtom("00001_Simple Ligand", 1, "1")
    self.assertEqual(n_at.charge, 0)
    
if __name__== '__main__':
  from ost import testutils
  testutils.RunTests()


 
