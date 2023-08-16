import unittest

from ost import mol, conop


class TestCompound(unittest.TestCase):
    def setUp(self):
        self.compound_lib=conop.GetDefaultLib()

    def testFindCompound(self):
        compound=self.compound_lib.FindCompound('***')
        self.assertEqual(compound, None)
        compound=self.compound_lib.FindCompound('ALA')
        self.assertNotEqual(compound, None)
        self.assertEqual(compound.id, 'ALA')
        self.assertEqual(compound.three_letter_code, 'ALA')
        self.assertEqual(compound.one_letter_code, 'A')
        self.assertTrue(compound.IsPeptideLinking())
        self.assertEqual(compound.dialect, 'PDB')
        self.assertEqual(compound.formula, 'C3 H7 N O2')
        self.assertEqual(compound.chem_class, mol.L_PEPTIDE_LINKING)
        self.assertEqual(compound.inchi,
                        "InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1")
        self.assertEqual(compound.inchi_key, "QNAYBMKLOCPYGJ-REOHCLBHSA-N")
        self.assertEqual(compound.smiles, "C[C@@H](C(=O)O)N"  )

    def testFindCompoundBySMILES(self):
        """ Test FindCompound by="smiles"."""
        compound = self.compound_lib.FindCompound('O', by="smiles")
        self.assertNotEqual(compound, None)
        self.assertEqual(compound.smiles, 'O')

        # Now we should prefer a non-obsolete compound
        # Default ordering has DIS as first pick but FindCompound should sort
        # active compounds first.
        # This assumes there are non-obsolete O/HOH compounds in the compound
        # lib, which should always be the case.
        self.assertFalse(compound.obsolete)

    def testFindCompoundByInChI(self):
        """ Test FindCompound by="inchi_code|key"."""
        inchi_code = "InChI=1/H2O/h1H2"
        inchi_key = "XLYOFNOQVPJJNP-UHFFFAOYAF"
        compound = self.compound_lib.FindCompound(inchi_code, by="inchi_code")
        self.assertNotEqual(compound, None)
        self.assertEqual(compound.inchi, inchi_code)
        self.assertEqual(compound.inchi_key, inchi_key)

     
if __name__=='__main__':
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound library available. Ignoring test_compound.py tests.')