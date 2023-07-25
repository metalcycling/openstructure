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

     
if __name__=='__main__':
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound library available. Ignoring test_compound.py tests.')