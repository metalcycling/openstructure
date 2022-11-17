import unittest, os, sys
import ost
from ost import io, mol, settings, conop, seq
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg import stereochemistry
except ImportError:
    print("Failed to import stereochemistry. Happens when numpy missing. " \
          "Ignoring test_stereochemistry.py tests.")
    sys.exit(0)


def _LoadFile(file_name):
    """Helper to avoid repeating input path over and over."""
    return io.LoadPDB(os.path.join('testfiles', file_name))

class TestStereochemistry(unittest.TestCase):

    def test_GetResidueType(self):
        """ internal function in stereochemistry module, nevertheless crucial...
        """
        ent = _LoadFile("10b2.pdb")

        # get some residues
        ile_a_199 = ent.FindResidue("A", mol.ResNum(199))
        self.assertEqual(ile_a_199.GetName(), "ILE")
        pro_a_200 = ent.FindResidue("A", mol.ResNum(200))
        self.assertEqual(pro_a_200.GetName(), "PRO")
        glu_a_201 = ent.FindResidue("A", mol.ResNum(201))
        self.assertEqual(glu_a_201.GetName(), "GLU")
        pro_a_202 = ent.FindResidue("A", mol.ResNum(202))
        self.assertEqual(pro_a_202.GetName(), "PRO")
        ser_a_221 = ent.FindResidue("A", mol.ResNum(221))
        self.assertEqual(ser_a_221.GetName(), "SER")
        gly_a_222 = ent.FindResidue("A", mol.ResNum(222))
        self.assertEqual(gly_a_222.GetName(), "GLY")
        arg_a_223 = ent.FindResidue("A", mol.ResNum(223))
        self.assertEqual(arg_a_223.GetName(), "ARG")
        c_b_75 = ent.FindResidue("B", mol.ResNum(75))
        self.assertEqual(c_b_75.GetName(), "C")
        a_b_76 = ent.FindResidue("B", mol.ResNum(76))
        self.assertEqual(a_b_76.GetName(), "A")

        # check PRO_CIS
        res_type = stereochemistry._GetResidueType([ile_a_199.FindAtom("C"),
                                                    pro_a_200.FindAtom("N")])
        self.assertEqual(res_type, "PRO_CIS")
        res_type = stereochemistry._GetResidueType([pro_a_200.FindAtom("N"),
                                                    ile_a_199.FindAtom("C")])
        self.assertEqual(res_type, "PRO_CIS")
        res_type = stereochemistry._GetResidueType([ile_a_199.FindAtom("O"),
                                                    ile_a_199.FindAtom("C"),
                                                    pro_a_200.FindAtom("N")])
        self.assertEqual(res_type, "PRO_CIS")
        res_type = stereochemistry._GetResidueType([pro_a_200.FindAtom("N"),
                                                    ile_a_199.FindAtom("C"),
                                                    ile_a_199.FindAtom("O")])
        self.assertEqual(res_type, "PRO_CIS")

        # check PRO_TRANS
        res_type = stereochemistry._GetResidueType([glu_a_201.FindAtom("C"),
                                                    pro_a_202.FindAtom("N")])
        self.assertEqual(res_type, "PRO_TRANS")
        res_type = stereochemistry._GetResidueType([pro_a_202.FindAtom("N"),
                                                    glu_a_201.FindAtom("C")])
        self.assertEqual(res_type, "PRO_TRANS")
        res_type = stereochemistry._GetResidueType([glu_a_201.FindAtom("O"),
                                                    glu_a_201.FindAtom("C"),
                                                    pro_a_202.FindAtom("N")])
        self.assertEqual(res_type, "PRO_TRANS")
        res_type = stereochemistry._GetResidueType([pro_a_202.FindAtom("N"),
                                                    glu_a_201.FindAtom("C"),
                                                    glu_a_201.FindAtom("O")])
        self.assertEqual(res_type, "PRO_TRANS")

        # check GLY
        res_type = stereochemistry._GetResidueType([ser_a_221.FindAtom("C"),
                                                    gly_a_222.FindAtom("N")])
        self.assertEqual(res_type, "GLY")
        res_type = stereochemistry._GetResidueType([gly_a_222.FindAtom("N"),
                                                    ser_a_221.FindAtom("C")])
        self.assertEqual(res_type, "GLY")
        res_type = stereochemistry._GetResidueType([ser_a_221.FindAtom("O"),
                                                    ser_a_221.FindAtom("C"),
                                                    gly_a_222.FindAtom("N")])
        self.assertEqual(res_type, "GLY")
        res_type = stereochemistry._GetResidueType([gly_a_222.FindAtom("N"),
                                                    ser_a_221.FindAtom("C"),
                                                    ser_a_221.FindAtom("O")])
        self.assertEqual(res_type, "GLY")

        # check PEPTIDE
        res_type = stereochemistry._GetResidueType([arg_a_223.FindAtom("C"),
                                                    gly_a_222.FindAtom("N")])
        self.assertEqual(res_type, "PEPTIDE")
        res_type = stereochemistry._GetResidueType([gly_a_222.FindAtom("N"),
                                                    arg_a_223.FindAtom("C")])
        self.assertEqual(res_type, "PEPTIDE")
        res_type = stereochemistry._GetResidueType([arg_a_223.FindAtom("O"),
                                                    arg_a_223.FindAtom("C"),
                                                    gly_a_222.FindAtom("N")])
        self.assertEqual(res_type, "PEPTIDE")
        res_type = stereochemistry._GetResidueType([gly_a_222.FindAtom("N"),
                                                    arg_a_223.FindAtom("C"),
                                                    arg_a_223.FindAtom("O")])
        self.assertEqual(res_type, "PEPTIDE")

        # check NA
        res_type = stereochemistry._GetResidueType([a_b_76.FindAtom("P"),
                                                    c_b_75.FindAtom("O3'")])
        self.assertEqual(res_type, "NA")
        res_type = stereochemistry._GetResidueType([c_b_75.FindAtom("O3'"),
                                                    a_b_76.FindAtom("P")])
        self.assertEqual(res_type, "NA")
        res_type = stereochemistry._GetResidueType([a_b_76.FindAtom("O5'"),
                                                    a_b_76.FindAtom("P"),
                                                    c_b_75.FindAtom("O3'")])
        self.assertEqual(res_type, "NA")
        res_type = stereochemistry._GetResidueType([c_b_75.FindAtom("O3'"),
                                                    a_b_76.FindAtom("P"),
                                                    a_b_76.FindAtom("O5'")])
        self.assertEqual(res_type, "NA")

    def test_FullParameters(self):
        """ Tests whether we get parameters for each bond/angle
        """
        ent = _LoadFile("10b2.pdb")
        data = stereochemistry.GetDefaultStereoData()
        link_data = stereochemistry.GetDefaultStereoLinkData()
        # do bonds
        for b in ent.bonds:
            param = stereochemistry.GetBondParam(b.first, b.second,
                                                 stereo_data=data,
                                                 stereo_link_data = link_data)
            self.assertTrue(None not in param)
        # do angles
        angles = stereochemistry._GetAngles(ent.bonds)
        for a in angles:
            param = stereochemistry.GetAngleParam(a[0], a[1], a[2],
                                                  stereo_data=data,
                                                  stereo_link_data = link_data)
            self.assertTrue(None not in param)

    def test_StereoCheck(self):
        """ Test StereoCheck

        GetClashes, GetBadBonds and GetBadAngles get tested implicitely
        """
        ent = _LoadFile("10b2.pdb")
        sterochecked_ent, clashes, bad_bonds, bad_angles = \
        stereochemistry.StereoCheck(ent)

        self.assertEqual(len(clashes), 0)
        self.assertEqual(len(bad_bonds), 0)
        self.assertEqual(len(bad_angles), 0)

        phe_A_218 = ent.FindResidue("A", mol.ResNum(218))
        c_B_75 = ent.FindResidue("B", mol.ResNum(75))

        phe_at = phe_A_218.FindAtom("CD2")
        c_at = c_B_75.FindAtom("C4'")
        ed = ent.EditXCS()
        ed.SetAtomPos(phe_at, 0.5*(phe_at.GetPos() + c_at.GetPos()))

        stereochecked_ent, clashes, bad_bonds, bad_angles = \
        stereochemistry.StereoCheck(ent)

        self.assertEqual(len(clashes), 2)
        self.assertEqual(len(bad_bonds), 2)
        self.assertEqual(len(bad_angles), 2)

        # c_B_75 gets completely removed as two backbone atoms are involved
        # in clashes. phe_A_218 has the sidechain removed
        self.assertEqual(len(ent.residues)-1,
                         len(stereochecked_ent.residues))
        c_B_75 = stereochecked_ent.FindResidue("B", mol.ResNum(75))
        self.assertFalse(c_B_75.IsValid())
        phe_A_218 = stereochecked_ent.FindResidue("A", mol.ResNum(218))
        self.assertTrue(phe_A_218.IsValid())
        phe_atoms = [a.GetName() for a in phe_A_218.atoms]
        self.assertEqual(sorted(["N", "CA", "C", "O"]),
                         sorted(phe_atoms))

if __name__ == "__main__":
    from ost import testutils
    if testutils.SetDefaultCompoundLib():
        testutils.RunTests()
    else:
        print('No compound library available. Ignoring test_stereochemistry.py tests.')
