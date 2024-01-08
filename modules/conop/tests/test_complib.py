import unittest, os, sys
import ost
from ost import conop
import subprocess
import tempfile
import warnings


def CreateComplib(compound_dict_path, chemlib_out_path, extra_args=None):
    prefix_path = ost.GetPrefixPath()
    chemdict_tool_path = os.path.join(prefix_path, "bin", "chemdict_tool")
    if not os.path.exists(chemdict_tool_path):
        raise RuntimeError("Expect chemdict_tool:", chemdict_tool_path)
    cmd = [chemdict_tool_path, "create", compound_dict_path, chemlib_out_path]
    if extra_args:
        cmd += extra_args
    subprocess.run(cmd, stdout=subprocess.PIPE)


class TestCompLib(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = tempfile.TemporaryDirectory()
        compound_dict_path = os.path.join("testfiles", "test_compounds.cif")
        complib_path = os.path.join(cls.tmp_dir.name, "test_complib.dat")
        CreateComplib(compound_dict_path, complib_path)
        cls.complib = conop.CompoundLib.Load(complib_path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp_dir.cleanup()

    def test_three_vs_five_letter_code(self):
        complib = self.complib

        comp_001 = complib.FindCompound("001")
        comp_A1LU6 = complib.FindCompound("A1LU6")
        comp_yolo = complib.FindCompound("yolo")

        self.assertFalse(comp_001 is None)
        self.assertFalse(comp_A1LU6 is None)
        self.assertTrue(comp_yolo is None)

    def test_smiles(self):
        complib = self.complib
        comp_001 = complib.FindCompound("001")
        self.assertTrue(comp_001.smiles == "COc1cc(cc(c1OC)OC)C(C(=O)N2CCCC[C@H]2C(=O)O[C@@H](CCCc3ccccc3)CCCc4cccnc4)(F)F")

    def test_charges(self):
        complib = self.complib
        comp_nh4 = complib.FindCompound("NH4")
        self.assertTrue(comp_nh4.atom_specs[0].charge == 1)
        self.assertTrue(comp_nh4.atom_specs[1].charge == 0)

    def test_obsolete(self):
        complib = self.complib
        comp_ox = complib.FindCompound("OX")
        # First test that we do get data
        self.assertTrue(comp_ox.smiles == "O")
        # Now the obsolete part
        self.assertTrue(comp_ox.obsolete)
        self.assertTrue(comp_ox.replaced_by == "O")
        # Ensure not everything is obsolete
        comp_001 = complib.FindCompound("001")
        self.assertFalse(comp_001.obsolete)

    def test_default_lib_version(self):
        compound_lib = conop.GetDefaultLib()
        if compound_lib is None:
            warnings.warn("Compound library not available. Some functionality may not work as expected.")
        else:
            lib_version = compound_lib.GetOSTVersionUsed()
            if lib_version < ost.__version__:
                warnings.warn("Using old version of the compound library: %s" % lib_version)

    def test_ignore_reserved(self):
        compound_dict_path = os.path.join("testfiles", "test_compounds.cif")
        complib_no_reserved_path = os.path.join(self.tmp_dir.name, "test_complib_no_reserved.dat")
        CreateComplib(compound_dict_path, complib_no_reserved_path, ["-i"])
        complib_no_reserved = conop.CompoundLib.Load(complib_no_reserved_path)

        # 01-98 are reserved
        assert self.complib.FindCompound("98") is not None
        assert complib_no_reserved.FindCompound("98") is None

        # DRG, INH and LIG are reserved
        assert self.complib.FindCompound("DRG") is not None
        assert complib_no_reserved.FindCompound("DRG") is None
        assert self.complib.FindCompound("INH") is not None
        assert complib_no_reserved.FindCompound("INH") is None
        assert self.complib.FindCompound("LIG") is not None
        assert complib_no_reserved.FindCompound("LIG") is None

        # OX is obsolete but not reserved
        assert complib_no_reserved.FindCompound("OX") is not None

        # 00, 000, 001, 010, 986, 98B are not reserved
        assert complib_no_reserved.FindCompound("00") is not None
        assert complib_no_reserved.FindCompound("000") is not None
        assert complib_no_reserved.FindCompound("001") is not None
        assert complib_no_reserved.FindCompound("010") is not None
        assert complib_no_reserved.FindCompound("986") is not None
        assert complib_no_reserved.FindCompound("98B") is not None

    def test_ignore_obsolete(self):
        compound_dict_path = os.path.join("testfiles", "test_compounds.cif")
        complib_no_obsolete_path = os.path.join(self.tmp_dir.name, "test_complib_no_obsolete.dat")
        CreateComplib(compound_dict_path, complib_no_obsolete_path, ["-o"])
        complib_no_obsolete = conop.CompoundLib.Load(complib_no_obsolete_path)

        # 01-98, DRG, INH and LIG are reserved but not obsolete
        assert complib_no_obsolete.FindCompound("98") is not None
        assert complib_no_obsolete.FindCompound("DRG") is not None
        assert complib_no_obsolete.FindCompound("INH") is not None
        assert complib_no_obsolete.FindCompound("LIG") is not None

        # OX is obsolete
        assert complib_no_obsolete.FindCompound("OX") is None


if __name__ == "__main__":
    from ost import testutils
    testutils.RunTests()
