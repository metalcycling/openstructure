import unittest, os, sys
import ost
from ost import conop
import subprocess
import tempfile
import warnings


class TestCompLib(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        prefix_path = ost.GetPrefixPath()
        chemdict_tool_path = os.path.join(prefix_path, "bin", "chemdict_tool")
        if not os.path.exists(chemdict_tool_path):
            raise RuntimeError("Expect chemdict_tool:", chemdict_tool_path)
        cls.tmp_dir = tempfile.TemporaryDirectory()
        compounds_path = os.path.join("testfiles", "test_compounds.cif")
        complib_path = os.path.join(cls.tmp_dir.name, "test_complib.dat")
        cmd = [chemdict_tool_path, "create", compounds_path, complib_path]
        subprocess.run(cmd)
        cls.complib = conop.CompoundLib.Load(complib_path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp_dir.cleanup()

    def test_three_vs_five_letter_code(self):
        complib = self.complib

        comp_001 = complib.FindCompound("001")
        comp_hello = complib.FindCompound("hello")
        comp_yolo = complib.FindCompound("yolo")

        self.assertFalse(comp_001 is None)
        self.assertFalse(comp_hello is None)
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


if __name__ == "__main__":
    from ost import testutils
    testutils.RunTests()
