import unittest, os, sys
import ost
from ost import conop
import subprocess
import tempfile


class TestCompLib(unittest.TestCase):

    def test_three_vs_five_letter_code(self):

        prefix_path = ost.GetPrefixPath()
        chemdict_tool_path = os.path.join(prefix_path, "bin", "chemdict_tool")
        if not os.path.exists(chemdict_tool_path):
            raise RuntimeError("Expect chemdict_tool:", chemdict_tool_path)
        tmp_dir = tempfile.TemporaryDirectory()
        compounds_path = os.path.join("testfiles", "test_compounds.cif")
        complib_path = os.path.join(tmp_dir.name, "test_complib.dat")
        cmd = [chemdict_tool_path, "create", compounds_path, complib_path]
        subprocess.run(cmd)

        complib = conop.CompoundLib.Load(complib_path)

        comp_001 = complib.FindCompound("001")
        comp_hello = complib.FindCompound("hello")
        comp_yolo = complib.FindCompound("yolo")

        self.assertFalse(comp_001 is None)
        self.assertFalse(comp_hello is None)
        self.assertTrue(comp_yolo is None)

if __name__ == "__main__":
    from ost import testutils
    if testutils.SetDefaultCompoundLib():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_complib tests.')