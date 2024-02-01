''' Unit tests for the Naccess wrapper
'''

import unittest
import sys
import shutil
import tempfile

from ost.bindings import naccess
from ost import io, settings

class TestNaccessBindings(unittest.TestCase):
    def testTempDirWithDot(self):
        # naccess itself has problems with '.' in the path. So we need to test
        # that it works when data is used which has a correpsonding path.
        na_tmp_dir = tempfile.mkdtemp(prefix="ih.")
        def cleanup():
            shutil.rmtree(na_tmp_dir)
        self.addCleanup(cleanup)
        ost_ent = io.LoadPDB('testfiles/testprotein.pdb')
        excp_raised = False
        try:
            sasa = naccess.CalculateSurfaceArea(ost_ent, 
                                                scratch_dir=na_tmp_dir)
        except:
            excp_raised = True
            raise
        self.assertEqual(excp_raised, False, msg="Naccess raised an "+
                          "exception on a path containing a '.'. This is not "+
                          "supposed to happen.")

if __name__ == "__main__":
    try:
        settings.Locate("naccess")
    except:
        print("Could not find NACCESS, could not test binding...")
        sys.exit(0)
    from ost import testutils
    testutils.RunTests()
