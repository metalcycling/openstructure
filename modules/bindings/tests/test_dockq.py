''' Unit tests for the DockQ wrapper
'''
import sys
import unittest

import ost
from ost import settings
from ost import io
from ost.bindings import dockq


class TestDockQBinding(unittest.TestCase):

    def testDimerExample(self):
        # runs Dimer example from Bjoerns github repo and checks result
        mdl = io.LoadPDB("testfiles/dockq_model.pdb")
        ref = io.LoadPDB("testfiles/dockq_native.pdb")
        result = dockq.DockQ(settings.Locate("DockQ.py"), mdl, ref,
                             "A", "B", "A", "B")
        self.assertEqual(result.Fnat, 0.533)
        self.assertEqual(result.native_contacts, 60)
        self.assertEqual(result.Fnonnat, 0.238)
        self.assertEqual(result.iRMS, 1.232)
        self.assertEqual(result.LRMS, 1.516)
        self.assertEqual(result.DockQ, 0.700)

    def testMultichainExample(self):
        # multichain means one or both interface partner (ligand/receptor)
        # consist of multiple chains. DockQ provides such functionality and the
        # interface of the binding can deal with it. However, I just had no time
        # for proper testing. The binding therefore raises a NotImplementedError
        # for such cases.
        mdl = io.LoadPDB("testfiles/dockq_model.pdb")
        ref = io.LoadPDB("testfiles/dockq_native.pdb")
        with self.assertRaises(NotImplementedError):
            result = dockq.DockQ(settings.Locate("DockQ.py"), mdl, ref,
                                 ["A", "B"], "B", ["A", "B"], "B")


if __name__ == "__main__":

    try:
        settings.Locate("DockQ.py")
    except:
        print("Could not find DockQ.py, could not test binding...")
        sys.exit(0)

    from ost import testutils
    testutils.RunTests()
