'''Unit tests for the Antechamber submodule.'''

import unittest
import os, sys

from ost.mol import mm

class TestAntechamber(unittest.TestCase):
    ###########################################################################
    # HELPERS
    def _CompareInteractions(self, int_new, int_ref):
        self.assertEqual(int_new.IsParametrized(), int_ref.IsParametrized())
        if int_new.IsParametrized():
            self.assertEqual(int_new.GetNames(), int_ref.GetNames())
            self.assertEqual(int_new.GetTypes(), int_ref.GetTypes())
            self.assertEqual(int_new.GetFuncType(), int_ref.GetFuncType())
            params_new = int_new.GetParam()
            params_ref = int_ref.GetParam()
            self.assertEqual(len(params_new), len(params_ref))
            for p_new, p_ref in zip(params_new, params_ref):
                self.assertAlmostEqual(p_new, p_ref)

    def _CompareInteractionLists(self, int_list_new, int_list_ref):
        self.assertEqual(len(int_list_new), len(int_list_ref))
        for int_new, int_ref in zip(int_list_new, int_list_ref):
            self._CompareInteractions(int_new, int_ref)

    def _CompareBuildingBlocks(self, ff_new, ff_ref, res_name):
        # get BuildingBlocks
        bb_new = ff_new.GetBuildingBlock(res_name)
        bb_ref = ff_ref.GetBuildingBlock(res_name)
        # get atoms and LJs
        self.assertEqual(len(bb_new.GetAtoms()), len(bb_ref.GetAtoms()))
        for aname, aname_ref in zip(bb_new.GetAtoms(), bb_ref.GetAtoms()):
            self.assertEqual(aname, aname_ref)
            atype = bb_new.GetType(aname)
            self.assertEqual(atype, bb_ref.GetType(aname))
            self.assertEqual(bb_new.GetMass(aname), bb_ref.GetMass(aname))
            self.assertEqual(ff_new.GetMass(atype), ff_ref.GetMass(atype))
            self.assertEqual(bb_new.GetCharge(aname), bb_ref.GetCharge(aname))
            self._CompareInteractions(ff_new.GetLJ(atype), ff_ref.GetLJ(atype))
        # compare all types of interactions
        self._CompareInteractionLists(bb_new.GetBonds(), bb_ref.GetBonds())
        self._CompareInteractionLists(bb_new.GetAngles(), bb_ref.GetAngles())
        self._CompareInteractionLists(bb_new.GetDihedrals(),
                                      bb_ref.GetDihedrals())
        self._CompareInteractionLists(bb_new.GetImpropers(),
                                      bb_ref.GetImpropers())
        # those below are expected to be 0 btw
        self._CompareInteractionLists(bb_new.GetCMaps(), bb_ref.GetCMaps())
        self._CompareInteractionLists(bb_new.GetConstraints(),
                                      bb_ref.GetConstraints())
        self._CompareInteractionLists(bb_new.GetExclusions(),
                                      bb_ref.GetExclusions())
    ###########################################################################

    def testParseFromFiles(self):
        # get new ff and compare with ref.
        ff_ref = mm.Forcefield.Load('TYR/ff_AA.dat')
        # self check
        self._CompareBuildingBlocks(ff_ref, ff_ref, 'TYR')
        # check parsed files
        ff = mm.Forcefield()
        mm.antechamber.AddFromFiles(ff, 'TYR/frcmod', 'TYR/out.mpdb')
        self._CompareBuildingBlocks(ff, ff_ref, 'TYR')

    def testParseFromPath(self):
        # check parsed path
        ff = mm.Forcefield()
        mm.antechamber.AddFromPath(ff, 'TYR')
        ff_ref = mm.Forcefield.Load('TYR/ff_AA.dat')
        self._CompareBuildingBlocks(ff, ff_ref, 'TYR')

    def testAddIon(self):
        # add pseudo iron and see
        ff = mm.Forcefield()
        mm.antechamber.AddIon(ff, 'FE2', 'FE', 55, 2, 1, 0)
        bb = ff.GetBuildingBlock('FE2')
        # check block
        self.assertEqual(bb.GetAtoms(), ['FE'])
        self.assertEqual(bb.GetTypes(), ['FE'])
        self.assertEqual(bb.GetMasses(), [55])
        self.assertEqual(bb.GetCharges(), [2])
        # nothing else...
        self.assertEqual(len(bb.GetAngles()), 0)
        self.assertEqual(len(bb.GetCMaps()), 0)
        self.assertEqual(len(bb.GetConstraints()), 0)
        self.assertEqual(len(bb.GetImpropers()), 0)
        self.assertEqual(len(bb.GetDihedrals()), 0)
        self.assertEqual(len(bb.GetBonds()), 0)
        self.assertEqual(len(bb.GetExclusions()), 0)
        # check force field
        lj = ff.GetLJ('FE')
        self.assertTrue(lj.IsParametrized())
        self.assertEqual(lj.GetTypes(), ['FE'])
        self.assertEqual(lj.GetFuncType(), mm.FuncType.LJ)
        self.assertEqual(lj.GetParam(), [1, 0])

if __name__ == "__main__":
    from ost import testutils
    testutils.RunTests()
