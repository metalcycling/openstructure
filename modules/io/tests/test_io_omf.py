import unittest
import math

from ost import geom
from ost import io

def compare_atoms(a1, a2, occupancy_thresh = 0.01, bfactor_thresh = 0.01,
                  dist_thresh = 0.001):
    if abs(a1.occupancy - a2.occupancy) > occupancy_thresh:
        return False
    if abs(a1.b_factor - a2.b_factor) > bfactor_thresh:
        return False
    if geom.Distance(a1.GetPos(), a2.GetPos()) > dist_thresh:
        return False
    if a1.is_hetatom != a2.is_hetatom:
        return False
    if a1.element != a2.element:
        return False
    return True

def compare_residues(r1, r2, at_occupancy_thresh = 0.01,
                     at_bfactor_thresh = 0.01, at_dist_thresh = 0.001,
                     skip_ss = False, skip_rnums=False):
    if r1.GetName() != r2.GetName():
        return False
    if skip_rnums is False:
        if r1.GetNumber() != r2.GetNumber():
            return False
    if skip_ss is False:
        if str(r1.GetSecStructure()) != str(r2.GetSecStructure()):
            return False
    if r1.one_letter_code != r2.one_letter_code:
        return False
    if r1.chem_type != r2.chem_type:
        return False
    if r1.chem_class != r2.chem_class:
        return False
    anames1 = [a.GetName() for a in r1.atoms]
    anames2 = [a.GetName() for a in r2.atoms]
    if sorted(anames1) != sorted(anames2):
        return False
    anames = anames1
    for aname in anames:
        a1 = r1.FindAtom(aname)
        a2 = r2.FindAtom(aname)
        if not compare_atoms(a1, a2,
                             occupancy_thresh = at_occupancy_thresh,
                             bfactor_thresh = at_bfactor_thresh,
                             dist_thresh = at_dist_thresh):
            return False
    return True

def compare_chains(ch1, ch2, at_occupancy_thresh = 0.01,
                   at_bfactor_thresh = 0.01, at_dist_thresh = 0.001,
                   skip_ss=False, skip_rnums=False):
    if len(ch1.residues) != len(ch2.residues):
        return False
    for r1, r2 in zip(ch1.residues, ch2.residues):
        if not compare_residues(r1, r2,
                                at_occupancy_thresh = at_occupancy_thresh,
                                at_bfactor_thresh = at_bfactor_thresh,
                                at_dist_thresh = at_dist_thresh,
                                skip_ss = skip_ss, skip_rnums=skip_rnums):
            return False
    return True

def compare_bonds(ent1, ent2):
    bonds1 = list()
    for b in ent1.bonds:
        bond_partners = [str(b.first), str(b.second)]
        bonds1.append([min(bond_partners), max(bond_partners), b.bond_order])
    bonds2 = list()
    for b in ent2.bonds:
        bond_partners = [str(b.first), str(b.second)]
        bonds2.append([min(bond_partners), max(bond_partners), b.bond_order])
    return sorted(bonds1) == sorted(bonds2)

def compare_ent(ent1, ent2, at_occupancy_thresh = 0.01,
                at_bfactor_thresh = 0.01, at_dist_thresh = 0.001,
                skip_ss=False, skip_cnames = False, skip_bonds = False,
                skip_rnums=False):
    chain_names_one = [ch.GetName() for ch in ent1.chains]
    chain_names_two = [ch.GetName() for ch in ent2.chains]
    if skip_cnames:
        # only check whether we have the same number of chains
        if len(chain_names_one) != len(chain_names_two):
            return False
    else:
        if chain_names_one != chain_names_two:
            return False
    for ch1, ch2 in zip(ent1.chains, ent2.chains):
        if not compare_chains(ch1, ch2,
                              at_occupancy_thresh = at_occupancy_thresh,
                              at_bfactor_thresh = at_bfactor_thresh,
                              at_dist_thresh = at_dist_thresh,
                              skip_ss=skip_ss, skip_rnums=skip_rnums):
            return False
    if not skip_bonds:
        if not compare_bonds(ent1, ent2):
            return False
    return True

class TestOMF(unittest.TestCase):

    def setUp(self):
        ent, seqres, info = io.LoadMMCIF("testfiles/mmcif/3T6C.cif.gz", 
                                         seqres=True,
                                         info=True)
        self.ent = ent
        self.seqres = seqres
        self.info = info

    def test_AU(self):
        omf = io.OMF.FromMMCIF(self.ent, self.info)
        omf_bytes = omf.ToBytes()
        loaded_omf = io.OMF.FromBytes(omf_bytes)
        loaded_ent = loaded_omf.GetAU()
        self.assertTrue(compare_ent(self.ent, loaded_ent))

    def test_default_peplib(self):
        omf = io.OMF.FromMMCIF(self.ent, self.info)
        omf_bytes = omf.ToBytes()
        omf_def_pep = io.OMF.FromMMCIF(self.ent, self.info,
                                       io.OMFOption.DEFAULT_PEPLIB)
        omf_def_pep_bytes = omf_def_pep.ToBytes()
        loaded_omf_def_pep = io.OMF.FromBytes(omf_def_pep_bytes)
        loaded_ent = loaded_omf_def_pep.GetAU()

        self.assertTrue(len(omf_def_pep_bytes) < len(omf_bytes))
        self.assertTrue(compare_ent(self.ent, loaded_ent))

    def test_lossy(self):
        omf = io.OMF.FromMMCIF(self.ent, self.info)
        omf_bytes = omf.ToBytes()
        omf_lossy = io.OMF.FromMMCIF(self.ent, self.info,
                                     io.OMFOption.LOSSY)
        omf_lossy_bytes = omf_lossy.ToBytes()
        loaded_omf_lossy = io.OMF.FromBytes(omf_lossy_bytes)
        loaded_ent = loaded_omf_lossy.GetAU()

        self.assertTrue(len(omf_lossy_bytes) < len(omf_bytes))
        self.assertFalse(compare_ent(self.ent, loaded_ent))
        max_dist = math.sqrt(3*0.05*0.05)
        self.assertTrue(compare_ent(self.ent, loaded_ent,
                                    at_dist_thresh=max_dist))

    def test_avg_bfactors(self):
        omf = io.OMF.FromMMCIF(self.ent, self.info)
        omf_bytes = omf.ToBytes()
        omf_avg_bfac = io.OMF.FromMMCIF(self.ent, self.info,
                                        io.OMFOption.AVG_BFACTORS)
        omf_avg_bfac_bytes = omf_avg_bfac.ToBytes()
        loaded_omf_avg_bfac = io.OMF.FromBytes(omf_avg_bfac_bytes)
        loaded_ent = loaded_omf_avg_bfac.GetAU()

        self.assertTrue(len(omf_avg_bfac_bytes) < len(omf_bytes))
        self.assertFalse(compare_ent(self.ent, loaded_ent))
        # just give a huge slack for bfactors and check averaging manually
        self.assertTrue(compare_ent(self.ent, loaded_ent,
                                    at_bfactor_thresh=1000))

        self.assertEqual(len(self.ent.residues), len(loaded_ent.residues))
        for r_ref, r in zip(self.ent.residues, loaded_ent.residues):
            exp_bfac = sum([a.b_factor for a in r_ref.atoms])
            exp_bfac /= r_ref.atom_count
            for a in r.atoms:
                self.assertTrue(abs(a.b_factor - exp_bfac) < 0.008)

    def test_round_bfactors(self):
        omf = io.OMF.FromMMCIF(self.ent, self.info)
        omf_bytes = omf.ToBytes()
        omf_round_bfac = io.OMF.FromMMCIF(self.ent, self.info,
                                        io.OMFOption.ROUND_BFACTORS)
        omf_round_bfac_bytes = omf_round_bfac.ToBytes()
        loaded_omf_round_bfac = io.OMF.FromBytes(omf_round_bfac_bytes)
        loaded_ent = loaded_omf_round_bfac.GetAU()

        self.assertTrue(len(omf_round_bfac_bytes) < len(omf_bytes))
        self.assertFalse(compare_ent(self.ent, loaded_ent))
        self.assertTrue(compare_ent(self.ent, loaded_ent,
                                    at_bfactor_thresh=0.5))

    def test_skip_ss(self):
        omf = io.OMF.FromMMCIF(self.ent, self.info)
        omf_bytes = omf.ToBytes()
        omf_skip_ss = io.OMF.FromMMCIF(self.ent, self.info,
                                          io.OMFOption.SKIP_SS)
        omf_skip_ss_bytes = omf_skip_ss.ToBytes()
        loaded_omf_skip_ss = io.OMF.FromBytes(omf_skip_ss_bytes)
        loaded_ent = loaded_omf_skip_ss.GetAU()

        self.assertTrue(len(omf_skip_ss_bytes) < len(omf_bytes))
        self.assertFalse(compare_ent(self.ent, loaded_ent))
        self.assertTrue(compare_ent(self.ent, loaded_ent, skip_ss=True))

    def test_infer_pep_bonds(self):
        omf = io.OMF.FromMMCIF(self.ent, self.info)
        omf_bytes = omf.ToBytes()
        omf_infer_pep_bonds = io.OMF.FromMMCIF(self.ent, self.info,
                                               io.OMFOption.INFER_PEP_BONDS)
        omf_infer_pep_bonds_bytes = omf_infer_pep_bonds.ToBytes()
        loaded_omf_infer_pep_bonds = io.OMF.FromBytes(omf_infer_pep_bonds_bytes)
        loaded_ent = loaded_omf_infer_pep_bonds.GetAU()

        self.assertTrue(len(omf_infer_pep_bonds_bytes) < len(omf_bytes))
        self.assertTrue(compare_ent(self.ent, loaded_ent))

    def test_multiple_BU(self):
        ent, seqres, info = io.LoadMMCIF("testfiles/mmcif/3imj.cif.gz", 
                                         seqres=True,
                                         info=True)

        omf = io.OMF.FromMMCIF(ent, info)
        omf_bytes = omf.ToBytes()
        omf_loaded = io.OMF.FromBytes(omf_bytes)

        # there are quite some discrepancies between PDBize and OMF
        # - chain names: PDBize has specific chain names for ligands and
        #                water etc. OMF just iterates A, B, C, D, ...
        # - skip_bonds: Thats qualified atom name based. PDBize used rnums
        #               and insertion codes for waters...
        # - skip_rnums: Again, insertion codes for waters...
        self.assertTrue(compare_ent(omf_loaded.GetBU(0),
                                    info.GetBioUnits()[0].PDBize(ent),
                                    skip_cnames=True, skip_bonds=True,
                                    skip_rnums=True))

        self.assertTrue(compare_ent(omf_loaded.GetBU(1),
                                    info.GetBioUnits()[1].PDBize(ent),
                                    skip_cnames=True, skip_bonds=True,
                                    skip_rnums=True))

        # no check for the full guy... problem: PDBize throws all water
        # molecules in the same chain, whereas OMF keeps them separate
        # as in the chains from the assymetric unit... maybe needs some
        # thinking on how to resolve discrepancies between PDBize and OMF
        #self.assertTrue(compare_ent(omf_loaded.GetBU(2),
        #                            info.GetBioUnits()[2].PDBize(ent),
        #                            skip_cnames=True, skip_bonds=True,
        #                            skip_rnums=True))

if __name__== '__main__':
    from ost import testutils
    if testutils.SetDefaultCompoundLib():
        testutils.RunTests()
    else:
        print('No compound library available. Ignoring test_stereochemistry.py tests.')
