import unittest, os, sys
from functools import lru_cache

import numpy as np

from ost import io, mol, geom
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg.ligand_scoring import *
    from ost.mol.alg import ligand_scoring
except ImportError:
    print("Failed to import ligand_scoring.py. Happens when numpy, scipy or "
          "networkx is missing. Ignoring test_ligand_scoring.py tests.")
    sys.exit(0)


def _GetTestfilePath(filename):
    """Get the path to the test file given filename"""
    return os.path.join('testfiles', filename)


@lru_cache(maxsize=None)
def _LoadMMCIF(filename):
    path = _GetTestfilePath(filename)
    ent, seqres = io.LoadMMCIF(path, seqres=True)
    return ent


@lru_cache(maxsize=None)
def _LoadPDB(filename):
    path = _GetTestfilePath(filename)
    ent = io.LoadPDB(path)
    return ent


@lru_cache(maxsize=None)
def _LoadEntity(filename):
    path = _GetTestfilePath(filename)
    ent = io.LoadEntity(path)
    return ent


class TestLigandScoring(unittest.TestCase):

    def test_extract_ligands_mmCIF(self):
        """Test that we can extract ligands from mmCIF files.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        sc = LigandScorer(mdl, trg, None, None)

        assert len(sc.target_ligands) == 7
        assert len(sc.model_ligands) == 1
        assert len([r for r in sc.target.residues if r.is_ligand]) == 7
        assert len([r for r in sc.model.residues if r.is_ligand]) == 1

    def test_init_given_ligands(self):
        """Test that we can instantiate the scorer with ligands contained in
        the target and model entity and given in a list.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        # Pass entity views
        trg_lig = [trg.Select("rname=MG"), trg.Select("rname=G3D")]
        mdl_lig = [mdl.Select("rname=G3D")]
        sc = LigandScorer(mdl, trg, mdl_lig, trg_lig)

        assert len(sc.target_ligands) == 4
        assert len(sc.model_ligands) == 1
        # IsLigand flag should still be set even on not selected ligands
        assert len([r for r in sc.target.residues if r.is_ligand]) == 7
        assert len([r for r in sc.model.residues if r.is_ligand]) == 1

        # Ensure the residues are not copied
        assert len(sc.target.Select("rname=MG").residues) == 2
        assert len(sc.target.Select("rname=G3D").residues) == 2
        assert len(sc.model.Select("rname=G3D").residues) == 1

        # Pass residue handles
        trg_lig = [trg.FindResidue("F", 1), trg.FindResidue("H", 1)]
        mdl_lig = [mdl.FindResidue("L_2", 1)]
        sc = LigandScorer(mdl, trg, mdl_lig, trg_lig)

        assert len(sc.target_ligands) == 2
        assert len(sc.model_ligands) == 1

        # Ensure the residues are not copied
        assert len(sc.target.Select("rname=ZN").residues) == 1
        assert len(sc.target.Select("rname=G3D").residues) == 2
        assert len(sc.model.Select("rname=G3D").residues) == 1

    def test_init_sdf_ligands(self):
        """Test that we can instantiate the scorer with ligands from separate SDF files.

        In order to setup the ligand SDF files, the following code was used:
        for prefix in [os.path.join('testfiles', x) for x in ["1r8q", "P84080_model_02"]]:
            trg = io.LoadMMCIF("%s.cif.gz" % prefix)
            trg_prot = trg.Select("protein=True")
            io.SavePDB(trg_prot, "%s_protein.pdb.gz" % prefix)
            lig_num = 0
            for chain in trg.chains:
                if chain.chain_type == mol.ChainType.CHAINTYPE_NON_POLY:
                    lig_sel = trg.Select("cname=%s" % chain.name)
                    lig_ent = mol.CreateEntityFromView(lig_sel, False)
                    io.SaveEntity(lig_ent, "%s_ligand_%d.sdf" % (prefix, lig_num))
                    lig_num += 1
        """
        mdl = _LoadPDB("P84080_model_02_nolig.pdb")
        mdl_ligs = [_LoadEntity("P84080_model_02_ligand_0.sdf")]
        trg = _LoadPDB("1r8q_protein.pdb.gz")
        trg_ligs = [_LoadEntity("1r8q_ligand_%d.sdf" % i) for i in range(7)]

        # Pass entities
        sc = LigandScorer(mdl, trg, mdl_ligs, trg_ligs)

        assert len(sc.target_ligands) == 7
        assert len(sc.model_ligands) == 1
        # Ensure we set the is_ligand flag
        assert len([r for r in sc.target.residues if r.is_ligand]) == 7
        assert len([r for r in sc.model.residues if r.is_ligand]) == 1

        # Pass residues
        mdl_ligs_res = [mdl_ligs[0].residues[0]]
        trg_ligs_res = [res for ent in trg_ligs for res in ent.residues]

        sc = LigandScorer(mdl, trg, mdl_ligs_res, trg_ligs_res)

        assert len(sc.target_ligands) == 7
        assert len(sc.model_ligands) == 1

    def test_init_reject_duplicate_ligands(self):
        """Test that we reject input if multiple ligands with the same chain
         name/residue number are given.
        """
        mdl = _LoadPDB("P84080_model_02_nolig.pdb")
        mdl_ligs = [_LoadEntity("P84080_model_02_ligand_0.sdf")]
        trg = _LoadPDB("1r8q_protein.pdb.gz")
        trg_ligs = [_LoadEntity("1r8q_ligand_%d.sdf" % i) for i in range(7)]

        # Reject identical model ligands
        with self.assertRaises(RuntimeError):
            sc = LigandScorer(mdl, trg, [mdl_ligs[0], mdl_ligs[0]], trg_ligs)

        # Reject identical target ligands
        lig0 = trg_ligs[0].Copy()
        lig1 = trg_ligs[1].Copy()
        ed1 = lig1.EditXCS()
        ed1.RenameChain(lig1.chains[0], lig0.chains[0].name)
        ed1.SetResidueNumber(lig1.residues[0], lig0.residues[0].number)
        with self.assertRaises(RuntimeError):
            sc = LigandScorer(mdl, trg, mdl_ligs, [lig0, lig1])

    def test__ResidueToGraph(self):
        """Test that _ResidueToGraph works as expected
        """
        mdl_lig = _LoadEntity("P84080_model_02_ligand_0.sdf")

        graph = ligand_scoring._ResidueToGraph(mdl_lig.residues[0])
        assert len(graph.edges) == 34
        assert len(graph.nodes) == 32
        # Check an arbitrary node
        assert [a for a in graph.adj["14"].keys()] == ["13", "29"]

        graph = ligand_scoring._ResidueToGraph(mdl_lig.residues[0], by_atom_index=True)
        assert len(graph.edges) == 34
        assert len(graph.nodes) == 32
        # Check an arbitrary node
        assert [a for a in graph.adj[13].keys()] == [12, 28]

    def test__ComputeSymmetries(self):
        """Test that _ComputeSymmetries works.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        trg_mg1 = trg.FindResidue("E", 1)
        trg_g3d1 = trg.FindResidue("F", 1)
        trg_afb1 = trg.FindResidue("G", 1)
        trg_g3d2 = trg.FindResidue("J", 1)
        mdl_g3d = mdl.FindResidue("L_2", 1)

        sym = ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1)
        assert len(sym) == 72

        sym = ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1, by_atom_index=True)
        assert len(sym) == 72

        # Test that we can match ions read from SDF
        sdf_lig = _LoadEntity("1r8q_ligand_0.sdf")
        sym = ligand_scoring._ComputeSymmetries(trg_mg1, sdf_lig.residues[0], by_atom_index=True)
        assert len(sym) == 1

        # Test that it works with views and only consider atoms in the view
        # Skip PA, PB and O[1-3]A and O[1-3]B in target and model
        # We assume atom index are fixed and won't change
        trg_g3d1_sub = trg_g3d1.Select("aindex>6019").residues[0]
        mdl_g3d_sub = mdl_g3d.Select("aindex>1447").residues[0]

        sym = ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1_sub)
        assert len(sym) == 6

        sym = ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1_sub, by_atom_index=True)
        assert len(sym) == 6

        # Substructure matches
        sym = ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1_sub, substructure_match=True)
        assert len(sym) == 6

        # Missing atoms only allowed in target, not in model
        with self.assertRaises(NoSymmetryError):
            ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1, substructure_match=True)

    def test_SCRMSD(self):
        """Test that SCRMSD works.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        trg_mg1 = trg.FindResidue("E", 1)
        trg_g3d1 = trg.FindResidue("F", 1)
        trg_afb1 = trg.FindResidue("G", 1)
        trg_g3d2 = trg.FindResidue("J", 1)
        mdl_g3d = mdl.FindResidue("L_2", 1)

        rmsd = SCRMSD(mdl_g3d, trg_g3d1)
        self.assertAlmostEqual(rmsd, 2.21341e-06, 10)
        rmsd = SCRMSD(mdl_g3d, trg_g3d2)
        self.assertAlmostEqual(rmsd, 61.21325, 4)

        # Ensure we raise a NoSymmetryError if the ligand is wrong
        with self.assertRaises(NoSymmetryError):
            SCRMSD(mdl_g3d, trg_mg1)
        with self.assertRaises(NoSymmetryError):
            SCRMSD(mdl_g3d, trg_afb1)

        # Assert that transform works
        trans = geom.Mat4(-0.999256, 0.00788487, -0.0377333, -15.4397,
                          0.0380652, 0.0473315, -0.998154, 29.9477,
                          -0.00608426, -0.998848, -0.0475963, 28.8251,
                          0, 0, 0, 1)
        rmsd = SCRMSD(mdl_g3d, trg_g3d2, transformation=trans)
        self.assertAlmostEqual(rmsd, 0.293972, 5)

        # Assert that substructure matches work
        trg_g3d1_sub = trg_g3d1.Select("aindex>6019").residues[0] # Skip PA, PB and O[1-3]A and O[1-3]B.
        # mdl_g3d_sub = mdl_g3d.Select("aindex>1447").residues[0] # Skip PA, PB and O[1-3]A and O[1-3]B.
        with self.assertRaises(NoSymmetryError):
            SCRMSD(mdl_g3d, trg_g3d1_sub)  # no full match

        # But partial match is OK
        rmsd = SCRMSD(mdl_g3d, trg_g3d1_sub, substructure_match=True)
        self.assertAlmostEqual(rmsd, 2.2376232209353475e-06, 8)

        # Ensure it doesn't work the other way around - ie incomplete model is invalid
        with self.assertRaises(NoSymmetryError):
            SCRMSD(trg_g3d1_sub, mdl_g3d)  # no full match

    def test__compute_scores(self):
        """Test that _compute_scores works.
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")
        mdl_lig = io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))
        sc = LigandScorer(mdl, trg, [mdl_lig], None)

        # Note: expect warning about Binding site of H.ZN1 not mapped to the model
        sc._compute_scores()

        # Check RMSD
        assert sc.rmsd_matrix.shape == (7, 1)
        np.testing.assert_almost_equal(sc.rmsd_matrix, np.array(
            [[np.nan],
            [0.04244993],
            [np.nan],
            [np.nan],
            [np.nan],
            [0.29399303],
            [np.nan]]), decimal=5)

        # Check lDDT-PLI
        self.assertEqual(sc.lddt_pli_matrix.shape, (7, 1))
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[0, 0]))
        self.assertAlmostEqual(sc.lddt_pli_matrix[1, 0], 0.99843, 5)
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[2, 0]))
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[3, 0]))
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[4, 0]))
        self.assertAlmostEqual(sc.lddt_pli_matrix[5, 0], 1.0)
        self.assertTrue(np.isnan(sc.lddt_pli_matrix[6, 0]))

    def test_check_resnames(self):
        """Test check_resname argument works
        """
        # 4C0A has mismatching sequence and fails with check_resnames=True
        mdl_1r8q = _LoadMMCIF("1r8q.cif.gz")
        trg_4c0a = _LoadMMCIF("4c0a.cif.gz")

        mdl = mdl_1r8q.Select("cname=D or cname=F")
        trg = trg_4c0a.Select("cname=C or cname=I")

        with self.assertRaises(RuntimeError):
            sc = LigandScorer(mdl, trg, [mdl.FindResidue("F", 1)], [trg.FindResidue("I", 1)], check_resnames=True)
            sc._compute_scores()

        sc = LigandScorer(mdl, trg, [mdl.FindResidue("F", 1)], [trg.FindResidue("I", 1)], check_resnames=False)
        sc._compute_scores()

    def test__scores(self):
        """Test that the scores are computed correctly
        """
        # 4C0A has more ligands
        trg = _LoadMMCIF("1r8q.cif.gz")
        trg_4c0a = _LoadMMCIF("4c0a.cif.gz")
        sc = LigandScorer(trg, trg_4c0a, None, None, check_resnames=False)

        expected_keys = {"J", "F"}
        self.assertFalse(expected_keys.symmetric_difference(sc.rmsd.keys()))
        self.assertFalse(expected_keys.symmetric_difference(sc.rmsd_details.keys()))
        self.assertFalse(expected_keys.symmetric_difference(sc.lddt_pli.keys()))
        self.assertFalse(expected_keys.symmetric_difference(sc.lddt_pli_details.keys()))

        # rmsd
        self.assertAlmostEqual(sc.rmsd["J"][mol.ResNum(1)], 0.8016608357429504, 5)
        self.assertAlmostEqual(sc.rmsd["F"][mol.ResNum(1)], 0.9286373257637024, 5)
        # rmsd_details
        self.assertEqual(sc.rmsd_details["J"][mol.ResNum(1)]["chain_mapping"], {'F': 'D', 'C': 'C'})
        self.assertEqual(len(sc.rmsd_details["J"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(len(sc.rmsd_details["J"][mol.ResNum(1)]["bs_ref_res_mapped"]), 15)
        self.assertEqual(len(sc.rmsd_details["J"][mol.ResNum(1)]["bs_mdl_res_mapped"]), 15)
        self.assertEqual(sc.rmsd_details["J"][mol.ResNum(1)]["target_ligand"].qualified_name, 'I.G3D1')
        self.assertEqual(sc.rmsd_details["J"][mol.ResNum(1)]["model_ligand"].qualified_name, 'J.G3D1')
        self.assertEqual(sc.rmsd_details["F"][mol.ResNum(1)]["chain_mapping"], {'B': 'B', 'G': 'A'})
        self.assertEqual(len(sc.rmsd_details["F"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(len(sc.rmsd_details["F"][mol.ResNum(1)]["bs_ref_res_mapped"]), 15)
        self.assertEqual(len(sc.rmsd_details["F"][mol.ResNum(1)]["bs_mdl_res_mapped"]), 15)
        self.assertEqual(sc.rmsd_details["F"][mol.ResNum(1)]["target_ligand"].qualified_name, 'K.G3D1')
        self.assertEqual(sc.rmsd_details["F"][mol.ResNum(1)]["model_ligand"].qualified_name, 'F.G3D1')

        # lddt_pli
        self.assertAlmostEqual(sc.lddt_pli["J"][mol.ResNum(1)], 0.9127105666156202, 5)
        self.assertAlmostEqual(sc.lddt_pli["F"][mol.ResNum(1)], 0.915929203539823, 6)
        # lddt_pli_details
        self.assertAlmostEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["rmsd"], 0.8016608357429504, 4)
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["lddt_pli_n_contacts"], 5224)
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["chain_mapping"], {'F': 'D', 'C': 'C'})
        self.assertEqual(len(sc.lddt_pli_details["J"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(len(sc.lddt_pli_details["J"][mol.ResNum(1)]["bs_ref_res_mapped"]), 15)
        self.assertEqual(len(sc.lddt_pli_details["J"][mol.ResNum(1)]["bs_mdl_res_mapped"]), 15)
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["target_ligand"].qualified_name, 'I.G3D1')
        self.assertEqual(sc.lddt_pli_details["J"][mol.ResNum(1)]["model_ligand"].qualified_name, 'J.G3D1')
        self.assertAlmostEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["rmsd"], 0.9286373257637024, 4)
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["lddt_pli_n_contacts"], 5424)
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["chain_mapping"], {'B': 'B', 'G': 'A'})
        self.assertEqual(len(sc.lddt_pli_details["F"][mol.ResNum(1)]["bs_ref_res"]), 15)
        self.assertEqual(len(sc.lddt_pli_details["F"][mol.ResNum(1)]["bs_ref_res_mapped"]), 15)
        self.assertEqual(len(sc.lddt_pli_details["F"][mol.ResNum(1)]["bs_mdl_res_mapped"]), 15)
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["target_ligand"].qualified_name, 'K.G3D1')
        self.assertEqual(sc.lddt_pli_details["F"][mol.ResNum(1)]["model_ligand"].qualified_name, 'F.G3D1')

    def test_global_chain_mapping(self):
        """Test that the global and local chain mappings works.

        For RMSD, A: A results in a better chain mapping. However, C: A is a
        better global chain mapping from an lDDT perspective (and lDDT-PLI).
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        # Local by default
        sc = LigandScorer(mdl, trg, None, None)
        assert sc.rmsd_details["L_2"][1]["chain_mapping"] == {'A': 'A'}
        assert sc.lddt_pli_details["L_2"][1]["chain_mapping"] == {'C': 'A'}

        # Global
        sc = LigandScorer(mdl, trg, None, None, global_chain_mapping=True)
        assert sc.rmsd_details["L_2"][1]["chain_mapping"] == {'C': 'A'}
        assert sc.lddt_pli_details["L_2"][1]["chain_mapping"] == {'C': 'A'}

    def test_rmsd_assignment(self):
        """Test that the RMSD-based assignment works.

        For RMSD, A: A results in a better chain mapping. However, C: A is a
        better global chain mapping from an lDDT perspective (and lDDT-PLI).
        """
        trg = _LoadMMCIF("1r8q.cif.gz")
        mdl = _LoadMMCIF("P84080_model_02.cif.gz")

        # By default, assignment differs between RMSD and lDDT-PLI in this
        # specific test case, so we can first ensure it does.
        # For now we skip as this is slow
        # sc = LigandScorer(mdl, trg, None, None)
        # assert sc.rmsd_details["L_2"][1]["target_ligand"] != sc.lddt_pli_details["L_2"][1]["target_ligand"]

        # RMSD assignment forces the same assignment
        sc = LigandScorer(mdl, trg, None, None, rmsd_assignment=True)
        assert sc.rmsd_details["L_2"][1]["target_ligand"] == sc.lddt_pli_details["L_2"][1]["target_ligand"]


if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_ligand_scoring.py tests.')
