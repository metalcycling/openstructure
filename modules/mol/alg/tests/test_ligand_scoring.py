import unittest, os, sys

from ost import io, mol, geom
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg.ligand_scoring import *
    import ost.mol.alg.ligand_scoring
except ImportError:
    print("Failed to import ligand_scoring.py. Happens when numpy, scipy or "
          "networkx is missing. Ignoring test_ligand_scoring.py tests.")
    sys.exit(0)


class TestLigandScoring(unittest.TestCase):

    def test_extract_ligands_mmCIF(self):
        """Test that we can extract ligands from mmCIF files.
        """
        trg, trg_seqres = io.LoadMMCIF(os.path.join('testfiles', "1r8q.cif.gz"), seqres=True)
        mdl, mdl_seqres = io.LoadMMCIF(os.path.join('testfiles', "P84080_model_02.cif.gz"), seqres=True)

        sc = LigandScorer(mdl, trg, None, None)

        assert len(sc.target_ligands) == 7
        assert len(sc.model_ligands) == 1
        assert len([r for r in sc.target.residues if r.is_ligand]) == 7
        assert len([r for r in sc.model.residues if r.is_ligand]) == 1

    def test_init_given_ligands(self):
        """Test that we can instantiate the scorer with ligands contained in
        the target and model entity and given in a list.
        """
        trg, trg_seqres = io.LoadMMCIF(os.path.join('testfiles', "1r8q.cif.gz"), seqres=True)
        mdl, mdl_seqres = io.LoadMMCIF(os.path.join('testfiles', "P84080_model_02.cif.gz"), seqres=True)

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
        mdl = io.LoadPDB(os.path.join('testfiles', "P84080_model_02_nolig.pdb"))
        mdl_ligs = [io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))]
        trg = io.LoadPDB(os.path.join('testfiles', "1r8q_protein.pdb.gz"))
        trg_ligs = [io.LoadEntity(os.path.join('testfiles', "1r8q_ligand_%d.sdf" % i)) for i in range(7)]

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
        mdl = io.LoadPDB(os.path.join('testfiles', "P84080_model_02_nolig.pdb"))
        mdl_ligs = [io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))]
        trg = io.LoadPDB(os.path.join('testfiles', "1r8q_protein.pdb.gz"))
        trg_ligs = [io.LoadEntity(os.path.join('testfiles', "1r8q_ligand_%d.sdf" % i)) for i in range(7)]

        # Reject identical model ligands
        with self.assertRaises(RuntimeError):
            sc = LigandScorer(mdl, trg, [mdl_ligs[0], mdl_ligs[0]], trg_ligs)

        # Reject identical target ligands
        lig0 = trg_ligs[0]
        lig1 = trg_ligs[1]
        ed1 = lig1.EditXCS()
        ed1.RenameChain(lig1.chains[0], lig0.chains[0].name)
        ed1.SetResidueNumber(lig1.residues[0], lig0.residues[0].number)
        with self.assertRaises(RuntimeError):
            sc = LigandScorer(mdl, trg, mdl_ligs, [lig0, lig1])

    def test_ResidueToGraph(self):
        """Test that ResidueToGraph works as expected
        """
        mdl_lig = io.LoadEntity(os.path.join('testfiles', "P84080_model_02_ligand_0.sdf"))

        graph = ResidueToGraph(mdl_lig.residues[0])
        assert len(graph.edges) == 34
        assert len(graph.nodes) == 32
        # Check an arbitrary node
        assert [a for a in graph.adj["14"].keys()] == ["13", "29"]

    def test__ComputeSymmetries(self):
        """Test that _ComputeSymmetries works.
        """
        trg, trg_seqres = io.LoadMMCIF(os.path.join('testfiles', "1r8q.cif.gz"), seqres=True)
        mdl, mdl_seqres = io.LoadMMCIF(os.path.join('testfiles', "P84080_model_02.cif.gz"), seqres=True)

        trg_mg1 = trg.FindResidue("E", 1)
        trg_g3d1 = trg.FindResidue("F", 1)
        trg_afb1 = trg.FindResidue("G", 1)
        trg_g3d2 = trg.FindResidue("J", 1)
        mdl_g3d = mdl.FindResidue("L_2", 1)

        sym = ost.mol.alg.ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1)
        assert len(sym) == 72

        sym = ost.mol.alg.ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1, by_atom_index=True)
        assert len(sym) == 72

        # Test that it works with views and only consider atoms in the view
        # Skip PA, PB and O[1-3]A and O[1-3]B in target and model
        # We assume atom index are fixed and won't change
        trg_g3d1_sub = trg_g3d1.Select("aindex>6019").residues[0]
        mdl_g3d_sub = mdl_g3d.Select("aindex>1447").residues[0]

        sym = ost.mol.alg.ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1_sub)
        assert len(sym) == 6

        sym = ost.mol.alg.ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1_sub, by_atom_index=True)
        assert len(sym) == 6

        # Substructure matches
        self.skipTest("Substructure matches don't work yet")
        sym = ost.mol.alg.ligand_scoring._ComputeSymmetries(mdl_g3d, trg_g3d1_sub, substructure_match=True)
        assert len(sym) == 6

        # Missing atoms only allowed in target, not in model
        with self.assertRaises(NoSymmetryError):
            ost.mol.alg.ligand_scoring._ComputeSymmetries(mdl_g3d_sub, trg_g3d1, substructure_match=True)

    def test_LigandRMSD(self):
        """Test that LigandRMSD works.
        """
        trg, trg_seqres = io.LoadMMCIF(os.path.join('testfiles', "1r8q.cif.gz"), seqres=True)
        mdl, mdl_seqres = io.LoadMMCIF(os.path.join('testfiles', "P84080_model_02.cif.gz"), seqres=True)

        trg_mg1 = trg.FindResidue("E", 1)
        trg_g3d1 = trg.FindResidue("F", 1)
        trg_afb1 = trg.FindResidue("G", 1)
        trg_g3d2 = trg.FindResidue("J", 1)
        mdl_g3d = mdl.FindResidue("L_2", 1)

        rmsd = LigandRMSD(mdl_g3d, trg_g3d1)
        self.assertAlmostEqual(rmsd, 2.21341e-06, 10)
        rmsd = LigandRMSD(mdl_g3d, trg_g3d2)
        self.assertAlmostEqual(rmsd, 61.21325, 4)

        # Ensure we raise a NoSymmetryError if the ligand is wrong
        with self.assertRaises(NoSymmetryError):
            LigandRMSD(mdl_g3d, trg_mg1)
        with self.assertRaises(NoSymmetryError):
            LigandRMSD(mdl_g3d, trg_afb1)

        # Assert that transform works
        trans = geom.Mat4(-0.999256, 0.00788487, -0.0377333, -15.4397,
                          0.0380652, 0.0473315, -0.998154, 29.9477,
                          -0.00608426, -0.998848, -0.0475963, 28.8251,
                          0, 0, 0, 1)
        rmsd = LigandRMSD(mdl_g3d, trg_g3d2, transformation=trans)
        self.assertAlmostEqual(rmsd, 0.293972, 5)

        # Assert that substructure matches work
        self.skipTest("Substructure matches don't work yet")
        trg_g3d1_sub = trg_g3d1.Select("aindex>6019").residues[0] # Skip PA, PB and O[1-3]A and O[1-3]B.
        mdl_g3d_sub = mdl_g3d.Select("aindex>1447").residues[0] # Skip PA, PB and O[1-3]A and O[1-3]B.
        with self.assertRaises(NoSymmetryError):
            LigandRMSD(mdl_g3d, trg_g3d1_sub)  # no full match
        rmsd = LigandRMSD(mdl_g3d, trg_g3d1_sub, substructure_match=True)


if __name__ == "__main__":
    from ost import testutils
    if testutils.SetDefaultCompoundLib():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_ligand_scoring.py tests.')
