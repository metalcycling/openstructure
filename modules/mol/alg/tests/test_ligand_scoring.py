import unittest, os, sys

from ost import io, mol
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg.ligand_scoring import *
except ImportError:
    print("Failed to import ligand_scoring.py. Happens when numpy or scipy " \
          "missing. Ignoring test_ligand_scoring.py tests.")
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



if __name__ == "__main__":
    from ost import testutils
    if testutils.SetDefaultCompoundLib():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_ligand_scoring.py tests.')
