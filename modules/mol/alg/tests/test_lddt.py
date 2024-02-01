import unittest, os, sys
import ost
from ost import io, mol, settings, conop, seq
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg.lddt import *
    from ost.mol.alg.scoring import *
except ImportError:
    print("Failed to import qsscoring. Happens when numpy or scipy missing. " \
          "Ignoring test_lddt.py tests.")
    sys.exit(0)

def _LoadFile(file_name):
    """Helper to avoid repeating input path over and over."""
    return io.LoadPDB(os.path.join('testfiles', file_name))


class TestlDDT(unittest.TestCase):

    # compare monomers to lDDT C++ reference implementation
    def test_lDDT_monomer(self):

        # do 7SGN
        model = _LoadFile("7SGN_C_model.pdb")
        target = _LoadFile("7SGN_C_target.pdb")

        # do awesome implementation
        scorer = lDDTScorer(target)
        aws_score, aws_per_res_scores = scorer.lDDT(model)

        # do reference implementation
        dl = mol.alg.CreateDistanceList(target.CreateFullView(), 15.0)
        classic_score = mol.alg.LDDTHA(model.CreateFullView(), dl)
        classic_per_res_scores = list()
        for r in model.residues:
            if r.HasProp("locallddt"):
                classic_per_res_scores.append(r.GetFloatProp("locallddt"))
            else:
                classic_per_res_scores.append(None)

        self.assertAlmostEqual(aws_score, classic_score, places=5)

        self.assertEqual(len(aws_per_res_scores), len(classic_per_res_scores))
        for a,b in zip(aws_per_res_scores, classic_per_res_scores):
            if a is None and b is None:
                continue
            self.assertAlmostEqual(a, b, places = 5)

        # do 7W1F_B
        model = _LoadFile("7W1F_B_model.pdb")
        target = _LoadFile("7W1F_B_target.pdb")

        # do awesome implementation
        scorer = lDDTScorer(target)
        aws_score, aws_per_res_scores = scorer.lDDT(model)

        # do reference implementation
        dl = mol.alg.CreateDistanceList(target.CreateFullView(), 15.0)
        classic_score = mol.alg.LDDTHA(model.CreateFullView(), dl)
        classic_per_res_scores = list()
        for r in model.residues:
            if r.HasProp("locallddt"):
                classic_per_res_scores.append(r.GetFloatProp("locallddt"))
            else:
                classic_per_res_scores.append(None)

        self.assertAlmostEqual(aws_score, classic_score, places=5)

        self.assertEqual(len(aws_per_res_scores), len(classic_per_res_scores))
        for a,b in zip(aws_per_res_scores, classic_per_res_scores):
            if a is None and b is None:
                continue
            self.assertAlmostEqual(a, b, places = 5)

    # check oligo functionality
    def test_lDDT_oligo(self):

        ent_full = _LoadFile("4br6.1.pdb")
        model = ent_full.Select('peptide=true')
        target = ent_full.Select('peptide=true and cname=A,B')
        # hardcoded chain mapping
        chain_mapping = {"A": "A", "B": "B"}
        lddt_scorer = lDDTScorer(target)
        
        score, per_res_scores = lddt_scorer.lDDT(model, 
          chain_mapping=chain_mapping)
        self.assertAlmostEqual(score, 1.0, places=5)

        score, per_res_scores = lddt_scorer.lDDT(model, 
          chain_mapping=chain_mapping, no_interchain=True)
        self.assertAlmostEqual(score, 1.0, places=5)

        score, per_res_scores = lddt_scorer.lDDT(model, 
          chain_mapping=chain_mapping, no_interchain=False,
          penalize_extra_chains=True)
        self.assertAlmostEqual(score, 0.52084655, places=5)

        score, per_res_scores = lddt_scorer.lDDT(model, 
          chain_mapping=chain_mapping, no_interchain=True,
          penalize_extra_chains=True)
        self.assertAlmostEqual(score, 0.499570048, places=5)

    def test_lDDT_custom_resmapping(self):

        ent_full = _LoadFile("4br6.1.pdb")
        model = ent_full.Copy().Select('peptide=true')
        target = ent_full.Select('peptide=true and cname=A,B')

        # shift residue numbers in model
        ed = model.handle.EditXCS()
        for ch in model.chains:
            ed.RenumberChain(ch.handle, 42, True)

        # hardcoded chain mapping
        chain_mapping = {"A": "A", "B": "B"}
        lddt_scorer = lDDTScorer(target)

        # naively running lDDT will fail, as residue-residue mapping happens
        # with resnums. Since we shifted that stuff above we'll get an error
        # complaining about residue name mismatch
        with self.assertRaises(RuntimeError):
            score, per_res_scores = lddt_scorer.lDDT(model, 
              chain_mapping=chain_mapping, no_interchain=False,
              penalize_extra_chains=True)

        # we can rescue that with alignments
        res_map = dict()
        for mdl_ch_name, trg_ch_name in chain_mapping.items():
            mdl_ch = model.FindChain(mdl_ch_name)
            trg_ch = target.FindChain(trg_ch_name)
            mdl_seq = ''.join([r.one_letter_code for r in mdl_ch.residues])
            mdl_seq = seq.CreateSequence(mdl_ch_name, mdl_seq)
            trg_seq = ''.join([r.one_letter_code for r in trg_ch.residues])
            trg_seq = seq.CreateSequence(trg_ch_name, trg_seq)
            aln = seq.alg.GlobalAlign(trg_seq, mdl_seq, seq.alg.BLOSUM62)[0]
            res_map[mdl_ch_name] = aln

        score, per_res_scores = lddt_scorer.lDDT(model, 
              chain_mapping=chain_mapping, no_interchain=False,
              penalize_extra_chains=True, residue_mapping=res_map)
        self.assertAlmostEqual(score, 0.52084655, places=5)

    def test_lDDT_seqsep(self):
        target = _LoadFile("7SGN_C_target.pdb")
        with self.assertRaises(NotImplementedError):
            scorer = lDDTScorer(target, sequence_separation=42)
        scorer = lDDTScorer(target, sequence_separation=0)

    def test_bb_only(self):
        model = _LoadFile("7SGN_C_model.pdb")
        target = _LoadFile("7SGN_C_target.pdb")

        # do scoring and select aname=CA
        scorer = lDDTScorer(target.Select("aname=CA"))
        score_one, per_res_scores_one = scorer.lDDT(model)
        score_two, per_res_scores_two = scorer.lDDT(model.Select("aname=CA"))

        # no selection, just setting bb_only flag should give the same
        scorer = lDDTScorer(target, bb_only=True)
        score_three, per_res_scores_three = scorer.lDDT(model)

        # check
        self.assertAlmostEqual(score_one, score_two, places=5)
        self.assertAlmostEqual(score_one, score_three, places=5)
        for a,b in zip(per_res_scores_one, per_res_scores_two):
            self.assertAlmostEqual(a, b, places=5)
        for a,b in zip(per_res_scores_one, per_res_scores_three):
            self.assertAlmostEqual(a, b, places=5)

    def test_resname_match(self):
        model = _LoadFile("7SGN_C_model.pdb")
        target = _LoadFile("7SGN_C_target.pdb")

        # introduce name mismatch
        ed = model.handle.EditXCS()
        ed.RenameResidue(model.residues[42], "asdf")

        # do scoring and select aname=CA
        scorer = lDDTScorer(target.Select("aname=CA"))

        with self.assertRaises(RuntimeError):
            scorer.lDDT(model)

        scorer.lDDT(model, check_resnames=False)

    def test_intra_interchain(self):
        ent_full = _LoadFile("4br6.1.pdb")
        model = ent_full.Select('peptide=true and cname=A,B')
        target = ent_full.Select('peptide=true and cname=A,B')
        chain_mapping = {"A": "A", "B": "B"}

        lddt_scorer = lDDTScorer(target)

        # do lDDT only on interchain contacts (ic)
        lDDT_ic, per_res_lDDT_ic, lDDT_tot_ic, lDDT_cons_ic, \
        res_indices_ic, per_res_exp_ic, per_res_conserved_ic =\
        lddt_scorer.lDDT(model, no_intrachain=True, 
                         chain_mapping = chain_mapping,
                         return_dist_test = True)

        # do lDDT only on intrachain contacts (sc for single chain)
        lDDT_sc, per_res_lDDT_sc, lDDT_tot_sc, lDDT_cons_sc, \
        res_indices_sc, per_res_exp_sc, per_res_conserved_sc =\
        lddt_scorer.lDDT(model, no_interchain=True,
                         chain_mapping = chain_mapping,
                         return_dist_test = True)

        # do lDDT on everything
        lDDT, per_res_lDDT, lDDT_tot, lDDT_cons, res_indices, per_res_exp, \
        per_res_conserved = lddt_scorer.lDDT(model,
                                             chain_mapping = chain_mapping,
                                             return_dist_test = True)

        # sum of lDDT_tot_ic and lDDT_tot_sc should be equal to lDDT_tot
        self.assertEqual(lDDT_tot_ic + lDDT_tot_sc, lDDT_tot)

        # same for the conserved contacts
        self.assertEqual(lDDT_cons_ic + lDDT_cons_sc, lDDT_cons)


class TestlDDTBS(unittest.TestCase):

    def test_basic(self):
        mdl = _LoadFile("lddtbs_mdl.pdb")
        ref = _LoadFile("lddtbs_ref_1r8q.1.pdb")

        lddtbs_scorer = lDDTBSScorer(reference=ref, model=mdl)
        bs_repr = lddtbs_scorer.ScoreBS(ref.Select("rname=AFB"), radius = 5.0,
                                        lddt_radius = 12.0)

        # select residues manually from reference
        for at in ref.Select("rname=AFB").atoms:
            close_atoms = ref.FindWithin(at.GetPos(), 5.0)
            for close_at in close_atoms:
                close_at.GetResidue().SetIntProp("asdf", 1)

        ref_bs = ref.Select("grasdf:0=1")
        ref_bs = ref_bs.Select("peptide=true")
        ref_bs_names = [r.GetQualifiedName() for r in ref_bs.residues]
        repr_bs_names = [r.GetQualifiedName() for r in bs_repr.ref_residues]
        self.assertEqual(sorted(ref_bs_names), sorted(repr_bs_names))


        # everything below basically computes lDDTBS manually and
        # compares with the result we got above 

        # select residues manually from model
        fancy_mapping = {"A":"B", "B":"A"} # hardcoded chain mapping...
        mdl_bs = mdl.CreateEmptyView()
        for r in ref_bs.residues:
            mdl_res = mdl.FindResidue(fancy_mapping[r.GetChain().GetName()],
                                      r.GetNumber())
            mdl_bs.AddResidue(mdl_res, mol.ViewAddFlag.INCLUDE_ALL)

        # put that stuff in single chain structures
        sc_ref_bs = mol.CreateEntity()
        ed = sc_ref_bs.EditXCS()
        ch = ed.InsertChain("A")
        for r in ref_bs.residues:
            added_r = ed.AppendResidue(ch, r.GetName())
            for a in r.atoms:
                ed.InsertAtom(added_r, a.GetName(), a.GetPos())

        sc_mdl_bs = mol.CreateEntity()
        ed = sc_mdl_bs.EditXCS()
        ch = ed.InsertChain("A")
        for r in mdl_bs.residues:
            added_r = ed.AppendResidue(ch, r.GetName())
            for a in r.atoms:
                ed.InsertAtom(added_r, a.GetName(), a.GetPos())

        # compute and compare
        lddt_scorer = lDDTScorer(sc_ref_bs, inclusion_radius=12.0)
        self.assertAlmostEqual(bs_repr.lDDT, lddt_scorer.lDDT(sc_mdl_bs)[0])


if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound library available. Ignoring test_lddt.py tests.')
