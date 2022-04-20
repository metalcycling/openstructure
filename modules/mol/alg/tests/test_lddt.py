import unittest, os, sys
import ost
from ost import io, mol, settings, conop, seq
# check if we can import: fails if numpy or scipy not available
try:
    from ost.mol.alg.qsscoring import *
    from ost.mol.alg.lddt import *
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
        # we use functionality from QS-scorer to derive a mapping
        qs_scorer = QSscorer(model, target)
        lddt_scorer = lDDTScorer(target)

        score, per_res_scores = lddt_scorer.lDDT(model, 
          chain_mapping=qs_scorer.chain_mapping)
        self.assertAlmostEqual(score, 1.0, places=5)

        score, per_res_scores = lddt_scorer.lDDT(model, 
          chain_mapping=qs_scorer.chain_mapping, no_interchain=True)
        self.assertAlmostEqual(score, 1.0, places=5)

        score, per_res_scores = lddt_scorer.lDDT(model, 
          chain_mapping=qs_scorer.chain_mapping, no_interchain=False,
          penalize_extra_chains=True)
        self.assertAlmostEqual(score, 0.52084655, places=5)

        score, per_res_scores = lddt_scorer.lDDT(model, 
          chain_mapping=qs_scorer.chain_mapping, no_interchain=True,
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

        # we use functionality from QS-scorer to derive a mapping
        qs_scorer = QSscorer(model, target)
        lddt_scorer = lDDTScorer(target)

        # naively running lDDT will fail, as residue-residue mapping happens
        # with resnums. Since we shifted that stuff above we'll get an error
        # complaining about residue name mismatch
        with self.assertRaises(RuntimeError):
            score, per_res_scores = lddt_scorer.lDDT(model, 
              chain_mapping=qs_scorer.chain_mapping, no_interchain=False,
              penalize_extra_chains=True)

        # we can rescue that with alignments from qsscorer
        res_map = dict()
        for aln in qs_scorer.alignments:
            model_chain_name = aln.GetSequence(0).GetName()
            # we need to inverse the direction... qsscorer
            # has first model sequence and then target sequence
            # (at least the way we set it up above...)
            new_aln = seq.CreateAlignment()
            new_aln.AddSequence(aln.GetSequence(1))
            new_aln.AddSequence(aln.GetSequence(0))
            res_map[model_chain_name] = new_aln

        score, per_res_scores = lddt_scorer.lDDT(model, 
              chain_mapping=qs_scorer.chain_mapping, no_interchain=False,
              penalize_extra_chains=True, residue_mapping=res_map)
        self.assertAlmostEqual(score, 0.52084655, places=5)

    def test_lDDT_seqsep(self):
        target = _LoadFile("7SGN_C_target.pdb")
        with self.assertRaises(NotImplementedError):
            scorer = lDDTScorer(target, sequence_separation=42)
        scorer = lDDTScorer(target, sequence_separation=0)

    def test_calpha(self):
        model = _LoadFile("7SGN_C_model.pdb")
        target = _LoadFile("7SGN_C_target.pdb")

        # do scoring and select aname=CA
        scorer = lDDTScorer(target.Select("aname=CA"))
        score_one, per_res_scores_one = scorer.lDDT(model)
        score_two, per_res_scores_two = scorer.lDDT(model.Select("aname=CA"))

        # no selection, just setting calpha flag should give the same
        scorer = lDDTScorer(target, calpha=True)
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

if __name__ == "__main__":
    from ost import testutils
    if testutils.SetDefaultCompoundLib():
        testutils.RunTests()
    else:
        print('No compound library available. Ignoring test_lddt.py tests.')
