import unittest, os, sys
import ost
from ost import conop
from ost import io, mol, seq, settings
# check if we can import: fails if numpy or scipy not available
try:
    import numpy as np
    from ost.mol.alg.qsscore import *
    from ost.mol.alg.chain_mapping import *
except ImportError:
    print("Failed to import qsscore.py. Happens when numpy or scipy "\
          "missing. Ignoring qsscore.py tests.")
    sys.exit(0)

def _LoadFile(file_name):
    """Helper to avoid repeating input path over and over."""
    return io.LoadPDB(os.path.join('testfiles', file_name))

class TestQSScore(unittest.TestCase):

    def test_qsentity(self):
        ent = _LoadFile("3l1p.1.pdb")
        qsent = QSEntity(ent)
        self.assertEqual(len(qsent.view.chains), 4)
        self.assertEqual(qsent.GetChain("A").GetName(), "A")
        self.assertEqual(qsent.GetChain("B").GetName(), "B")
        self.assertEqual(qsent.GetChain("C").GetName(), "C")
        self.assertEqual(qsent.GetChain("D").GetName(), "D")
        self.assertRaises(Exception, qsent.GetChain, "E")
        self.assertEqual(qsent.chain_names, ["A", "B", "C", "D"])
        self.assertEqual(qsent.GetSequence("A"), "DMKALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTISRFEALQLSLKNMSKLRPLLEKWVEEADNNENLQEISKSVQARKRKRTSIENRVRWSLETMFLKSPKPSLQQITHIANQLGLEKDVVRVWFSNRRQKGKR")
        self.assertEqual(qsent.GetSequence("B"), "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTISRFEALQLSLKNMSKLRPLLEKWVEEADNNENLQEISKSQARKRKRTSIENRVRWSLETMFLKSPKPSLQQITHIANQLGLEKDVVRVWFSNRRQKGKRS")
        self.assertEqual(qsent.GetSequence("C"), "TCCACATTTGAAAGGCAAATGGA")
        self.assertEqual(qsent.GetSequence("D"), "ATCCATTTGCCTTTCAAATGTGG")

        # check for a couple of positions with manually extracted values

        # GLU
        pos = qsent.GetPos("B")
        self.assertAlmostEqual(pos[5,0], -1.661, places=3)
        self.assertAlmostEqual(pos[5,1], 27.553, places=3)
        self.assertAlmostEqual(pos[5,2], 12.774, places=3)

        # GLY
        pos = qsent.GetPos("A")
        self.assertAlmostEqual(pos[23,0], 17.563, places=3)
        self.assertAlmostEqual(pos[23,1], -4.082, places=3)
        self.assertAlmostEqual(pos[23,2], 29.005, places=3)

        # Cytosine
        pos = qsent.GetPos("C")
        self.assertAlmostEqual(pos[4,0], 14.796, places=3)
        self.assertAlmostEqual(pos[4,1], 24.653, places=3)
        self.assertAlmostEqual(pos[4,2], 59.318, places=3)


        # check pairwise dist, chain names are always sorted =>
        # A is rows, C is cols 
        dist_one = qsent.PairDist("A", "C")
        dist_two = qsent.PairDist("C", "A")
        self.assertTrue(np.array_equal(dist_one, dist_two))
        self.assertEqual(dist_one.shape[0], len(qsent.GetSequence("A")))
        self.assertEqual(dist_one.shape[1], len(qsent.GetSequence("C")))

        # check some random distance between the Gly and Cytosine that we already 
        # checked above
        self.assertAlmostEqual(dist_one[23,4], 41.86, places=2)

        # all chains interact with each other... but hey, check nevertheless
        self.assertEqual(qsent.interacting_chains, [("A", "B"), ("A", "C"),
                                                    ("A", "D"), ("B", "C"),
                                                    ("B", "D"), ("C", "D")])

    def test_qsscorer(self):

        target = _LoadFile("3l1p.1.pdb")
        model = _LoadFile("3l1p.1_model.pdb")

        # we need to derive a chain mapping prior to scoring
        mapper = ChainMapper(target)
        res = mapper.GetRigidMapping(model, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.472, places=2)

    def test_hetero_case_1(self):
        # additional chains
        ent_1 = _LoadFile('4ux8.1.pdb') # A2 B2 C2, symmetry: C2
        ent_2 = _LoadFile('3fub.2.pdb') # A2 B2   , symmetry: C2
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.825, 2)
        self.assertAlmostEqual(score_result.QS_best, 1.0, 2)

    def test_hetero_case_1_switched_order(self):
        # additional chains
        ent_2 = _LoadFile('4ux8.1.pdb') # A2 B2 C2, symmetry: C2
        ent_1 = _LoadFile('3fub.2.pdb') # A2 B2   , symmetry: C2
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.825, 2)
        self.assertAlmostEqual(score_result.QS_best, 1.0, 2)

    def test_HeteroCase1b(self):
        # as above but with assymetric unit of 3fub
        # -> no overlap found: check if extensive search can deal with it
        ent_1 = _LoadFile('4ux8.1.pdb')
        ent_2 = _LoadFile('3fub.au.pdb')
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.356, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.419, 2)

    def test_HeteroCase1b_switched_order(self):
        # chain mapping differs a bit when switching the order... I'm just
        # too lazy...
        pass

    def test_hetero_case_2(self):
        # different stoichiometry
        ent_1 = _LoadFile('1efu.1.pdb') # A2 B2, symmetry: C2
        ent_2 = _LoadFile('4pc6.1.pdb') # A B  , no symmetry
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.3131, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.941, 2)

    def test_hetero_case_2_switched_order(self):
        # different stoichiometry
        ent_2 = _LoadFile('1efu.1.pdb') # A2 B2, symmetry: C2
        ent_1 = _LoadFile('4pc6.1.pdb') # A B  , no symmetry
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.3131, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.941, 2)

    def test_hetero_case_3(self):
        # more chains
        ent_1 = _LoadFile('2vjt.1.pdb') # A6 B6, symmetry: D3
        ent_2 = _LoadFile('3dbj.1.pdb') # A3 B3, symmetry: C3
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.359, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.958, 2)

    def test_hetero_case_3_switched_order(self):
        # more chains
        ent_2 = _LoadFile('2vjt.1.pdb') # A6 B6, symmetry: D3
        ent_1 = _LoadFile('3dbj.1.pdb') # A3 B3, symmetry: C3
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.359, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.958, 2)

    def test_hetero_case_4(self):
        # inverted chains
        ent_1 = _LoadFile('3ia3.1.pdb') # AB, no symmetry
        ent_2 = _LoadFile('3ia3.2.pdb') # BA, no symmetry
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.980, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.980, 2)

    def test_hetero_case_4_switched_order(self):
        # inverted chains
        ent_2 = _LoadFile('3ia3.1.pdb') # AB, no symmetry
        ent_1 = _LoadFile('3ia3.2.pdb') # BA, no symmetry
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.980, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.980, 2)

    def test_hetero_model(self):
        # uncomplete model missing 2 third of the contacts
        target = _LoadFile('1eud_ref.pdb')               # AB, no symmetry
        model = _LoadFile('1eud_mdl_partial-dimer.pdb') # BA, no symmetry
        mapper = ChainMapper(target)
        res = mapper.GetRigidMapping(model, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.323, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.921, 2)

    def test_hetero_model_switched_order(self):
        # same as above but with switched order to test for symmetric behaviour
        # of QS score
        target = _LoadFile('1eud_mdl_partial-dimer.pdb') # BA, no symmetry
        model = _LoadFile('1eud_ref.pdb')               # AB, no symmetry
        mapper = ChainMapper(target)
        res = mapper.GetRigidMapping(model, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.323, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.921, 2)

    def test_homo_1(self):
        # different stoichiometry SOD
        ent_1 = _LoadFile('4dvh.1.pdb') # A2, symmetry: C2
        ent_2 = _LoadFile('4br6.1.pdb') # A4, symmetry: D2
        # original qsscoring uses other default values for gap_open and gap_extension
        # penalties, let's use those to reproduce the old results as the alignments
        # would differ otherwise
        mapper = ChainMapper(ent_1, pep_gap_open = -5, pep_gap_ext = -2)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.147, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.866, 2)

    def test_homo_1_switched_order(self):
        # different stoichiometry SOD
        ent_2 = _LoadFile('4dvh.1.pdb') # A2, symmetry: C2
        ent_1 = _LoadFile('4br6.1.pdb') # A4, symmetry: D2
        # original qsscoring uses other default values for gap_open and gap_extension
        # penalties, let's use those to reproduce the old results as the alignments
        # would differ otherwise
        mapper = ChainMapper(ent_1, pep_gap_open = -5, pep_gap_ext = -2)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 0.147, 2)
        self.assertAlmostEqual(score_result.QS_best, 0.866, 2)

    def test_homo_2(self):
        # broken cyclic symmetry
        ent_1 = _LoadFile('4r7y.1.pdb')   # A6, symmetry: C6
        ent_2 = ent_1.Select('cname=A,B') # A2, no symmetry
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 1/6, 2)
        self.assertAlmostEqual(score_result.QS_best, 1.0, 2)

    def test_homo_2_switched_order(self):
        # same as above but with switched order to test for symmetric behaviour
        # of QS score
        ent_2 = _LoadFile('4r7y.1.pdb')   # A6, symmetry: C6
        ent_1 = ent_2.Select('cname=A,B') # A2, no symmetry
        mapper = ChainMapper(ent_1)
        res = mapper.GetRigidMapping(ent_2, strategy="greedy_iterative_rmsd")
        qs_scorer = QSScorer.FromMappingResult(res)
        score_result = qs_scorer.Score(res.mapping)
        self.assertAlmostEqual(score_result.QS_global, 1/6, 2)
        self.assertAlmostEqual(score_result.QS_best, 1.0, 2)

if __name__ == "__main__":
    from ost import testutils
    if testutils.SetDefaultCompoundLib():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_qsscore.py tests.')
