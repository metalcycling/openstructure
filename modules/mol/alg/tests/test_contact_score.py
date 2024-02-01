import unittest, os, sys
import ost
from ost import conop
from ost import io, mol, seq, settings
# check if we can import: fails if numpy or scipy not available
try:
    import numpy as np
    from ost.mol.alg.contact_score import *
    from ost.mol.alg.chain_mapping import *
except ImportError:
    print("Failed to import contact_score.py. Happens when numpy or scipy "\
          "missing. Ignoring contact_score.py tests.")
    sys.exit(0)

def _LoadFile(file_name):
    """Helper to avoid repeating input path over and over."""
    return io.LoadPDB(os.path.join('testfiles', file_name))

class TestContactScore(unittest.TestCase):
    def test_ContactEntity(self):
        self.maxDiff = None
        ent = _LoadFile("3l1p.1.pdb")
        cent = ContactEntity(ent)
        self.assertEqual(cent.GetChain("A").GetName(), "A")
        self.assertEqual(cent.GetChain("B").GetName(), "B")
        self.assertEqual(cent.GetChain("C").GetName(), "C")
        self.assertEqual(cent.GetChain("D").GetName(), "D")
        self.assertRaises(Exception, cent.GetChain, "E")
        self.assertEqual(cent.chain_names, ["A", "B", "C", "D"])
        self.assertEqual(cent.GetSequence("A"), "DMKALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTISRFEALQLSLKNMSKLRPLLEKWVEEADNNENLQEISKSVQARKRKRTSIENRVRWSLETMFLKSPKPSLQQITHIANQLGLEKDVVRVWFSNRRQKGKR")
        self.assertEqual(cent.GetSequence("B"), "KALQKELEQFAKLLKQKRITLGYTQADVGLTLGVLFGKVFSQTTISRFEALQLSLKNMSKLRPLLEKWVEEADNNENLQEISKSQARKRKRTSIENRVRWSLETMFLKSPKPSLQQITHIANQLGLEKDVVRVWFSNRRQKGKRS")
        self.assertEqual(cent.GetSequence("C"), "TCCACATTTGAAAGGCAAATGGA")
        self.assertEqual(cent.GetSequence("D"), "ATCCATTTGCCTTTCAAATGTGG")
        self.assertEqual(cent.contact_mode, "aa")
        self.assertEqual(cent.contact_d, 5.0)
        self.assertEqual(cent.interacting_chains, [('A', 'B'), ('A', 'D'),
                                                   ('A', 'C'), ('B', 'C'),
                                                   ('B', 'D'), ('C', 'D')])
        exp_contacts = sorted(list(cent.contacts[('A', 'C')]))
        self.assertEqual(exp_contacts, [(40, 9), (41, 8), (41, 9), (42, 8),
                                        (42, 9), (42, 10), (43, 12), (44, 9),
                                        (44, 10), (44, 11), (45, 8), (45, 9),
                                        (48, 8), (48, 9), (54, 8), (55, 6),
                                        (55, 7), (57, 7), (58, 7), (58, 8),
                                        (62, 8), (91, 8), (91, 9), (91, 10),
                                        (93, 8), (93, 9), (93, 10), (95, 10),
                                        (95, 11), (113, 2), (113, 3), (115, 2),
                                        (134, 1), (139, 5), (141, 2), (141, 3),
                                        (142, 4), (142, 5), (142, 6), (145, 4)])

    def test_ContactScorer(self):
        target = _LoadFile("3l1p.1.pdb")
        model = _LoadFile("3l1p.1_model.pdb")

        # we need to derive a chain mapping prior to scoring
        mapper = ChainMapper(target)
        res = mapper.GetRigidMapping(model, strategy="greedy_iterative_rmsd")
        contact_scorer = ContactScorer.FromMappingResult(res)
        score_result = contact_scorer.ScoreICS(res.mapping)
        self.assertAlmostEqual(score_result.precision, 0.583, places=2)
        self.assertAlmostEqual(score_result.recall, 0.288, places=2)
        self.assertAlmostEqual(score_result.ics, 0.386, places=2)

        score_result = contact_scorer.ScoreIPS(res.mapping)
        self.assertAlmostEqual(score_result.precision, 0.779, places=2)
        self.assertAlmostEqual(score_result.recall, 0.493, places=2)
        self.assertAlmostEqual(score_result.ips, 0.432, places=2)

if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring contact_score.py tests.')




