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

    def test_extract_ligands(self):
        """Test that we can extract ligands from mmCIF files.
        """

        trg = io.LoadMMCIF(os.path.join('testfiles', "1r8q.cif.gz"))
        mdl = io.LoadMMCIF(os.path.join('testfiles', "P84080_model_02.cif.gz"))

        sc = LigandScorer(mdl, trg, None, None)

        assert len(sc.target_ligands) == 7
        assert len(sc.model_ligands) == 1


if __name__ == "__main__":
    from ost import testutils
    if testutils.SetDefaultCompoundLib():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_chain_mapping.py tests.')
