import unittest
from ost import geom, mol


class TestInvalidHandleOrView(unittest.TestCase):
    """ Test some basic behavior of invalid (chain|residue|atom|bond)
    (view|handle)s. in Python. This test checks the behavior
    of valid and IsValid() as well as presence and behavior of hash_code
    and GetHashCode() on relevant handles."""

    def setUp(self):
        """ Setup some test data."""
        # Create an entity with 3 atoms and 1 bond.
        self.entity = mol.CreateEntity()
        ed = self.entity.EditXCS()
        ch = ed.InsertChain("A")
        res = ed.AppendResidue(ch, "ALA")
        self.at0 = ed.InsertAtom(res, "A", geom.Vec3(0, 0, 0))
        self.at1 = ed.InsertAtom(res, "B", geom.Vec3(1, 0, 0))
        self.at2 = ed.InsertAtom(res, "C", geom.Vec3(2, 0, 0))
        ed.Connect(self.at0, self.at1)
        ed.UpdateICS()
        self.entity_view = self.entity.Select("cname=A")

    def assertValidHandleOrView(self, x):
        """Assertion for valid handles: valid and IsValid() are True,
        hash_code and GetHashCode() work and return a positive integer.

        This test also ensures the hash code is > 0. Although we make no
        promise about the exact value of the hash, an earlier behavior was to
        return a hash_code of 0 (null pointer) for some invalid handles. If you
        get a hash code of 0, your installation of OpenStructure is most likely
        buggy.

        """
        self.assertTrue(x.valid)
        self.assertTrue(x.IsValid())
        self.assertIsInstance(x.GetHashCode(), int)
        self.assertTrue(x.hash_code, int)
        self.assertTrue(x.GetHashCode() > 0)
        self.assertTrue(x.hash_code > 0)

    def assertInvalidHandleOrView(self, x):
        """Assertion for invalid handles: valid and IsValid() are False,
        hash_code and GetHashCode() raise an InvalidHandle error (which
        is converted to an Exception in Python with a specific error
        message).
        """
        self.assertFalse(x.valid)
        self.assertFalse(x.IsValid())
        with self.assertRaises(Exception, msg="Can not access invalid handle or view"):
            x.hash_code
        with self.assertRaises(Exception, msg="Can not access invalid handle or view"):
            x.GetHashCode()

    def test_InvalidChainHandle(self):
        """ Test behavior of chain handles"""

        # Valid
        ch = self.entity.FindChain("A")
        self.assertValidHandleOrView(ch)

        # Invalid
        ch = self.entity.FindChain("X")
        self.assertInvalidHandleOrView(ch)

    def test_InvalidChainView(self):
        """ Test behavior of chain views"""

        # Valid
        ch = self.entity_view.FindChain("A")
        self.assertValidHandleOrView(ch)

        # Invalid
        ch = self.entity_view.FindChain("X")
        self.assertInvalidHandleOrView(ch)

    def test_InvalidResidueHandle(self):
        """ Test behavior of residue handles"""

        # Valid
        res = self.entity.FindResidue("A", 1)
        self.assertValidHandleOrView(res)

        # Invalid
        res = self.entity.FindResidue("X", 1)
        self.assertInvalidHandleOrView(res)

    def test_InvalidResidueView(self):
        """ Test behavior of residue views"""

        # Valid
        res = self.entity_view.FindResidue("A", 1)
        self.assertValidHandleOrView(res)

        # Invalid
        res = self.entity_view.FindResidue("X", 1)
        self.assertInvalidHandleOrView(res)

    def test_InvalidAtomHandle(self):
        """ Test behavior of atom handles"""

        # Valid
        at = self.entity.FindAtom("A", 1, "A")
        self.assertValidHandleOrView(at)

        # Invalid
        at = self.entity.FindAtom("X", 1, "A")
        self.assertInvalidHandleOrView(at)

    def test_InvalidAtomView(self):
        """ Test behavior of atom views"""

        # Valid
        at = self.entity_view.FindAtom("A", 1, "A")
        self.assertValidHandleOrView(at)

        # Invalid
        at = self.entity_view.FindAtom("X", 1, "A")
        self.assertInvalidHandleOrView(at)

    def test_InvalidBondHandle(self):
        """ Test behavior of atom handles"""

        # Valid
        bond = self.at0.FindBondToAtom(self.at1)
        self.assertValidHandleOrView(bond)

        # Invalid
        bond = self.at0.FindBondToAtom(self.at0)
        self.assertInvalidHandleOrView(bond)


if __name__ == "__main__":
    unittest.main()

