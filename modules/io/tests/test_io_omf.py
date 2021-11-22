import unittest
from ost import geom
from ost import io

def compare_atoms(a1, a2):
    if abs(a1.occupancy - a2.occupancy) > 0.01:
        return False
    if abs(a1.b_factor - a2.b_factor) > 0.01:
        return False
    if geom.Distance(a1.GetPos(), a2.GetPos()) > 0.001:
        return False
    if a1.is_hetatom != a2.is_hetatom:
        return False
    if a1.element != a2.element:
        return False
    return True

def compare_residues(r1, r2):
    if r1.GetName() != r2.GetName():
        return False
    if r1.GetNumber() != r2.GetNumber():
        return False
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
        if not compare_atoms(a1, a2):
            return False
    return True

def compare_chains(ch1, ch2):
    if len(ch1.residues) != len(ch2.residues):
        return False
    for r1, r2 in zip(ch1.residues, ch2.residues):
        if not compare_residues(r1, r2):
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

def compare_ent(ent1, ent2):
    chain_names_one = [ch.GetName() for ch in ent1.chains]
    chain_names_two = [ch.GetName() for ch in ent2.chains]
    if not sorted(chain_names_one) == sorted(chain_names_two):
        return False
    chain_names = chain_names_one
    for chain_name in chain_names:
        ch1 = ent1.FindChain(chain_name)
        ch2 = ent2.FindChain(chain_name)
        if not compare_chains(ch1, ch2):
            return False
    if not compare_bonds(ent1, ent2):
        return False
    return True

class TestOMF(unittest.TestCase):
    def test_AU(self):
        ent, seqres, info = io.LoadMMCIF("testfiles/mmcif/3T6C.cif.gz", 
                                         seqres=True,
                                         info=True)
        omf = io.OMF.FromMMCIF(ent, info)
        omf_bytes = omf.ToBytes()
        loaded_omf = io.OMF.FromBytes(omf_bytes)
        loaded_ent = loaded_omf.GetAU()
        self.assertTrue(compare_ent(ent, loaded_ent))

if __name__== '__main__':
  from ost import testutils
  testutils.RunTests()
