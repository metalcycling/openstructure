from ost import io, mol, geom, seq
import unittest
import os
import random

class TestPDBize(unittest.TestCase):
  def test_numbers_water_molecules_with_ins_codes(self):
    m = mol.CreateEntity()
    e = m.EditXCS(mol.BUFFERED_EDIT)
    c = e.InsertChain("A");
    e.SetChainType(c, mol.CHAINTYPE_WATER)
    for i in range(27):
      e.AppendResidue(c, "HOH")
    pdbizer = mol.alg.PDBize()
    transformations = geom.Mat4List()
    transformations.append(geom.Mat4())
    seqs = seq.CreateSequenceList()
    pdbizer.Add(m.Select(''), transformations, seqs)
    pdbized = pdbizer.Finish()
    self.assertEqual([c.name for c in pdbized.chains], ["-"])
    residues = pdbized.residues
    for i in range(26):
      self.assertEqual(residues[i].number.num, 1)
      self.assertEqual(residues[i].number.ins_code, chr(ord('A')+i))
    self.assertEqual(residues[26].number.num, 2)
    self.assertEqual(residues[26].number.ins_code, 'A')

  def test_starts_from_last_water_rnum(self):
    m = mol.CreateEntity()
    e = m.EditXCS(mol.BUFFERED_EDIT)
    c = e.InsertChain("A");
    e.SetChainType(c, mol.CHAINTYPE_WATER)
    e.AppendResidue(c, "HOH")
    pdbizer = mol.alg.PDBize()
    transformations = geom.Mat4List()
    transformations.append(geom.Mat4())
    seqs = seq.CreateSequenceList()
    pdbizer.Add(m.Select(''), transformations,seqs)
    pdbizer.Add(m.Select(''), transformations,seqs)
    pdbized = pdbizer.Finish()
    self.assertEqual([c.name for c in pdbized.chains], ["-"])
    residues = pdbized.residues
    self.assertEqual([r.number for r in residues],
                      [mol.ResNum(1, 'A'), mol.ResNum(1, 'B')])

  def test_large_poly_ligand_inscode(self):
    """Make sure polymeric ligands with 26+ residues get split up and do not
    run into an overflow with the insertion codes.
    """
    m = mol.CreateEntity()
    e = m.EditXCS(mol.BUFFERED_EDIT)
    c = e.InsertChain("A");
    e.SetChainType(c, mol.CHAINTYPE_OLIGOSACCHARIDE)
    for i in range(27):
      e.AppendResidue(c, "NAG")
    pdbizer = mol.alg.PDBize(min_polymer_size=28)
    transformations = geom.Mat4List()
    transformations.append(geom.Mat4())
    seqs = seq.CreateSequenceList()
    pdbizer.Add(m.Select(''), transformations, seqs)
    pdbized = pdbizer.Finish()
    self.assertEqual([c.name for c in pdbized.chains], ["_"])
    residues = pdbized.residues
    for i in range(26):
      self.assertEqual(residues[i].number.num, 1)
      self.assertEqual(residues[i].number.ins_code, chr(ord('A')+i))
    self.assertEqual(residues[26].number.num, 2)
    self.assertEqual(residues[26].number.ins_code, 'A')

  def test_large_poly_ligand_resnum(self):
    """Make sure polymeric ligands with 26 residues do not increase the residue
    number by 2. Otherwise the next ligand is off by 1, like NAG26Z, ZN28.
    """
    m = mol.CreateEntity()
    e = m.EditXCS(mol.BUFFERED_EDIT)
    c = e.InsertChain("A");
    e.SetChainType(c, mol.CHAINTYPE_OLIGOSACCHARIDE)
    for i in range(26):
      e.AppendResidue(c, "NAG")
    c = e.InsertChain("B");
    e.AppendResidue(c, "ZN")
    pdbizer = mol.alg.PDBize(min_polymer_size=28)
    transformations = geom.Mat4List()
    transformations.append(geom.Mat4())
    seqs = seq.CreateSequenceList()
    pdbizer.Add(m.Select(''), transformations, seqs)
    pdbized = pdbizer.Finish()
    self.assertEqual([c.name for c in pdbized.chains], ["_"])
    residues = pdbized.residues
    for i in range(26):
      self.assertEqual(residues[i].number.num, 1)
      self.assertEqual(residues[i].number.ins_code, chr(ord('A')+i))
    self.assertEqual(residues[26].number.num, 2)
    self.assertEqual(residues[26].number.ins_code, '\0')

  def _CheckMinSize(self, ost_ent, seq_list, chn_nm_lst, **kwargs):
    """Check effects of the *_min_size parameter.

    :param ost_ent: OST entity to be PDBized.
    :type ost_ent: :class:`~ost.mol.EntityHandle`
    :param seq_list: Sequence list for the chains in ost_ent.
    :type seq_list: :class:`~ost.seq.SequenceList`
    :param chn_nm_lst: List of expected chain names in PDBized entity.
    :type chn_nm_lst: :class:`list` of :class:`str`
    """
    if "saccharide_min_size" not in kwargs:
      kwargs["saccharide_min_size"] = 10
    if "nucleicacid_min_size" not in kwargs:
      kwargs["nucleicacid_min_size"] = 10
    if "peptide_min_size" not in kwargs:
      kwargs["peptide_min_size"] = 10
    transformations = geom.Mat4List()
    transformations.append(geom.Mat4())
    pdbizer = mol.alg.PDBize(**kwargs)
    pdbizer.Add(ost_ent.Select(''), transformations, seq_list)
    pdbized = pdbizer.Finish()
    self.assertEqual(len(pdbized.chains), len(chn_nm_lst))
    for i in range(0, len(chn_nm_lst)):
      self.assertEqual(pdbized.chains[i].name, chn_nm_lst[i])
    return pdbized

  def test_peptide_min_size(self):
    """Make sure the peptide_min_size parameter works, place a polypeptide in
    chain '_'.
    """
    m = mol.CreateEntity()
    e = m.EditXCS(mol.BUFFERED_EDIT)
    c = e.InsertChain("A");
    e.SetChainType(c, mol.CHAINTYPE_POLY_PEPTIDE_L)
    for i in range(10):
      e.AppendResidue(c, "ALA")
    seqs = seq.CreateSequenceList()
    seqs.AddSequence(seq.CreateSequence("LotsOfAlanin", "AAAAAAAAAA"))

    # test that small peptide chains end up in the ligand chain "_"
    self._CheckMinSize(m, seqs, ["_"], peptide_min_size=11)

    # test again with two small peptide chains
    c = e.InsertChain("B");
    e.SetChainType(c, mol.CHAINTYPE_POLY_PEPTIDE_L)
    for i in range(15):
      e.AppendResidue(c, "ALA")
    seqs.AddSequence(seq.CreateSequence("MoreAlanin", "AAAAAAAAAAAAAAA"))
    self._CheckMinSize(m, seqs, ["_"], peptide_min_size=16)

    # test one peptide in ligand chain, second as polymer chain
    self._CheckMinSize(m, seqs, ["_", "A"], peptide_min_size=11)

    # actually disabling min. polymer size
    self._CheckMinSize(m, seqs, ["A", "B"], peptide_min_size=0)

  def test_nucleicacid_min_size(self):
    """Make sure the nucleicacid_min_size parameter works, place a
    polynucleotide in chain '_'.
    """
    m = mol.CreateEntity()
    e = m.EditXCS(mol.BUFFERED_EDIT)
    c = e.InsertChain("A");
    e.SetChainType(c, mol.CHAINTYPE_POLY_DN)
    for i in range(10):
      e.AppendResidue(c, "DA")
    seqs = seq.CreateSequenceList()
    seqs.AddSequence(seq.CreateSequence("LotsOfAdenine", "AAAAAAAAAA"))

    # test that small nucleotide chains end up in the ligand chain "_"
    self._CheckMinSize(m, seqs, ["_"], nucleicacid_min_size=11)

    # test again with two small nucleic acid chains
    c = e.InsertChain("B");
    e.SetChainType(c, mol.CHAINTYPE_POLY_DN)
    for i in range(15):
      e.AppendResidue(c, "DA")
    seqs.AddSequence(seq.CreateSequence("MoreAdenine", "AAAAAAAAAAAAAAA"))
    self._CheckMinSize(m, seqs, ["_"], nucleicacid_min_size=16)

    # test one nucleic acid in ligand chain, second as polymer chain
    self._CheckMinSize(m, seqs, ["_", "A"], nucleicacid_min_size=11)

    # actually disabling min. polymer size
    self._CheckMinSize(m, seqs, ["A", "B"], nucleicacid_min_size=0)

  def test_saccharide_min_size(self):
    """Make sure the saccharide_min_size parameter works, place an
    oligosaccharide in chain '_'.
    """
    m = mol.CreateEntity()
    e = m.EditXCS(mol.BUFFERED_EDIT)
    c = e.InsertChain("A");
    e.SetChainType(c, mol.CHAINTYPE_OLIGOSACCHARIDE)
    for i in range(10):
      e.AppendResidue(c, "NAG")
    seqs = seq.CreateSequenceList()

    # test that small oligosaccharides end up in the ligand chain "_"
    self._CheckMinSize(m, seqs, ["_"], saccharide_min_size=11)

    # test again with two small oligosaccharide chains
    c = e.InsertChain("B");
    e.SetChainType(c, mol.CHAINTYPE_OLIGOSACCHARIDE)
    for i in range(15):
      e.AppendResidue(c, "NAG")
    self._CheckMinSize(m, seqs, ["_"], saccharide_min_size=16)

    # test one oligosaccharide in ligand chain, second as polymer chain
    self._CheckMinSize(m, seqs, ["_", "A"], saccharide_min_size=11)

    # actually disabling min. polymer size

  def test_peptide_nucleicacid_saccharide_min_sizes(self):
    """Make sure that all thre thresholds play well together.
    """
    m = mol.CreateEntity()
    e = m.EditXCS(mol.BUFFERED_EDIT)
    c = e.InsertChain("A");
    e.SetChainType(c, mol.CHAINTYPE_POLY_PEPTIDE_L)
    for i in range(10):
      e.AppendResidue(c, "ALA")
    seqs = seq.CreateSequenceList()
    seqs.AddSequence(seq.CreateSequence("LotsOfAlanin", "AAAAAAAAAA"))
    c = e.InsertChain("B");
    e.SetChainType(c, mol.CHAINTYPE_POLY_DN)
    for i in range(10):
      e.AppendResidue(c, "DA")
    seqs.AddSequence(seq.CreateSequence("LotsOfAdenine", "AAAAAAAAAA"))
    c = e.InsertChain("C");
    e.SetChainType(c, mol.CHAINTYPE_OLIGOSACCHARIDE)
    for i in range(10):
      e.AppendResidue(c, "NAG")

    # Check branched entities can be abandoned in the ligand chain while
    # peptides and nucleic acids live in their own chains.
    pdbized = self._CheckMinSize(m, seqs, ["A", "B", "_"],
                                 saccharide_min_size=11,
                                 peptide_min_size=0,
                                 nucleicacid_min_size=0)
    self.assertTrue(pdbized.chains[0].IsPolypeptide())
    self.assertTrue(pdbized.chains[1].IsPolynucleotide())
    self.assertEqual(pdbized.chains[2].residues[0].GetStringProp("type"),
                     "oligosaccharide")

    # test to store a short polynucleotide and sugar in the ligand chain but keep
    # longer polynucleotide and the peptide outside of the ligand chain.
    c = e.InsertChain("D");
    e.SetChainType(c, mol.CHAINTYPE_POLY_DN)
    for i in range(5):
      e.AppendResidue(c, "DG")
    seqs.AddSequence(seq.CreateSequence("LotsOfGuanine", "GGGGG"))
    pdbized = self._CheckMinSize(m, seqs, ["A", "B", "_"],
                                 saccharide_min_size=11,
                                 peptide_min_size=0,
                                 nucleicacid_min_size=6)
    self.assertTrue(pdbized.chains[0].IsPolypeptide())
    self.assertTrue(pdbized.chains[1].IsPolynucleotide())
    self.assertEqual(pdbized.chains[2].residues[0].GetStringProp("type"),
                     "oligosaccharide")
    self.assertEqual(pdbized.chains[2].residues[-1].GetStringProp("type"),
                     "polydeoxyribonucleotide")

    # test to add a small peptide to the ligand chain
    c = e.InsertChain("E");
    e.SetChainType(c, mol.CHAINTYPE_POLY_PEPTIDE_L)
    for i in range(5):
      e.AppendResidue(c, "ALA")
    seqs.AddSequence(seq.CreateSequence("SomeAlanin", "AAAAA"))
    pdbized = self._CheckMinSize(m, seqs, ["A", "B", "_", "C"],
                                 saccharide_min_size=11,
                                 peptide_min_size=6,
                                 nucleicacid_min_size=3)
    self.assertTrue(pdbized.chains[0].IsPolypeptide())
    self.assertTrue(pdbized.chains[1].IsPolynucleotide())
    self.assertEqual(pdbized.chains[2].residues[0].GetStringProp("type"),
                     "oligosaccharide")
    self.assertEqual(pdbized.chains[2].residues[-1].GetStringProp("type"),
                     "polypeptide(L)")


if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()
