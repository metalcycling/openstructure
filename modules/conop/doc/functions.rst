Conop Functions
=======================================================================

.. currentmodule:: ost.conop

.. function:: AssignBackboneTorsions(prev, res, next)
              AssignBackboneTorsions(chain)
              AssignBackboneTorsions(residues)
   
  Assigns the backbone torsions PHI, PSI and OMEGA. The backbone atoms are
  required to be connected for the torsions to be added. In addition, only
  residues for which :attr:`~ost.mol.ResidueHandle.peptide_linking` is True are
  considered.

  The first signature assigns the torsions to *res*, assuming prev is 
  the amino acid before, and *next* is the amino acid next to *res*.  
  Both *next* and *prev* may be invalid residues.

  The second and third signature assign the torsions to all peptidic 
  residues of the chain/the residue list. The residues in the chain, 
  the residue list are thought to run from N to C terminus.

  :param prev: The amino acid before *res*
  :type prev:  :class:`~ost.mol.ResidueHandle`
  :param res: The central amino acid. Must be a valid handle.
  :type res:  :class:`~ost.mol.ResidueHandle`
  :param next: The amino acid after *res*
  :type next:  :class:`~ost.mol.ResidueHandle`
  :param chain: Chain to process.
  :type chain:  :class:`~ost.mol.ChainHandle`
  :param residues: List of residues to process.
  :type residues:  :class:`~ost.mol.ResidueHandleList`

.. function:: IsBondFeasible(atom_a, atom_b)

  :return: True, if *atom_a* and *atom_b* are within a reasonable distance for a
           bond to be present. Depends on :attr:`~ost.mol.AtomHandle.radius` of
           atoms and heuristic formulas.
  :rtype:  :class:`bool`
  :param atom_a: Atom to be checked.
  :type atom_a:  :class:`~ost.mol.AtomHandle`
  :param atom_b: Atom to be checked.
  :type atom_b:  :class:`~ost.mol.AtomHandle`

.. function:: GetUnknownAtomsOfResidue(residue, compound, strict_hydrogens=False)

  Returns the list of atoms present in *residue* that are not part of the 
  atom specifications in *compound*.

  :param residue: The residue to analyze.
  :type residue:  :class:`~ost.mol.ResidueHandle`
  :param compound: Chemical compound to compare with.
  :type compound:  :class:`Compound`
  :param strict_hydrogens: When set to True, hydrogen atoms are checked as well.
  :type strict_hydrogens:  :class:`bool`