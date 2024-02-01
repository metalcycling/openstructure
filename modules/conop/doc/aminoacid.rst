Functions and classes for standard amino acids
================================================================================

.. currentmodule:: ost.conop

This document describes functions and classes to work with the 20 standard amino acids. The functions convert between different representations, e.g. one-letter-code, three-letter-code or the AminoAcid *enum*


The AminoAcid enum
--------------------------------------------------------------------------------

.. class:: AminoAcid

  The amino acid enum enumerates all 20 standard amino acid and the special
  value XXX, to signify an unknown amino acid. The amino acid enum supports the
  following values:

    ALA, ARG, ASN, ASP, GLN, GLU, LYS, SER, CYS, MET
    TRP, TYR, THR, VAL, ILE, LEU, GLY, PRO, HIS, PHE
    XXX

Converter functions
--------------------------------------------------------------------------------
.. function:: ResidueToAminoAcid(residue)
              ResidueNameToAminoAcid(rname) 
              OneLetterCodeToAminoAcid(olc) 

  :return: Amino acid from residue, residue name (three-letter-code)
           or one-letter-code. Returns XXX if residue, residue name or 
           one-letter-code is not one of the 20 standard amino acids.
  :rtype:  :class:`AminoAcid`

.. function:: OneLetterCodeToResidueName(olc)
              AminoAcidToResidueName(amino_acid)

  :return: Residue name from one-letter-code or amino_acid. For invalid 
           one-letter-codes or XXX, 'UNK' is returned.
  :rtype:  :class:`str`

  :param olc: One-letter-code
  :type olc:  :class:`str`
  :param amino_acid: Amino acid type
  :type amino_acid:  :class:`AminoAcid`


.. function:: ResidueNameToOneLetterCode(rname)

  :return: One-letter-code for the given residue name. Returns 'X' if residue
           name is not one of the 20 standard amino acids.
  :rtype:  :class:`str`

  :param rname: Residue name
  :type rname:  :class:`str`


.. class:: AminoAcidSet

  A set of amino acids, with constant-time access
  
  .. staticmethod:: CreatePolarSet()
                    CreateAromaticSet()
                    CreateApolarSet()

    Returns a set containing all polar, aromatic or apolar amino acids, 
    respectively.


  .. method:: Add(amino_acid)
  
    Add amino acid to the set.
    
  .. method:: Remove(amino_acid)
  
    Remove amino acid from the set
    
  .. method Contains(amino_acid)
  
    Whether the set contains the given amino acid. 
    
  .. method:: Empty()
  
    Whether the set is empty, i.e. doesn't contain any amino acids.