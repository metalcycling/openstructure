:mod:`~ost.bindings.tmtools` - Structural superposition
================================================================================

.. module:: ost.bindings.tmtools
  :synopsis: Sequence dependent and independent structure superposition

The :mod:`~ost.bindings.tmtools` module provides access to the structural 
superposition programs TMscore and Tmalign developed by Y. Zhang 
and J. Skolnick. These programs superpose a model onto a reference structure, 
using the positions of the Calpha atoms only. While at their core, these 
programs essentially use the same algorithm, they differ on how the Calphas are 
paired. TMscore pairs the Calpha atom based on the residue number, TMalign 
calculates an optimal pairing of Calpha atom based on heuristics.

Citation:

  Yang Zhang and Jeffrey Skolnick, Proteins 2004 57: 702-710
  Y. Zhang and J. Skolnick, Nucl. Acids Res. 2005 33, 2302-9

Besides using the standalone TM-align program, ost also provides wrappers 
around USalign as published in:

  Chengxin Zhang, Morgan Shine, Anna Marie Pyle, Yang Zhang
  (2022) Nat Methods

USalign can be used as external standalone tool. Alternatively, ost allows to
directly inject structural data in the USAlign c++ layer. The advantage of that
is that no intermediate files must be generated. 


Distance measures used by TMscore
--------------------------------------------------------------------------------

There are many different ways to describe the structural similarity of two 
protein structures at the Calpha level. TMscore calculate several of these 
measures. The most common is to describe the difference in terms of the root 
mean square deviation of the Calpha positions, the RMSD. Despite its common use, 
RMSD has several drawbacks when working with incomplete models. Since the RMSD 
highly depends on the set of included atoms, it is relatively easy to obtain a 
smaller RMSD by omitting flexible parts of a protein structure. This has lead to 
the introduction of the global distance test (GDT). A model is compared to a 
reference by calculating the fraction of Calpha atoms that can be superposed 
below a certain cutoff, e.g. 1Å. The fractions of several such cutoffs are 
combined into the GDT_TS (1, 2, 4 and 8Å) and GDT_HA (0.5, 1, 2, 4Å) and divided 
by four to obtain the final measure. In contrast to RSMD, GDT is an agreement 
measure. The higher the value, the more similar the two structures are. TM-score 
(not to be confused by TMscore, the program), additionally adds a size 
dependences to the GDT measure by taking the protein length into account. As 
with GDT, the bigger the value, the more similar the two structures are.

Common Usage
--------------------------------------------------------------------------------

The following example shows how to use TMscore to superpose two protein 
structures and print the RMSD as well as the GDT_TS and GDT_HA similarity measures.

.. code-block:: python

  from ost.bindings import tmtools
  
  pdb1=io.LoadPDB('1ake.pdb', restrict_chains='A')
  pdb2=io.LoadPDB('4ake.pdb', restrict_chains='A')
  result=tmtools.TMScore(pdb1, pdb2)
  print(result.rmsd_below_five) # 1.9
  print(result.gdt_ha) # 0.41
  print(result.gdt_ts) # 0.56

Usage of TMalign
--------------------------------------------------------------------------------

.. autofunction:: ost.bindings.tmtools.TMAlign


Usage of TMscore
--------------------------------------------------------------------------------

.. autofunction:: ost.bindings.tmtools.TMScore

.. autoclass:: ost.bindings.tmtools.TMScoreResult

Usage of USalign
--------------------------------------------------------------------------------

For higher order complexes, ost provides access to USalign. This corresponds to
calling USalign with the preferred way of comparing full biounits:

.. code-block:: bash

  USalign mdl.pdb ref.pdb -mm 1 -ter 0

.. autofunction:: ost.bindings.tmtools.USAlign


C++ wrappers
--------------------------------------------------------------------------------

.. currentmodule:: ost.bindings

Instead of calling the external executables, ost also provides a wrapper around
the USalign c++ implementation which is shipped with the ost source code.
The advantage is that no intermediate files need to be  generated.

.. code-block:: python

  from ost import bindings
  
  pdb1=io.LoadPDB('1ake.pdb').Select("peptide=true")
  pdb2=io.LoadPDB('4ake.pdb').Select("peptide=true")
  result = bindings.WrappedTMAlign(pdb1.chains[0], pdb2.chains[0], 
                                   fast=True)
  print(result.tm_score)
  print(result.alignment.ToString(80))


.. class:: TMAlignResult(rmsd, tm_score, aligned_length, transform, alignment)

  All parameters of the constructor are available as attributes of the class

  :param rmsd:          RMSD of the superposed residues
  :param tm_score:      TMScore of the superposed residues
  :param tm_score_swapped: TMScore when reference is swapped
  :param aligned_length: Number of superposed residues
  :param transform:     Transformation matrix to superpose first chain onto 
                        reference
  :param alignment:     The sequence alignment given the structural superposition
  :type rmsd:           :class:`float`
  :type tm_score:       :class:`float`
  :type aligned_length: :class:`int`
  :type transform:      :class:`geom.Mat4`
  :type alignment:      :class:`ost.seq.AlignmentHandle`

.. method:: WrappedTMAlign(chain1, chain2, [fast=False])

  Takes two chain views and runs TMalign from USAlign with *chain2* as
  reference. The positions and sequences are directly extracted from the chain
  residues for every residue that fulfills:
  
    * peptide linking and valid CA atom OR nucleotide linking and valid C3'
      atom
    * valid one letter code(no '?')

  The function automatically identifies whether the chains consist of peptide
  or RNA residues. An error is raised if the two types are mixed.

  :param chain1:        Chain from which position and sequence are extracted
                        to run TMalign.
  :param chain2:        Chain from which position and sequence are extracted
                        to run TMalign, this is the reference.
  :param fast:          Whether to apply the *fast* flag to TMAlign
  :type chain1:         :class:`ost.mol.ChainView`
  :type chain2:         :class:`ost.mol.ChainView`
  :type fast:           :class:`bool`
  :rtype:               :class:`ost.bindings.TMAlignResult`


.. method:: WrappedTMAlign(pos1, pos2, seq1, seq2 [fast=False, rna=False])

  Similar as described above, but directly feeding in raw data.

  :param pos1:          CA/C3' positions of the first chain
  :param pos2:          CA/C3' positions of the second chain, this is the reference.
  :param seq1:          Sequence of first chain
  :param seq2:          Sequence of second chain
  :param fast:          Whether to apply the *fast* flag to TMAlign
  :param rna:           Whether to treat as RNA
  :type pos1:           :class:`ost.geom.Vec3List`
  :type pos2:           :class:`ost.geom.Vec3List`
  :type seq1:           :class:`ost.seq.SequenceHandle`
  :type seq2:           :class:`ost.seq.SequenceHandle`
  :type fast:           :class:`bool`
  :type rna:            :class:`bool`
  :rtype:               :class:`ost.bindings.TMAlignResult`
  :raises:              :class:`ost.Error` if pos1 and seq1, pos2 and seq2 
                        respectively are not consistent in size.

For higher order complexes, ost provides access to the MMalign functionality
from USalign. 

.. class:: MMAlignResult(rmsd, tm_score, transform, aligned_length, alignments,\
                         ent1_mapped_chains, ent2_mapped_chains)

  All parameters of the constructor are available as attributes of the class

  :param rmsd:          RMSD of the superposed residues
  :param tm_score:      TMScore of the superposed residues
  :param tm_score_swapped: TMScore when reference is swapped
  :param aligned_length: Number of superposed residues
  :param transform:     Transformation matrix to superpose mdl onto reference 
  :param alignments:    Alignments of all mapped chains, with first sequence
                        being from ent1 and second sequence from ent2
  :param ent1_mapped_chains: All mapped chains from ent1
  :param ent2_mapped_chains: The respective mapped chains from ent2
  :type rmsd:           :class:`float`
  :type tm_score:       :class:`float`
  :type aligned_length: :class:`int`
  :type transform:      :class:`geom.Mat4`
  :type alignments:     :class:`ost.seq.AlignmentList`
  :type ent1_mapped_chains: :class:`ost.StringList` 
  :type ent2_mapped_chains: :class:`ost.StringList`

.. method:: WrappedMMAlign(ent1, ent2, [fast=False])

  Takes two entity views and runs MMalign with *ent2* as reference.
  The positions and sequences are directly extracted from the chain
  residues for every residue that fulfills:
  
    * peptide linking and valid CA atom OR nucleotide linking and valid C3'
      atom
    * valid one letter code(no '?')

  The function automatically identifies whether the chains consist of peptide
  or RNA residues. An error is raised if the two types are mixed in the same
  chain.

  :param ent1:          Entity from which position and sequence are extracted
                        to run MMalign.
  :param ent2:          Entity from which position and sequence are extracted
                        to run MMalign, this is the reference.
  :param fast:          Whether to apply the *fast* flag to MMAlign
  :type ent1:           :class:`ost.mol.EntityView`
  :type ent2:           :class:`ost.mol.EntityView`
  :type fast:           :class:`bool`
  :rtype:               :class:`ost.bindings.MMAlignResult`
