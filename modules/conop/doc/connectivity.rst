Connectivity
================================================================================

.. currentmodule:: ost.conop


Motivation
--------------------------------------------------------------------------------


The connectivity of atoms is notoriously difficult to come by for biological 
macromolecules. PDB files, the de facto standard exchange format for structural 
information allows bonds to be specified in CONECT records. However, they are not
mandatory. Many programs, especially the ones not requiring on connectivity of 
atoms, do not write CONECT records. As a result, programs and structural biology 
frameworks can't rely on connectivity information to be present. The connectivity
information needs to be derived in the program itself.


Loader heuristics are great if you are the one that implemented them but are 
problematic if you are just the user of a software that has them. As time goes 
on, these heuristics become buried in thousands of lines of code and they are 
often hard yet impossible to trace back.

Different clients of the framework have different requirements. A visualisation 
software wants to read in a PDB files as is without making any changes. A 
script in an automated pipeline, however, does want to either strictly reject 
files that are incomplete or fill-in missing structural features. All these 
aspects are implemented in the conop module, separated from the loading of the 
PDB file, giving clients a fine grained control over the loading process. The
conop logic can thus be reused in code requiring the presence of 

The conop module defines a :class:`Processor` interface, to run connectivity 
algorithms, that is to connect the atoms with bonds and perform basic clean up 
of erroneous structures. The clients of the conop module can specify how the 
Processor should treat unknown amino acids, missing atoms and chemically 
infeasible bonds.

Processors
--------------------------------------------------------------------------------

The exact behaviour for a processor is implementation-specific. So far, two
classes implement the processor interface: A heuristic and a rule-based
processor. The processor mainly differ in the source of their connectivity
information. The Heuristicprocessor uses a hard-coded heuristic connectivity
table for the 20  standard amino acids as well as nucleotides.For other
compounds such as ligands the HeuristicProcessor runs a distance-based
connectivity algorithm that connects two atoms if they are closer than a certain
threshold. The RuleBasedProcessor uses a connectivity library containing all
molecular components present in the PDB files on PDB.org. The library can easily
be extended with custom  connectivity information, if required. If a
:doc:`compound library <compoundlib>` is present, the
:class:`RuleBasedProcessor` is enabled by default, otherwise the
:class:`HeuristicProcessor` is used as a fallback.


.. class:: Processor

  .. attribute:: check_bond_feasibility

    Whether an additional bond feasibility check is performed. Disabled by
    default. If turned on, atoms are only connected by bonds if they are within
    a reasonable distance (as defined by :func:`IsBondFeasible`).

    :type: :class:`bool`

  .. attribute:: assign_torsions

    Whether backbone torsions should be added to the backbone. Enabled by
    default. If turned on, PHI, PSI and OMEGA torsions are assigned to the
    peptide residues. See also :func:`AssignBackboneTorsions`.

    :type: :class:`bool`

  .. attribute:: connect

    Whether to connect atoms by bonds. Enabled by default. Turn this off if you
    would like to speed up the loading process and do not require connectivity
    information to be present in your structures. Note though that
    :attr:`peptide_bonds` may be ignored if this is turned off.

    :type: :class:`bool`

  .. attribute:: peptide_bonds

    Whether to connect residues by peptide bonds. Enabled by default. This also
    sets the :attr:`~ost.mol.ResidueHandle.is_protein` property of residues when
    peptide bonds are created. Turn this off if you would like to create your
    own peptide bonds.

    :type: :class:`bool`

  .. attribute:: zero_occ_treatment

    Controls the behaviour of importing atoms with zero occupancy. By default,
    this is set to warn.

    :type: :class:`ConopAction`

  .. method:: Process(ent)
  
    Processess the entity *ent* according to the current options.


.. class:: HeuristicProcessor(check_bond_feasibility=False, \
                              assign_torsions=True, connect=True, \
                              peptide_bonds=True, zero_occ_treatment=CONOP_WARN)
   
  The :class:`HeuristicProcessor` implements the :class:`Processor` interface.
  Refer to its documentation for methods and accessors common to all processor.

  :param check_bond_feasibility: Sets :attr:`~Processor.check_bond_feasibility`
  :param assign_torsions: Sets :attr:`~Processor.assign_torsions`
  :param connect: Sets :attr:`~Processor.connect`
  :param peptide_bonds: Sets :attr:`~Processor.peptide_bonds`
  :param zero_occ_treatment: Sets :attr:`~Processor.zero_occ_treatment`


.. class:: RuleBasedProcessor(compound_lib, fix_elements=True, \
                              strict_hydrogens=False, \
                              unknown_res_treatment=CONOP_WARN, \
                              unknown_atom_treatment=CONOP_WARN, \
                              check_bond_feasibility=False, \
                              assign_torsions=True, connect=True, \
                              peptide_bonds=True, zero_occ_treatment=CONOP_WARN)
   
  The :class:`RuleBasedProcessor` implements the :class:`Processor` interface.
  Refer to its documentation for methods and accessors common to all processor.

  :param compound_lib: The compound library to use
  :type compound_lib:  :class:`CompoundLib`
  :param fix_elements: Sets :attr:`fix_elements`
  :param strict_hydrogens: Sets :attr:`strict_hydrogens`
  :param unknown_res_treatment: Sets :attr:`unk_atom_treatment`
  :param unknown_atom_treatment: Sets :attr:`unk_res_treatment`
  :param check_bond_feasibility: Sets :attr:`~Processor.check_bond_feasibility`
  :param assign_torsions: Sets :attr:`~Processor.assign_torsions`
  :param connect: Sets :attr:`~Processor.connect`
  :param peptide_bonds: Sets :attr:`~Processor.peptide_bonds`
  :param zero_occ_treatment: Sets :attr:`~Processor.zero_occ_treatment`

  .. attribute:: fix_elements

    Whether the element of the atom should be changed to the atom defined in the
    compound library. Enabled by default.

    :type: :class:`bool`

  .. attribute:: strict_hydrogens

    Whether to use strict hydrogen naming rules outlined in the compound library.
    Disabled by default.

    :type: :class:`bool`

  .. attribute:: unk_atom_treatment

    Treatment upon encountering an unknown atom. Warn by default.

    :type: :class:`ConopAction`

  .. attribute:: unk_res_treatment

    Treatment upon encountering an unknown residue. Warn by default.

    :type: :class:`ConopAction`


.. class:: ConopAction

  Defines actions to take when certain events happen during processing. Possible
  values:

    ``CONOP_WARN``, ``CONOP_SILENT``, ``CONOP_REMOVE``, ``CONOP_REMOVE_ATOM``,
    ``CONOP_REMOVE_RESIDUE``, ``CONOP_FATAL``
