:mod:`~ost.bindings.dssp` - Secondary structure assignment
================================================================================

.. module:: ost.bindings.dssp
  :synopsis: Interface to the DSSP command line utility

Introduction
--------------------------------------------------------------------------------

DSSP is a program developed by Wolfgang Kabsch and Chris Sander to assign
secondary structure states to protein structures. The assignment is based on
hydrogen bonding patterns and geometric features.

The program can be downloaded from `<http://swift.cmbi.ru.nl/gv/dssp/>`_.

.. warning::
  Support of the DSSP program has been deprecated. OpenStructure provides
  functionality to assign equivalent secondary structures
  (:func:`ost.mol.alg.AssignSecStruct`) and solvent accessibility
  (:func:`ost.mol.alg.Accessibility`). You're advised to use these
  algorithms.

  :func:`AssignDSSP` still exists and provides the "old" interface but
  internally uses the OpenStructure impmlementations and does not call
  an external program anymore.


Examples
--------------------------------------------------------------------------------

The following example assigns secondary structure states to an entity by using
the DSSP program.


.. code-block:: python

  from ost.bindings import dssp
  ent=io.LoadPDB('1ake.pdb')
  dssp.AssignDSSP(ent)


Now we fetch structure information plus solvent accessibility for an entity
using the mmCIF interface.


.. code-block:: python

  from ost.bindings import dssp
  ent=io.LoadMMCIF('1ake.cif')
  dssp.AssignDSSP(ent, extract_burial_status=True)
  for chain in ent.chains:
    if chain.is_polypeptide:
      for res in chain.residues:
        print(res.GetFloatProp('relative_solvent_accessibility'))


DSSP bindings Usage
--------------------------------------------------------------------------------

.. autofunction:: ost.bindings.dssp.AssignDSSP

.. autofunction:: ost.bindings.dssp.LoadDSSP

..  LocalWords:  dssp AssignDSSP ent GetFloatProp autofunction
