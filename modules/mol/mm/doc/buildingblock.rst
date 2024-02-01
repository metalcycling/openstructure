Blocks
================================================================================

.. currentmodule:: ost.mol.mm

The most basic type of residue description is the BuildingBlock. It contains
information on atom names and their corresponding types, charges and
optionally also their masses. Interactions for all particles can also
be defined even though they are optional, except the bonds giving
information about the connectivity. 
You often need to manipulate building blocks or
the residues they refer to in an automated manner. That's where the 
BlockModifiers come in, with the GromacsBlockModifier as a specific
implementation. As a special case there also exist HydrogenConstructors.


The BuildingBlock Class
--------------------------------------------------------------------------------

.. class:: BuildingBlock

  .. method:: Match(residue, [,match_connectivity=True])

    Checks, whether the given residue matches the atom names in the 
    BuildingBlock. The connectivity gets checked optionally.

    :param residue:    
    :param match_connectivity: If set to true, the function checks the bonds
                               in the residue with the defined bonds in the
                               BuildingBlock 

    :type residue:    :class:`ResidueHandle`
    :type match_connectivity:    :class:`bool`

    :returns: :class:`bool`

  .. method:: Connect(residue, editor)

    Connects atoms of **residue** based on the bond definitions of the 
    BuildingBlock.

    :param residue:     Residue to be connected
    :param editor:      Editor associated to the residue's entity

    :type residue:      :class:`ResidueHandle`
    :type editor:       :class:`XCSEditor`

    :raises:            :class:`RuntimeError` when required atom can not be found in 
                        residue


  .. method:: AddAtom(name, type, charge,[,mass = None])

    :param name:        Name of atom
    :param type:        Its corresponding forcefield type
    :param charge:      Its charge
    :param mass:        Its mass

    :type name:         :class:`str`
    :type type:         :class:`str`
    :type charge:       :class:`float`
    :type mass:         :class:`float`



  .. method:: RemoveAtom(name)

    Removes atom from BuildingBlock with all its associated values and 
    interactions

    :param name:        Name of atom to be removed

    :type name:         :class:`str`



  .. method:: ReplaceAtom(name, new_name, new_type, new_charge, [,new_mass = None])

    Replace given atom by resetting name, type, charge and mass. All interactions
    containing that atom get adapted as well

    :param name:        Name of atom to be replaced
    :param new_name:    New name of atom
    :param new_type:    New type of atom
    :param new_charge:  New charge of atom
    :param new_mass:    New mass of atom

    :type name:         :class:`str`
    :type new_name:     :class:`str`
    :type new_type:     :class:`str`
    :type new_charge:   :class:`float`
    :type new_mass:     :class:`float`



  .. method:: RemoveInteractionsToPrev()

    Removes all interactions associated to an atom of a previous residues.
    This gets indicated in the :class:`BuildingBlock` by an atom name prefixed by a *-*


  .. method:: RemoveInteractionsToNext()

    Removes all interactions associated to an atom of a next residues.
    This gets indicated in the :class:`BuildingBlock` by an atom name prefixed by a *+*


  .. method:: AddBond(bond, [,replace_existing = False])

    :param bond:        Bond to be added
    :param replace_existing: Whether potentially already existing bond for the
                             same atoms should be replaced.

    :type bond:         :class:`Interaction`
    :type replace_existing: :class:`bool`


  .. method:: AddAngle(angle, [,replace_existing = False])

    :param angle:       Angle to be added
    :param replace_existing: Whether a potentially already existing angle for the
                             same atoms should be replaced.

    :type angle:         :class:`Interaction`
    :type replace_existing: :class:`bool`


  .. method:: AddDihedral(dihedral, [,replace_existing = False])

    :param dihedral:    Dihedral to be added
    :param replace_existing: Whether potentially already existing dihedral for the
                             same atoms should be replaced.

    :type dihedral:     :class:`Interaction`
    :type replace_existing: :class:`bool`


  .. method:: AddImproper(improper, [,replace_existing = False])

    :param improper:    Improper to be added
    :param replace_existing: Whether potentially already existing improper for the
                             same atoms should be replaced.

    :type improper:     :class:`Interaction`
    :type replace_existing: :class:`bool`


  .. method:: AddExclusion(exclusion, [,replace_existing = False])

    :param exclusion:        Exclusion to be added
    :param replace_existing: Whether potentially already existing Exclusion for the
                             same atoms should be replaced.

    :type exclusion:         :class:`Interaction`
    :type replace_existing: :class:`bool`


  .. method:: AddCMap(cmap, [,replace_existing = False])

    :param cmap:        CMap to be added
    :param replace_existing: Whether potentially already existing cmap for the
                             same atoms should be replaced.

    :type cmap:         :class:`Interaction`
    :type replace_existing: :class:`bool`


  .. method:: AddConstraint(constraint, [,replace_existing = False])

    :param constraint:  Constraint to be added
    :param replace_existing: Whether potentially already existing constraint for the
                             same atoms should be replaced.

    :type constraint:         :class:`Interaction`
    :type replace_existing: :class:`bool` 


  .. method:: GetType(name)

    Gets forcefield type from atom with given name

    :param name:        Name of atom you want the type from

    :type name:         :class:`str`

    :returns: :class:`str`

    :raises:            :class:`RuntimeError` when atom can not be found in 
                        BuildingBlock    

  .. method:: GetCharge(name)

    Gets charge from atom with given name

    :param name:        Name of atom you want the charge from

    :type name:         :class:`str`

    :returns: :class:`float`

    :raises: :class:`RuntimeError` when atom can not be found in 
                                   BuildingBlock 


  .. method:: GetMass(name)

    Gets mass from atom with given name

    :param name:        Name of atom you want the mass from

    :type name:         :class:`str`

    :returns: :class:`float`
    
    :raises: :class:`RuntimeError` when atom can not be found in 
                                   BuildingBlock 

  .. method:: GetAtoms()

    :returns: :class:`list` of all atom names

  .. method:: GetTypes()

    :returns: :class:`list` of all atom types

  .. method:: GetCharges()

    :returns: :class:`list` of all charges

  .. method:: GetMasses()

    :returns: :class:`list` of all masses

  .. method:: GetBonds()

    :returns: :class:`list` of all bonds


  .. method:: GetAngles()

    :returns: :class:`list` of all angles


  .. method:: GetDihedrals()

    :returns: :class:`list` of all dihedrals


  .. method:: GetImpropers()

    :returns: :class:`list` of all impropers


  .. method:: GetCMaps()

    :returns: :class:`list` of all cmaps


  .. method:: GetExclusions()

    :returns: :class:`list` of all exlusions


  .. method:: GetConstraints()

    :returns: :class:`list` of all constraints

Block Modifiers
--------------------------------------------------------------------------------
.. class:: BlockModifier
  
  Basis class. Block modifiers are used to change building blocks or residues.
  See :class:`GromacsBlockModifier` for a specific example.

.. class:: GromacsBlockModifier

  .. method:: ApplyOnBuildingBlock(block)

    Applies all defined rules on the given :class:`BuildingBlock`

    :param block:       BuildingBlock to be modified
    :type block:        :class:`BuildingBlock`

  .. method:: ApplyOnResidue(residue)

    Applies all defined rules on the given :class:`ResidueHandle`

    :param residue:     Residue to be modified
    :type residue:      :class:`ResidueHandle`

  .. method:: AddReplaceRule(name, new_name, new_type, new_charge)

    Rule, that basically renames an atom and also resets its type and charge
    in a :class:`BuildingBlock`. A simple renaming occurs in a :class:`ResidueHandle`.

    :param name:        Name of the atom to be changed
    :param new_name:    Its new name
    :param new_type:    Its new type
    :param new_charge:  Its new charge

    :type name:         :class:`str`
    :type new_name:     :class:`str`
    :type new_type:     :class:`str`
    :type new_charge:   :class:`float`

  .. method:: AddAddRule(number, method, atom_names, anchors, type, charge)

    A rule to add new atoms the Gromacs way, see
    `Gromacs Manual <http://www.gromacs.org/Documentation/Manual>`_ 
    for the exact definition of the parameters. A :class:`BuildingBlock`
    gets modified by adding the new atom definitions and also the
    corresponding bonds describing the connectivity. In case of 
    :class:`ResidueHandle` the new Atoms with connectivity get added with 
    their positions as defined by the Gromacs adding rule. 

    :param number:      Number of atoms to be added
    :param method:      Gromacs adding rule
    :param atom_names:  Strings containing the new atom names
    :param anchors:     Strings containing atom names used as anchor
    :param type:        The type the atoms will have
    :param charge:      The charge the atoms will have

    :type number:      :class:`int`
    :type method:      :class:`int` 
    :type atom_names:  :class:`list`
    :type anchors:     :class:`list`
    :type type:        :class:`str`
    :type charge:      :class:`float`

  .. method:: AddDeleteAtom(name)

    Defines an atom that has to be removed. In case of the 
    :class:`BuildingBlock` this removes this particular atom plus all 
    interactions connected to it, in case of :class:`ResidueHandle`, the
    atom simply gets deleted.

    :param name:        Atom to be deleted

    :type name:         :class:`str`

  .. method:: AddBond(bond)

    Adds a bond, this only has effect on :class:`BuildingBlock`, not
    on :class:`ResidueHandle` when the corresponding Apply function gets 
    called

    :param bond:        Bond to be added
    :type bond:         :class:`Interaction`

  .. method:: AddAngle(angle)

    Adds an angle, this only has effect on :class:`BuildingBlock`, not
    on :class:`ResidueHandle` when the corresponding Apply function gets 
    called

    :param angle:        Angle to be added
    :type angle:        :class:`Interaction`

  .. method:: AddDihedral(dihedral)

    Adds a dihedral, this only has effect on :class:`BuildingBlock`, not
    on :class:`ResidueHandle` when the corresponding Apply function gets 
    called

    :param dihedral:    Dihedral to be added
    :type dihedral:     :class:`Interaction`

  .. method:: AddImproper(improper)

    Adds an improper, this only has effect on :class:`BuildingBlock`, not
    on :class:`ResidueHandle` when the corresponding Apply function gets 
    called

    :param improper:    Improper to be added
    :type improper:     :class:`Interaction`

  .. method:: AddCMap(cmap)

    Adds a cmap, this only has effect on :class:`BuildingBlock`, not
    on :class:`ResidueHandle` when the corresponding Apply function gets 
    called

    :param cmap:        CMap to be added
    :type cmap:        :class:`Interaction`


Hydrogen Constructors
--------------------------------------------------------------------------------

.. class:: HydrogenConstructor

  Basis class. Hydrogen constructors are used to add hydrogens to residues.


.. class:: GromacsHydrogenConstructor

  The :class:`GromacsHydrogenConstructor` is the Gromacs way of adding
  hydrogens to a structure.

  .. method:: ApplyOnBuildingBlock(block)

    Guess what it does: !!ABSOLUTELY NOTHING!! just there for consistency

    :param block:       Block that won't be changed at all!
                        Isn't that awesome?

    :type block:        :class:`BuildingBlock`

  .. method:: ApplyOnResidue(residue)

    Constructs hydrogens based on the defined hydrogen addition rules

    :param residue:     Residue to be modified

    :type residue:      :class:`ResidueHandle`

  .. method:: AddHydrogenRule(number, method, hydrogen_names, anchors)

    Adds a hydrogen building rule as defined in 
    `Gromacs Manual <http://www.gromacs.org/Documentation/Manual>`_ 
    
    :param number:      Number of hydrogens to be added
    :param method:      Gromacs adding rule
    :param hydrogen_names:  Strings containing the hydrogen names
    :param anchors:     Strings containing atom names used as anchor

    :type number:       :class:`int`
    :type method:       :class:`int`
    :type hydrogen_names:  :class:`list`
    :type anchors:      :class:`list`

.. class:: HeuristicHydrogenConstructor(block)

  As soon as we leave the well defined world of Gromacs residue definitions,
  we have to find new ways for constructing hydrogens. The 
  :class:`HeuristicHydrogenConstructor` takes a :class:`BuildingBlock` at
  initialization and builds heuristic rules to build hydrogens based on 
  the connecticity defined in the given block.

  :param block:         :class:`BuildingBlock` from which the connectivity
                        information for hydrogen construction is extracted

  .. method:: ApplyOnBuildingBlock(block)

    Guess what it does: !!ABSOLUTELY NOTHING!! Just there for consistency

    :param block:       Block that won't be changed at all!
                        Isn't that awesome?

    :type block:        :class:`BlockModifier`

  .. method:: ApplyOnResidue(residue)

    Constructs hydrogen based on heuristic rules

    :param residue:     Residue to be modified

    :type residue:      :class:`ResidueHandle`

