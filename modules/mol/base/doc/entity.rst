The Molecular Entity
================================================================================

.. currentmodule:: ost.mol

This document describes the :class:`EntityHandle` and the related classes.

The Handle Classes
--------------------------------------------------------------------------------

:class:`Entity <EntityHandle>`, :class:`residue <ResidueHandle>`,
:class:`atom <AtomHandle>` and :class:`bond <BondHandle>` handles can store
arbitrary :doc:`generic properties <../../base/generic>`.

.. function:: CreateEntity()

   Creates a new entity. The created entity is empty, that is, it does not
   contain any atoms, residues, chains, bonds or torsions. To populate the
   entity, use an :doc:`entity editor <editors>`.
   
   :returns: The newly created :class:`EntityHandle`

Entity Handle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: EntityHandle

  The entity class represents a molecular structure. Such a structure is in
  general made up of one or more chains of residues, which in turn are
  formed by one or more atoms.

  The interface of entities is tailored to biological macromolecules, but this
  does not prevent it to be used for molecules in general: An entity also
  represent a ligand or a collection of water molecules - hence the rather
  generic name.

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.

  .. attribute:: chains
    
    List of all chains of this entity. The chains are in the same 
    order they have been added, for example the order they have been listed 
    in a PDB file.
    
    Also available as :meth:`GetChainList`. To access a single chain, use 
    :meth:`FindChain`.
    
    This property is read-only.

    :type: :class:`ChainHandleList` (list of :class:`ChainHandle`)

  .. attribute:: chain_count
    
    Number of chains. Read-only. See :meth:`GetChainCount`.
    
    :type: :class:`int`
    
  .. attribute:: residues
    
    List of all residues of this entity. The returned residues are first 
    sorted by chain and then from N- to C-terminus.
    
    This property is read-only.
    
    **Example usage:**
    
    .. code-block:: python
    
      for res in ent.residues:
        print(res.name, res.atom_count)
    
    Also available as :meth:`GetResidueList`. To access a single residue, use 
    :meth:`FindResidue`. 
  
    :type: :class:`ResidueHandleList` (list of :class:`ResidueHandle`)
    
  .. attribute:: residue_count
    
    Number of residues. Read-only. See :meth:`GetResidueCount`.
    
    :type: :class:`int`
  
  .. attribute:: atoms
  
     Get list of all atoms of this entity. To access a single atom, use
     :meth:`FindAtom`.
     
     This property is read-only. Also available as :meth:`GetAtomList`
     
     :type: :class:`AtomHandleList` (list of :class:`AtomHandle`)
    
  .. attribute:: atom_count
    
    Number of atoms. Read-only. See :meth:`GetAtomCount`.
    
    :type: :class:`int`
  
  .. attribute:: bounds
  
    Axis-aligned bounding box of the entity. Read-only
    
    :type: :class:`ost.geom.AlignedCuboid`

  .. attribute:: mass
    
    The total mass of the entity in Dalton. Also available as :meth:`GetMass`.
    
    :type: float

  .. attribute:: center_of_atoms

    Center of atoms, that is the average atom position of the entity.
    Use :attr:`center_of_mass` for the the mass-weighted center of the entity.
    Also available as :meth:`GetCenterOfAtoms`.

    :type: :class:`~ost.geom.Vec3`

  .. attribute:: center_of_mass
  
    Center of mass of the entity.
    Also available as :meth:`GetCenterOfMass`
    
    :type: :class:`~ost.geom.Vec3`

  .. attribute:: geometric_center

    Mid-point of the axis aligned bounding box of the entity.
    Also available as :meth:`GetGeometricCenter`

    :type: Vec3

  .. attribute:: positions

    Equivalent to calling :meth:`GetPositions` with *sort_by_index = True*. This
    property is read-only and only available if OpenStructure was compiled with
    an enabled ``USE_NUMPY`` flag (see :ref:`here <cmake-flags>` for details).

    :type: :class:`numpy.array`

  .. attribute:: valid

    Validity of handle.

    :type: bool

  .. method:: GetName()

    :returns: Name associated to this entity.
    :rtype:   :class:`str`

  .. method:: SetName(name)

    :param name: Sets this as new name to be associated to this entity.
    :type name:  :class:`str`

  .. method:: FindChain(chain_name)
    
    Get chain by name. See also :attr:`chains`
    
    :param chain_name:  Chain identifier, e.g. "A"
    :type  chain_name:  str
    :returns:           A valid :class:`ChainHandle`, if the entity contains a
                        chain with the given name, an invalid
                        :class:`ChainHandle` otherwise.
                        
  .. method:: GetChainList()
    
    See :attr:`chains`

  .. method:: GetChainCount()

    See :attr:`chain_count`
    
  .. method:: FindResidue(chain_name, res_num)
    
    Get residue by chain name and residue number. See also 
    :meth:`GetResidueList`
    
    :param chain_name:  Chain identifier, e.g. "A"
    :type  chain_name:  str
    :param    res_num:  residue number
    :type     res_num:  :class:`ResNum`
    
    :returns:           A valid :class:`ResidueHandle` if the chain exists and
                        the chain contains a residue of the given residue
                        number, an invalid :class:`ResidueHandle` otherwise.
    
  .. method:: GetResidueList()
    
    See :attr:`residues`

  .. method:: GetResidueCount()
    
    See :attr:`residue_count`
    
  .. method:: FindAtom(chain_name, res_num, atom_name)
    
    Get atom by chain name, residue number and atom name. See also
    :attr:`atoms`
    
    :param chain_name:  Chain identifier, e.g. "A"
    :type  chain_name:  str
    :param    res_num:  residue number
    :type     res_num:  :class:`ResNum`
    :param  atom_name:  atom name, e.g. CA
    :type   atom_name:  str
    
    :returns:           A valid :class:`AtomHandle` if an atom matching the
                        parameters could be found, an invalid
                        :class:`AtomHandle` otherwise
    
  .. method:: GetAtomList()
    
    See :attr:`atoms`
    
  .. method:: GetAtomCount()
    
    See :attr:`atom_count`
    
  .. method:: EditXCS([edit_mode=mol.EditMode.UNBUFFERED_EDIT])
    
    Request :class:`XCSEditor` for editing the external coordinate system. This
    call will fail when there are pending changes of the internal coordinate
    system.
    
    :param edit_mode: Must be EditMode.BUFFERED_EDIT or
                      EditMode.UNBUFFERED_EDIT. For more details, see the
                      editor documentation.
    :type edit_mode: mol.EditMode
    
    :returns: :class:`XCSEditor`
    
  .. method:: EditICS([edit_mode=mol.EditMode.UNBUFFERED_EDIT])
    
    Request :class:`ICSEditor` for editing the internal coordinate system, such
    as torsions, bond lengths and angle between two bonds. This call will fail
    when there are pending changes of the external coordinate system.
    
    :param edit_mode: Must be EditMode.BUFFERED_EDIT or
                      EditMode.UNBUFFERED_EDIT. For more details, see the
                      editor documentation.
    :type edit_mode: mol.EditMode
    
    :returns: :class:`ICSEditor`
    
  .. method:: Select(query, flags=0)
    
    Perform a selection on the entity. The result of the selection is an 
    :class:`EntityView` which (usually) contains only a subset of chains,
    residues, atoms and bonds of the original entity.
    
    **Example Usage:**
    
    .. code-block:: python

      # select calpha atoms of peptides
      calphas = ent.Select('aname=CA and peptide=true')
      
      # select atoms in a box of size 10, centred at the origin
      in_box = ent.Select('x=-5:5 and y=-5:5 and z=-5:5')
    
    :param query: The query to be executed.
    :type  query: :class:`Query` / :class:`str`
    :param flags: An ORed together combination of :class:`QueryFlag`
    :type  flags: :class:`int` / :class:`QueryFlag`
    :returns:     An :class:`entity view <EntityView>`.
    :raises:      :class:`QueryError` when the query could not be executed due
                  to syntactic errors.
    
  .. method:: CreateFullView()
    
    Creates  an entity view containing all chains, residues, atoms and bonds of
    this entity.
    
    .. code-block:: python
    
      # these two lines are identical
      full=ent.Select('')
      full=ent.CreateFullView()
    
    :returns: :class:`EntityView`
    
  .. method:: CreateEmptyView()
    
    Creates an entity view pointing to this entity, but otherwise empty. This
    method is usually the starting point for manual entity view creation, e.g.:
    
    .. code-block:: python
    
      view=ent.CreateEmtpyView()
      for atom in ent.atoms:
        if ComplicatedPredicate(atom):
           view.AddAtom(atom)
    
    :returns: :class:`EntityView`

  .. method:: Copy()
    
    Creates a deep copy of the entity.
    
    :returns: A new :class:`EntityHandle` that is an exact copy of this entity.
    
    .. note::
      
      alternative atom positions are not handled yet.

  .. method:: GetMass()

    See :attr:`mass`

  .. method:: GetCenterOfAtoms()

    See :attr:`center_of_atoms`

  .. method:: GetCenterOfMass()

    See :attr:`center_of_mass`

  .. method:: GetGeometricCenter()

    See :attr:`geometric_center`

  .. method:: GetPositions(sort_by_index=True)

    :return: Array of atom positions for this entity.
    :rtype:  :class:`numpy.array` (shape [:attr:`atom_count`, 3])
    :param sort_by_index: If True, the atoms are sorted by their
                          :attr:`~AtomHandle.index`. Otherwise, they are sorted
                          as they appear in the :attr:`atoms` list.

    This method is only available if OpenStructure was compiled with an enabled
    ``USE_NUMPY`` flag (see :ref:`here <cmake-flags>` for details).
    
  .. method:: FindWithin(pos, radius)
  
    Find all atoms in sphere of given radius centred at *pos*. To turn the
    returned list of atoms into an :class:`EntityView`, use
    :func:`CreateViewFromAtomList`.
    
    :param pos: Center of sphere
    :type pos: :class:`~ost.geom.Vec3`
    :param radius: The radius of the sphere
    :type radius: float
    
    :returns: :class:`AtomHandleList` (list of :class:`AtomHandle`)

  .. method:: IsValid()
  
    See :attr:`valid`

Chain Handle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: ChainHandle

  A chain of one or more :class:`residues <ResidueHandle>`. Chains are always 
  part of an entity (see :ref:`note on memory management <memory_management>`).

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.

  .. attribute:: atoms

     Get list of all atoms of this chain. To access a single atom, use
     :meth:`FindAtom`.

     This property is read-only. Also available as :meth:`GetAtomList`

     :type: :class:`AtomHandleList` (list of :class:`AtomHandle`)

  .. attribute:: bounds

    Axis-aligned bounding box of the chain. Read-only

    :type: :class:`ost.geom.AlignedCuboid`

  .. attribute:: center_of_atoms

    Center of atoms (not mass weighted). Also available as
    :meth:`GetCenterOfAtoms`.

    :type: :class:`~ost.geom.Vec3`

  .. attribute:: center_of_mass

    Center of mass. Also available as :meth:`GetCenterOfMass`.

    :type: :class:`~ost.geom.Vec3`

  .. attribute:: geometric_center

    Mid-point of the axis aligned bounding box of the chain.
    Also available as :meth:`GetGeometricCenter`

    :type: Vec3

  .. attribute:: description

     Details about the chain. Not categorised, just text.

  .. attribute:: in_sequence

     Whether the residue numbers are in ascending order. For example:

     .. code-block:: python

       chain=ent.FindChain("A")
       print(chain.residues) # [A.GLY1, A.GLY2, A.GLY4A, A.GLY4B]
       print(chain.in_sequence) # prints true

       chain=ent.FindChain("B")
       print(chain.residues) # [B.GLY1, B.GLY4, B.GLY3]
       print(chain.in_sequence) # prints false

  .. attribute:: is_oligosaccharide

    Indicates if the chain is an oligosaccharide, a branched, non-linear entity
    of multiple sugars. Also available as :meth:`IsOligosaccharide`.

    :type: :class:`bool`

  .. attribute:: is_polymer

    Indicates if a chain is a polymer. True for polypeptides, polynucleotides,
    polysaccharides, oligosaccharides and branched chains. Also available as
    :meth:`IsPolymer`.

    :type: :class:`bool`

  .. attribute:: is_polynucleotide

    Indicates if a chain is a nucleic acid. Also available as
    :meth:`IsPolynucleotide`.

    :type: :class:`bool`

  .. attribute:: is_polypeptide

    Indicates if a chain is a protein. Also available as :meth:`IsPolypeptide`.

    :type: :class:`bool`

  .. attribute:: is_polysaccharide

    Indicates if a chain is a polysaccharide. Also available as
    :meth:`IsPolysaccharide()`.

    :type: :class:`bool`

  .. attribute:: mass

    The total mass of this chain in Dalton. Also available as :meth:`GetMass`

    :type: float

  .. attribute:: name
  
     The chain name. The name uniquely identifies the chain in the entity. In 
     most cases, the chain name is one character. This is restriction of the PDB 
     file format. However, you are free to use longer names as long as you don't 
     want to save them as PDB files
     
     This property is read-only. To change the name, use an :class:`XCSEditor`. 

     Also available as :meth:`GetName`
     
     :type: str

  .. attribute:: residue_count

    Number of residues. Read-only. See :meth:`GetResidueCount`.

    :type: :class:`int`

  .. attribute:: residues
   
     List of all residues of this chain. The residues are sorted from N- to 
     C-terminus. Usually the residue numbers are in ascending order 
     (see :attr:`in_sequence`).
   
     This property is read-only.
   
     **Example usage:**
   
     .. code-block:: python
     
       chain=ent.FindChain("A")
       for res in chain.residues:
         print(res.name, res.atom_count)
   
     Also available as :meth:`GetResidueList`. To access a single residue, use 
     :meth:`FindResidue`. 

     :type: :class:`ResidueHandleList` (list of :class:`ResidueHandle`)
  
  .. attribute:: type

     Describes the type of the chain.

     :type: :class:`ChainType`.

  .. attribute:: hash_code

    A unique identifier for this chain. Note that a deep copy of an entity (see
    :meth:`EntityHandle.Copy`) will have chains with differing identifiers.
    Shallow copies of the entity preserve the identifier.
    Also available as :meth:`GetHashCode`.

    :type: int

  .. attribute:: valid

    Validity of handle.

    :type: bool
    
  .. method:: FindResidue(res_num)
   
    Get residue by residue number. See also :attr:`residues`.
    
    :param    res_num:  residue number
    :type     res_num:  :class:`ResNum`
    
    :returns:           A valid :class:`ResidueHandle` if the chain contains
                        a residue with matching residue number, an invalid
                        :class:`ResidueHandle` otherwise.
                        
  .. method:: GetResidueList()

    See :attr:`residues`

  .. method:: GetResidueCount()

    See :attr:`residue_count`

  .. method:: FindAtom(res_num, atom_name)

    Get atom by residue number and atom name. See also :attr:`atoms`

    :param    res_num:  residue number
    :type     res_num:  :class:`ResNum`
    :param  atom_name:  atom name, e.g. CA
    :type   atom_name:  str

    :returns:           A valid :class:`AtomHandle` if an atom matching the
                        parameters could be found, an invalid
                        :class:`AtomHandle` otherwise

  .. method:: GetAtomList()
    
    See :attr:`atoms`

  .. method:: GetName()
  
    See :attr:`name`

  .. method:: GetType()

    See :attr:`type`

  .. method:: IsOligosaccharide()

    See :attr:`is_oligosaccharide`

  .. method:: IsPolymer()

    See :attr:`is_polymer`

  .. method:: IsPolynucleotide()

    See :attr:`is_polynucleotide`

  .. method:: IsPolypeptide()

    See :attr:`is_polypeptide`

  .. method:: IsPolysaccharide()

    See :attr:`is_polysaccharide`

  .. method:: GetDescription()

    See :attr:`description`

  .. method:: GetHashCode()

    See :attr:`hash_code`

    :rtype: int

  .. method:: IsValid()
  
    See :attr:`valid`

Residue Handle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: ResidueHandle

  The residue is either used to represent complete molecules or building blocks 
  in a polymer, e.g. in a protein, DNA or RNA. A residue consists of one or 
  more :class:`atoms <AtomHandle>`. Residues are always part of a 
  :class:`ChainHandle`, even if they are ligands or water molecules where the 
  concept of a chain does not apply, and of an entity (see
  :ref:`note on memory management <memory_management>`).

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.
  
  .. attribute:: name
  
    The residue name is usually a str of 3 characters, e.g. `GLY` for 
    glycine or `ALA` for alanine, but may be shorter, e.g. `G` for guanosine, 
    or longer for structures loaded from formats other than PDB.
    Also available as :meth:`GetName`.
    
    This property is read-only. To change the name of the residue, use
    :meth:`~EditorBase.RenameResidue`.

    :type: str

  .. attribute:: qualified_name

    The qualified name consists of a residue identifier and chain name.
    For a glycine with residue number 2 of chain A, the qualified name is
    "A.GLY2". Also available as :meth:`GetQualifiedName`.

    :type: str
  
  .. attribute:: number
  
    The number of this residue. The residue number has a numeric part and an
    insertion-code. This property is read-only. Also available as 
    :meth:`GetNumber`.
    
    :type: :class:`ResNum`
  
  .. attribute:: one_letter_code
  
    For amino acids, and nucleotides the `one_letter_code` is an alpha-numeric 
    character. For unknown or more *exotic* residues, the one letter code is set 
    to '?'.
    
    **Example**
    
    This code-snippet shows how to get the sequence string from a list of 
    residues.
    
    .. code-block:: python
    
      print(''.join([r.one_letter_code for r in chain.residues]))
    
    :type: str
    
  .. attribute:: atoms

     Get list of all atoms of this residue. To access a single atom, use
     :meth:`FindAtom`.
   
     This property is read-only. Also available as :meth:`GetAtomList`

     :type: :class:`AtomHandleList` (list of :class:`AtomHandle`)

  .. attribute:: atom_count

    Number of atoms. Read-only. See :meth:`GetAtomCount`.

    :type: :class:`int`

  .. attribute:: bond_count

    Number of bonds. Read-only. See :meth:`GetBondCount`.

    :type: :class:`int`

  .. attribute:: bounds
  
    Axis-aligned bounding box of the residue. Read-only.
    
    :type: :class:`ost.geom.AlignedCuboid`
    
  .. attribute:: mass
  
    The total mass of this residue in Dalton. Also available as :meth:`GetMass`.
  
    :type: float
  
  .. attribute:: center_of_atoms
    
    Center of atoms (not mass weighted). Also available as 
    :meth:`GetCenterOfAtoms`.
    
    :type: :class:`~ost.geom.Vec3`

  .. attribute:: center_of_mass

    Center of mass. Also available as :meth:`GetCenterOfMass`

    :type: :class:`~ost.geom.Vec3`

  .. attribute:: geometric_center

    Mid-point of the axis aligned bounding box of the residue.

    :type: Vec3

  .. attribute:: chain
  
    The chain this residue belongs to. Read-only. Also available as 
    :meth:`GetChain`
    
    :type: :class:`ChainHandle`
  
  .. attribute:: phi_torsion
  
    The PHI dihedral angle between this residue and the next. For residues 
    that are not amino acids, residues that do not have all atoms required or 
    residues that do not have bonds between the four atoms involved in the 
    torsion, the PHI torsion is an invalid handle.
    
    Read-only. Also available as :meth:`GetPhiTorsion`
    
    :type: :class:`TorsionHandle`
  
  .. attribute:: psi_torsion
  
    The PSI dihedral angle between this residue and the previous. For residues 
    that are not amino acids, residues that do not have all atoms required or 
    residues that do not have bonds between the four atoms involved in the 
    torsion, the PSI torsion is an invalid handle.
    
    Read-only. Also available as :meth:`GetPsiTorsion`

  .. attribute:: omega_torsion

    The OMEGA dihedral angle between this residue and the previous. For residues
    that are not amino acids, residues that do not have all atoms required or
    residues that do not have bonds between the four atoms involved in the
    torsion, the OMEGA torsion is an invalid handle.

    Read-only. Also available as :meth:`GetOmegaTorsion`
    
    :type: :class:`TorsionHandle`
  
  .. attribute:: chem_class
  
    The chemical class of the residue.

    :type: :class:`ChemClass`

  .. attribute:: chem_type

    The chemical type of the residue. The type is only properly set if a
    compound library is used.

    :type: :class:`ChemType`
  
  .. attribute:: sec_structure
  
    The secondary structure of the residue.

    :type: :class:`SecStructure`
  
  .. attribute:: is_ligand

    Whether the residue is a ligand. When loading PDB structures, this property
    is set based on the HET records. This also means, that this property will
    most likely not be set properly for all except PDB files coming from
    pdb.org. When loading MMCIF structures, this property is set based on
    `seqres` information and differs from PDB structures.
    Also available as :meth:`IsLigand`, :meth:`SetIsLigand`.

    :type: bool
  
  .. attribute:: is_protein
  
    Whether the residue is considered to be part of a protein. This is set when
    loading a structure if the residue forms a feasible peptide bond to the
    previous or next residue (see :meth:`~ost.conop.IsBondFeasible`). Also
    available as :meth:`IsProtein`, :meth:`SetIsProtein`. In contrast to
    :meth:`IsPeptideLinking` this excludes residues which are not connected to
    neighbouring residues such as CA-only residues or badly positioned ones.

    :type: bool

  .. attribute:: peptide_linking
  
    Whether residue can form peptide bonds. This is determined based on
    :attr:`chem_class` which is set when loading the structure.

    :type: :class:`bool`

  .. attribute:: nucleotide_linking
  
    Whether residue can form nucleotide bonds. This is determined based on
    :attr:`chem_class` which is set when loading the structure.

    :type: :class:`bool`

  .. attribute:: index

    Residue index (starting at 0) within chain.

  .. attribute:: central_atom

    Central atom used for rendering traces. For peptides, this is usually
    the CA atom. For nucleotides, this is usually the P atom.

    :type: :class:`AtomHandle`

  .. attribute:: central_normal

    Normal computed for :attr:`central_atom`. Only defined for peptides and
    nucleotides if all required atoms available. Otherwise, the (1,0,0) vector
    is returned.

    :type: :class:`~ost.geom.Vec3`

  .. attribute:: hash_code

    A unique identifier for this residue. Note that a deep copy of an entity (see
    :meth:`EntityHandle.Copy`) will have residues with differing identifiers.
    Shallow copies of the entity preserve the identifier. Also available as
    :meth:`GetHashCode`.

    :type: int

  .. attribute:: valid

    Validity of handle.

    :type: bool

  .. attribute:: next

    Residue after this one in the same chain. Invalid handle returned if there
    is no next residue. Residues are ordered as in :attr:`ChainHandle.residues`
    independently on whether they are connected or not (see
    :func:`InSequence` to check for connected residues).

    :type: :class:`ResidueHandle`

  .. attribute:: prev

    Residue before this one in the same chain. Otherwise same behaviour as
    :attr:`next`.

    :type: :class:`ResidueHandle`

  .. method:: FindAtom(atom_name)

    Get atom by atom name. See also :attr:`atoms`

    :param  atom_name:  atom name, e.g. CA
    :type   atom_name:  str

    :returns:           A valid :class:`AtomHandle` if an atom with the given
                        name could be found, an invalid :class:`AtomHandle`
                        otherwise

  .. method:: HasAltAtoms()

    Returns whether the residue has any alternative atom groups

    :rtype:             :class:`bool`

  .. method:: HasAltAtomGroup(group_name)

    Returns whether the residue has a certain alternative atom group

    :param group_name:  Group of interest
    :type group_name:   :class:`str`
    :rtype:             :class:`bool`

  .. method:: GetAltAtomGroupNames()

    Returns names of all alternative atom groups in residue

    :rtype:             :class:`list` of :class:`str`

  .. method:: GetCurrentAltGroupName()

    Returns the currently active alternative atom group, empty
    :class:`str` if residue has no alternative atom groups

    :rtype:             :class:`str`

  .. method:: SwitchAtomPos(group_name)

    Switches the atoms in residue to the specified alternative atom group

    :param group_name:  Group of interest
    :type group_name:   :class:`str`
    :rtype:             :class:`bool`
    :returns:           Whether the switch was successful (e.g. False if no such
                        group exists)

  .. method:: GetName()

    See :attr:`name`

  .. method:: GetQualifiedName()

    See :attr:`qualified_name`

  .. method:: GetNumber()

    See :attr:`number`

  .. method:: GetOneLetterCode()

    See :attr:`one_letter_code`

  .. method:: GetAtomList()

    See :attr:`atoms`

  .. method:: GetAtomCount()

    See :attr:`atom_count`

  .. method:: GetBondCount()

    See :attr:`bond_count`

  .. method:: GetBounds()

    See :attr:`bounds`

  .. method:: GetMass()

    See :attr:`mass`

  .. method:: GetCenterOfAtoms()

    See :attr:`center_of_atoms`

  .. method:: GetCenterOfMass()

    See :attr:`center_of_mass`

  .. method:: GetGeometricCenter()

    See :attr:`geometric_center`

  .. method:: GetChain()

    See :attr:`chain`

  .. method:: GetPhiTorsion()

    See :attr:`phi_torsion`

  .. method:: GetPsiTorsion()

    See :attr:`psi_torsion`

  .. method:: GetOmegaTorsion()

    See :attr:`omega_torsion`

  .. method:: GetChemClass()

    See :attr:`chem_class`

  .. method:: GetChemType()

    See :attr:`chem_type`

  .. method:: GetSecStructure()

    See :attr:`sec_structure`

  .. method:: IsLigand()

    See :attr:`is_ligand`

  .. method:: IsProtein()

    See :attr:`is_protein`

  .. method:: IsPeptideLinking()

    See :attr:`peptide_linking`

  .. method:: IsNucleotideLinking()

    See :attr:`nucleotide_linking`

  .. method:: GetIndex()
    
    See :attr:`index`

  .. method:: GetCentralAtom()
              SetCentralAtom()
    
    See :attr:`central_atom`

  .. method:: GetCentralNormal()
    
    See :attr:`central_normal`

  .. method:: GetHashCode()

    See :attr:`hash_code`

    :rtype: int

  .. method:: IsValid()
  
    See :attr:`valid`

  .. method:: SetIsLigand()

    Set the :meth:`IsLigand` flag explicitly.

    :param ligand: Whether this residue is a ligand or not
    :type  ligand: bool

  .. method:: SetIsProtein()

    Set the :meth:`IsProtein` flag explicitly.

    :param protein: Whether this residue is a protein or not
    :type  protein: bool

  .. method:: GetNext()

    See :attr:`next`

  .. method:: GetPrev()

    See :attr:`next`

Atom Handle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: AtomHandle

  Represents an atom in a molecular structure. Atoms are always part of a 
  residue and of an entity (see
  :ref:`note on memory management <memory_management>`)

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.
  
  .. attribute:: name
  
    The atom name, usually a capitalized string consisting of the element and
    an position indicator, e.g. `G1`. In general, the atom name uniquely 
    identifies an atom within a residue. However, this is not enforced and there
    may be more than one atoms with the same name. 
    
    :type: str
  
  .. attribute:: qualified_name
  
     The qualified name consists of the atom name as well as a residue
     identifier and chain name. For CA of a glycine with residue number 2 of
     chain A, the qualified name is "A.GLY2.CA".
     
     :type: str
     
  .. attribute:: element
  
    The atom's element. By convention in Openstructure, this is the chemical
    symbol in uppercase, but this is not strictly enforced and may be a non-
    existing element or an empty string. Also available as :meth:`GetElement`.
    Read-only.
    
    :type: str
    
  .. attribute:: mass
  
    The atom's mass in Dalton. Also available as :meth:`GetMass`. Read-only.
    
    :type: float
    
  .. attribute:: pos
  
    The atom's position in global coordinates, transformed by the entity 
    transformation. Also available as :meth:`GetPos`.
    
    Read-only.
    
    :type: :class:`~ost.geom.Vec3`
    
  .. attribute:: original_pos
  
    The atom's untransformed position in global coordinates. Also available as 
    :meth:`GetOriginalPos`.
    
    Read-only.
    
    :type: :class:`~ost.geom.Vec3`
    
  .. attribute:: radius
  
    The van-der-Waals radius of the atom. Also available as :meth:`GetRadius`. 
    Read/write.
    
    :type: float
  
  .. attribute:: occupancy
  
    The atom's occupancy in the range 0 to 1. Read/write. Also available as 
    :meth:`GetOccupancy`, :meth:`SetOccupancy`.

    :type: float
    
  .. attribute:: b_factor
  
    The atom's temperature factor. Read/write. Also available as 
    :meth:`GetBFactor`, :meth:`SetBFactor`.
    
    :type: float    

  .. attribute:: charge
    
    The atom's charge. Also available as :meth:`GetCharge`, :meth:`SetCharge`.
    
    :type: float

  .. attribute:: residue
  
    The residue this atom belongs to. Read-only.
    
    :type: :class:`ResidueHandle`
  
  .. attribute:: is_hetatom
  
    Indicates whether this atom is a *HETATM*. When loading a structure from a 
    PDB file, all non-standard aminoacid and nucleotide atoms as well as ligands 
    are marked as HETATM. 
    
  .. attribute:: bonds
    
    List of bonds this atom is involved in. The bonds are  in no particular 
    order. Also available as :meth:`GetBondList`.
    
    :type: list of :class:`bond handles<BondHandle>`
  
  .. attribute:: index

    Atom index (starting at 0) within entity.

    :type: int

  .. attribute:: hash_code

    A unique identifier for this atom. Note that a deep copy of an entity (see
    :meth:`EntityHandle.Copy`) will have atoms with differing identifiers.
    Shallow copies of the entity preserve the identifier. Atom views on a handle
    have different identifiers, but the atom view handles (see
    :attr:`AtomView.handle`) have the same identifier. Also available
    as :meth:`GetHashCode`.

    :type: int

  .. attribute:: valid

    Validity of handle.

    :type: bool

  .. method:: FindBondToAtom(other_atom)

    Finds and returns the bond formed between this atom and `other_atom`. If no 
    bond exists between the two atoms, an invalid :class:`BondHandle` is 
    returned.
    
    :param other_atom: The other atom
    :type  other_atom: :class:`AtomHandle`
    :rtype: :class:`BondHandle`

  .. method:: GetBondCount()
    
    :rtype: int

  .. method:: GetBondList()
    
    See :attr:`bonds`
    
    :rtype: :class:`BondHandleList` (list of :class:`BondHandle`)

  .. method:: GetBondPartners()
    
    Get list of atoms this atom is bonded with. This method exists as a shortcut 
    to determine all the bonding partners of an atom and avoids code like this: 
    
    .. code-block:: python
    
      bond_parters=[]
      for bond in atom.bonds:
        if bond.first==atom:
          bond_partners.append(bond.second)
        else:
          bond_partners.append(bond.first)

    :rtype: :class:`AtomHandleList`

  .. method:: GetCharge()
  
    See :attr:`charge`
  
    :rtype: float

  .. method:: GetElement()
    
    See :attr:`element`
    
    :rtype: str

  .. method:: GetEntity()
    
    The entity this atom belongs to
    
    :rtype: :class:`EntityHandle`

  .. method:: GetHashCode()
    
    See :attr:`hash_code`
    
    :rtype: int

  .. method:: GetIndex()
    
    See :attr:`index`
    
    :rtype: int

  .. method:: GetMass()
    
    See :attr:`mass`
    
    :rtype: float

  .. method:: GetName()

    See :attr:`name`
    
    :rtype: str

  .. method:: GetOriginalPos()

    :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetPos()

    See :attr:`pos`
    
    :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetQualifiedName()
    
    See :attr:`qualified_name`
    
    :rtype: str

  .. method:: GetRadius()
    
    See :attr:`radius`
    
    :rtype: float

  .. method:: GetResidue()
    
    See :attr:`residue`
    
    :rtype: :class:`ResidueHandle`

  .. method:: IsHetAtom()
    
    See :attr:`is_hetatom`
    
    :rtype: bool

  .. method:: IsValid()
  
    See :attr:`valid`

Bond Handle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: BondHandle

  Represents a chemical bond between two atoms (first and second).

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.

  .. attribute:: first
                 second

    Atoms involved in the bond. No assumptions about the order should be made.
    With the internal coordinate system enabled, first and second may even be
    swapped when rebuilding the internal connectivity tree. Also available as
    :meth:`GetFirst` and :meth:`GetSecond`.

    :type: :class:`AtomHandle`

  .. attribute:: pos

    Midpoint between the two atoms (transformed coordinates). Also available as
    :meth:`GetPos`.

    :type: :class:`~ost.geom.Vec3`

  .. attribute:: length

    Length of the bond. Also available as :meth:`GetLength`.

    :type: float

  .. attribute:: bond_order

    The bond order. Possible values:

      * ``1`` - single bond
      * ``2`` - double bond
      * ``3`` - triple bond
      * ``4`` - aromatic bond

    Also available as :meth:`GetBondOrder`.

    :type: int

  .. attribute:: hash_code

    A unique identifier for this bond handle.  Also available as :meth:`GetHashCode`.

    :type: int

  .. attribute:: valid

    Validity of handle.  Also available as :meth:`IsValid`.

    :type: bool

  .. method:: GetFirst()

    See :attr:`first`

  .. method:: GetSecond()

    See :attr:`second`

  .. method:: GetPos()

    See :attr:`pos`

  .. method:: GetLength()

    See :attr:`length`

  .. method:: GetBondOrder()

    See :attr:`bond_order`

  .. method:: GetOther(other_atom)

    Get the other atom. Returns the one of the two atoms that does not match
    the given one.

    :param other_atom: The other atom
    :type  other_atom: :class:`AtomHandle`
    :rtype: :class:`AtomHandle`

  .. method:: SetBondOrder(order)

    Set the bond order. See :meth:`GetBondOrder`.

    :param order: The bond order
    :type  order: :class:`int`

    See :attr:`bond_order`

  .. method:: IsValid()

    See :attr:`valid`

  .. method:: GetHashCode()

    See :attr:`hash_code`


The View Classes
--------------------------------------------------------------------------------

:class:`Entity <EntityView>`, :class:`residue <ResidueView>` and
:class:`atom <AtomView>` views can store arbitrary
:doc:`generic properties <../../base/generic>`.

Entity View
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: EntityView

  An entity view represents a structural subset of an :class:`EntityHandle` and
  contains :class:`ChainView`\s, :class:`ResidueView`\s, :class:`AtomView`\s and
  :class:`BondHandle`\s.
  For an introduction, see :doc:`../../intro-01`.

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.

  .. attribute:: handle

     The underlying :class:`handle <EntityHandle>` of the entity view. Also
     available as :meth:`GetHandle`.

     :type: :class:`EntityHandle`
  
  .. attribute:: chains
   
    List of all chains of this entity. The chains are in the same 
    order they have been added.
    
    Also available as :meth:`GetChainList`. To access a single chain, use 
    :meth:`FindChain`.
    
    This property is read-only.

    :type: :class:`ChainViewList` (list of :class:`ChainView`)

  .. attribute:: chain_count
    
    Number of chains. Read-only. See :meth:`GetChainCount`.
    
    :type: :class:`int`
    
  .. attribute:: residues
    
    List of all residues of this entity. The returned residues are first 
    sorted by chain and then from N- to C-terminus.
    
    This property is read-only.
    
    **Example usage:**
    
    .. code-block:: python
    
      for res in ent.residues:
        print(res.name, res.atom_count)
    
    Also available as :meth:`GetResidueList`. To access a single residue, use 
    :meth:`FindResidue`. 
    
    :type: :class:`ResidueViewList` (list of :class:`ResidueView`)
    
  .. attribute:: residue_count
    
    Number of residues. Read-only. See :meth:`GetResidueCount`.
    
    :type: :class:`int`

  .. attribute:: atoms

    Get list of all atoms of this entity. To access a single atom, use
    :meth:`FindAtom`.
    
    This property is read-only. Also available as :meth:`GetAtomList`
    
    :type: :class:`AtomViewList` (list of :class:`AtomView`)
    
  .. attribute:: atom_count
    
    Number of atoms. Read-only. See :meth:`GetAtomCount`.
    
    :type: :class:`int`

  .. attribute:: bounds
  
    See :attr:`EntityHandle.bounds`

  .. attribute:: mass

    See :attr:`EntityHandle.mass`

  .. attribute:: center_of_atoms

    See :attr:`EntityHandle.center_of_atoms`

  .. attribute:: center_of_mass

    See :attr:`EntityHandle.center_of_mass`

  .. attribute:: geometric_center

    See :attr:`EntityHandle.geometric_center`

  .. attribute:: valid

     Validity of view.

     :type: bool

  .. method:: GetName()

    :returns: :func:`~EntityHandle.GetName` of entity :attr:`handle`.
    :rtype:   :class:`str`

  .. method:: SetName(name)

    :param name: Passed to :func:`~EntityHandle.SetName` of :attr:`handle`.
    :type name:  :class:`str`

  .. method:: CreateEmptyView()
  
    See :meth:`EntityHandle.CreateEmptyView`

    :rtype: :class:`EntityView`

  .. method:: CreateFullView()

    Returns a copy of this entity. Provided for `duck-typing <http://en.wikipedia.org/wiki/Duck_typing>` purposes.
    
    :rtype: :class:`EntityView`

  .. method:: AddChain(chain_handle[, view_add_flags])

    Add chain to view. By default, only the chain is added to the view, but not 
    its residues and atoms. This behaviour can be changed by passing in an
    appropriate set of `view_add_flags`.

    The following example just adds a chain without residues and atoms:

    .. code-block:: python

        pdb = ost.io.LoadPDB(<PDB file name>)
        v = pdb.CreateEmptyView()
        v.AddChain(pdb.chains[0])

    To get the chain with residues and atoms, go like:

    .. code-block:: python

        pdb = ost.io.LoadPDB(<PDB file name>)
        v = pdb.CreateEmptyView()
        v.AddChain(pdb.chains[0], ost.mol.INCLUDE_ALL)

    Note that the view above still lacks bonds which can be added with the
    :meth:`AddAllInclusiveBonds` method.

    :param chain_handle: The chain handle to be added.
    :type  chain_handle: :class:`ChainHandle`
    :param view_add_flags: An ORed together combination of :class:`ViewAddFlag`
    :type view_add_flags: :class:`int` / :class:`ViewAddFlag`
    :rtype: :class:`ChainView`

  .. method:: AddResidue(residue_handle[, view_add_flags])

    Add residue to view. If the residue's chain is not already part of the
    view, it will be added. By default, only the residue is added, but not its
    atoms. This behaviour can be modified by passing in an appropriate
    combination of :class:`ViewAddFlag`.

    :param residue_handle: The residue handle to be added
    :type  residue_handle: :class:`ResidueHandle`
    :param view_add_flags: An ORed together combination of :class:`ViewAddFlag`
    :type  view_add_flags: :class:`int` / :class:`ViewAddFlag`
    :rtype: :class:`ResidueView`

  .. method:: AddAtom(atom_handle[, view_add_flags])

    Add the atom to view. The chain and residue of the atom are added to the 
    view if they haven't already been added.
    
    :param atom_handle: The atom handle
    :type  atom_handle: :class:`AtomHandle`
    :param view_add_flags: An ORed together combination of :class:`ViewAddFlag`
    :type  view_add_flags: :class:`int` / :class:`ViewAddFlag`
    :rtype: :class:`AtomView`

  .. method:: AddBond(bond_handle)

    Add bond to the view.
    
    :param bond_handle: Bond to be added
    :type  bond_handle: BondHandle
    :returns: True if the bond has been added, false if the bond was already 
       present in the view.
    :rtype: bool

  .. method:: AddAllInclusiveBonds()
  
     Convenience method to add all bonds to the view for which both of the bond 
     partners are in the view. This method is useful when manually assembling 
     views.
     
  .. method:: RemoveChain(chain)
    
    Remove chain and all its residues from the entity view.
    
    :param chain: The chain to be removed. May be invalid
    :type  chain: ChainView

  .. method:: RemoveResidue(residue)

    Remove residue from view. 
    
    :param residue: The residue to be removed. May be invalid
    
    :type  residue: ResidueView
    
  .. method:: RemoveAtom(atom)
  
    Remove atom from view

    :param atom: The atom to be removed. May be invalid
    
    :type  atom: AtomView

  .. method:: GetAngle(atom1, atom2, atom3)

    :param atom1:
    :type  atom1: AtomHandle
    :param atom2:
    :type  atom2: AtomHandle
    :param atom3:
    :type  atom3: AtomHandle
    :rtype: float

  .. method:: FindWithin(pos, radius)

    Find all atoms that are within radius of the given position. See 
    :meth:`EntityHandle.FindWithin`.
    
    :param pos: Center of sphere
    :type pos: :class:`~ost.geom.Vec3`
    :param radius: The radius of the sphere
    :type radius: float
    
    :returns: :class:`AtomHandleList` (list of :class:`AtomHandle`)

  .. method:: FindChain(chain_name)

    Find chain by name.
    
    :param chain_name:  Chain identifier, e.g. "A"
    :type  chain_name:  str
    :returns: The chain if present in the view, an invalid :class:`ChainView` 
       otherwise
    :rtype: :class:`ChainView`

  .. method:: FindResidue(residue)
    
    Deprecated. Use :meth:`ViewForHandle` instead.

    :param residue: Residue handle
    :type  residue: ResidueHandle
    :returns: The residue view pointing the handle, or an invalid handle if  the residue is not part of the view
    :rtype: :class:`ResidueView`

  .. method:: FindResidue(chain_name, res_num)
    :noindex:

    Find residue by chain name and residue number.

    :param chain_name:  Chain identifier, e.g. "A"
    :type  chain_name:  str
    :param    res_num:  residue number
    :type     res_num:  :class:`ResNum`
    :returns: The residue if present in the view, an invalid :class:`ResidueView`
       otherwise
    :rtype: :class:`ResidueView`

  .. method:: FindAtom(chain_name, res_num, atom_name)

    :param chain_name: The chain name
    :type  chain_name: str
    :param res_num: The residue number
    :type  res_num: :class:`ResNum` or :class:`int`
    :param atom_name: The name of the atom
    :type  atom_name: str
    :returns: The atom if present in the view, an invalid :class:`AtomView`
       otherwise
    :rtype: :class:`AtomView`

  .. method:: ViewForHandle(handle)

    Find chain view, residue view or atom view of pointing to the given handle.

    :param handle: handle to search for
    :type  residue: (Chain\|Residue\|Atom)Handle
    :returns: The view pointing the handle, or an invalid handle if the handle
         is not part of the view
    :rtype: (Chain\|Residue\|Atom)View

  .. method:: ExtendViewToResidues()

    Extend current view to include all atoms of each residue where
    at least one atom is selected currently.

    :returns: The extended view
    :rtype: :class:`EntityView`

  .. method:: ExtendViewToSurrounding(handle)

    Extend current view to include all atoms that are within the sum
    of their van-der-Waals radius + gap.

    This includes all atoms within: at1.GetRadius() + at2.GetRadius() + gap.

    :param gap: the gap between atoms
    :type  gap: float
    :returns: The extended view
    :rtype: :class:`EntityView`

  .. method:: Select(query, flags=0)

    Perform selection on entity view. See :meth:`EntityHandle.Select`.
    
    :param query: The query
    :type  query: :class:`Query` / :class:`str`
    :param flags: An ORed together combination of :class:`QueryFlag`
    :type  flags: :class:`int` / :class:`QueryFlag`
    :rtype: :class:`EntityView`

  .. method:: Copy()
    
    Returns a deep copy of the entity view. The effect is identical to calling
    
    .. code-block:: python
      
      the_copy = view.Select('')
    
    :rtype: :class:`EntityView`

  .. method:: Dump()

    Returns a string containing a human-readable summary of the entity view.

    :rtype: :class:`str`

  .. method:: GetMass()

    Get total mass of view.

    :rtype: float

  .. method:: GetCenterOfMass()

    Get the center of mass. For a non-mass-weighted center, see 
    :meth:`GetCenterOfAtoms`.
    
    :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetGeometricCenter()
    
    Get the geometric center (the center of the axis-aligned 
    bounding box).
    
    :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetGeometricStart()
    
    :rtype: :class:`~ost.geom.Vec3`
    
  .. method:: GetCenterOfAtoms()
    
    Calculates the center of atoms. For a mass-weighted center, use 
    :meth:`GetCenterOfMass`.
    
    :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetBondCount()

    Get number of bonds

    :rtype: :class:`int`

  .. method:: GetBondList()

    See :attr:`bonds`
    
    :rtype: :class:`BondHandleList` (list of :class:`BondHandle`)

  .. method:: GetHandle()
  
    See :attr:`handle`

    :rtype: :class:`EntityHandle`

  .. method:: GetGeometricEnd()

    :rtype: :class:`~ost.geom.Vec3`

  .. method:: GetChainList()

    See :attr:`chains`

  .. method:: GetChainCount()

    See :attr:`chain_count`

  .. method:: GetResidueList()

    See :attr:`residues`

  .. method:: GetResidueCount()
    
    See :attr:`residue_count`

  .. method:: GetAtomList()

    See :attr:`atoms`
    
  .. method:: GetAtomCount()
    
    See :attr:`atom_count`

  .. method:: IsValid()
  
    See :attr:`valid`

Chain View
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: ChainView

  A view representation of a :class:`ChainHandle`. Mostly, the same
  functionality is provided as for the handle.

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.

  .. attribute:: handle

    The chain handle this view points to. Also available as :meth:`GetHandle`.

    :type: :class:`ChainHandle`

  .. attribute:: hash_code

    A unique identifier for this chain view. Note, that this is not the same as
    for the chain handle (see :attr:`ChainHandle.hash_code`).
    Also available as :meth:`GetHashCode`.

    :type: int

  .. attribute:: name
  
     The chain name. The name uniquely identifies the chain in the entity. In
     most cases, the chain name is one character. This is a restriction of the
     PDB file format. However, you are free to use longer names as long as you
     don't  want to save them as PDB files.
     
     This property is read-only. To change the name, use an :class:`XCSEditor`.
     
     Also available as :meth:`GetName`
     
     :type: str

  .. attribute:: residues
   
     List of all residues of this chain. The residues are sorted from N- to
     C-terminus. Usually the residue numbers are in ascending order
     (see :attr:`in_sequence`).
   
     This property is read-only.
   
     **Example usage:**
   
     .. code-block:: python
     
       chain=ent.FindChain("A")
       for res in chain.residues:
         print(res.name, res.atom_count)
   
     Also available as :meth:`GetResidueList`. To access a single residue, use
     :meth:`FindResidue`.

     :type: :class:`ResidueViewList` (list of :class:`ResidueView`)
  
  .. attribute:: in_sequence
  
     Whether the residue numbers are in ascending order. For example:
     
     .. code-block:: python
     
       chain=ent.FindChain("A")
       print(chain.residues) # [A.GLY1, A.GLY2, A.GLY4A, A.GLY4B]
       print(chain.in_sequence) # prints true
       
       chain=ent.FindChain("B")
       print(chain.residues) # [B.GLY1, B.GLY4, B.GLY3]
       print(chain.in_sequence) # prints false

    Note that the value of  `in_sequence` is independent from the value of
    :attr:`ChainHandle.in_sequence`.

    :type: bool

  .. attribute:: atoms

     Get list of all atoms of this chain. To access a single atom, use
     :meth:`FindAtom`.
   
     This property is read-only. Also available as :meth:`GetAtomList`

     :type: :class:`AtomViewList` (list of :class:`AtomView`)

  .. attribute:: bounds

    Axis-aligned bounding box of the chain. Read-only
 
    :type: :class:`ost.geom.AlignedCuboid`
    
  .. attribute:: mass
  
    The total mass of this chain in Dalton. Also available as :meth:`GetMass`
  
    :type: float

  .. attribute:: center_of_atoms

    See :attr:`ChainHandle.center_of_atoms`

  .. attribute:: center_of_mass

    See :attr:`ChainHandle.center_of_mass`

  .. attribute:: geometric_center

    See :attr:`ChainHandle.geometric_center`

  .. attribute:: valid

     Validity of view.

     :type: bool

  .. method:: AddAtom(atom_handle[, view_add_flags])

    Add atom to the view. If the residue of the atom is not already part of the 
    view, it will be added. If the atom does not belong to chain, the result is
    undefined.
    
    :param atom_handle: The atom to be added
    :type  atom_handle: :class:`AtomHandle`
    :param view_add_flags: An ORed together combination of :class:`ViewAddFlag`
    :type  view_add_flags: :class:`int` / :class:`ViewAddFlag`
    :rtype: :class:`AtomView`

  .. method:: AddResidue(residue_handle[, view_add_flags])

    Add residue to the view. If the atom does not belong to chain, the result is
    undefined. By default, only the residue, but no atoms are added to the
    view. To change the behavior, pass in a suitable combination of
    :class:`ViewAddFlag`.
    
    :param residue_handle: The residue handle to be added.
    :type  residue_handle: :class:`ResidueHandle`
    :param view_add_flags: An ORed together combination of :class:`ViewAddFlag`
    :type  view_add_flags: :class:`int` / :class:`ViewAddFlag`
    :rtype: :class:`ResidueView`

  .. method:: FindAtom(res_num, atom_name)

    Find atom with the given residue number and atom name. 
    
    :param res_num: The residue number
    :type  res_num: :class:`ResNum`
    :param atom_name: The atom name
    :type  atom_name: str
    :returns: The atom view, or an invalid atom view if no atom with the given 
       parameters is part of the view.
    
    :rtype: :class:`AtomView`

  .. method:: FindResidue(res_num)

    Find residue by residue number.
    
    :param res_num:
    :type  res_num: :class:`ResNum`
    :rtype: :class:`ResidueView`
    :returns: The residue view, or an invalid residue view if no residue with 
       the given residue number is in the view.

  .. method:: ViewForHandle(handle)

    Find residue view or atom view of pointing to the given handle.

    :param handle: handle to search for
    :type  residue: (Residue\|Atom)Handle
    :returns: The view pointing the handle, or an invalid handle if the handle
         is not part of the view
    :rtype: (Residue\|Atom)View

  .. method:: GetCenterOfAtoms()

    See :attr:`center_of_atoms`
    
  .. method:: GetCenterOfMass()

    See :attr:`center_of_mass`
    
  .. method:: GetEntity()

    The parent entity.

  .. method:: GetGeometricCenter()
    
    See :attr:`geometric_center`

  .. method:: GetHandle()

    See :attr:`handle`

  .. method:: GetHashCode()
    
    See :attr:`hash_code`

  .. method:: GetMass()

    See :attr:`mass`

  .. method:: GetName()

    See :attr:`name`

  .. method:: GetResidueByIndex(index)
    
    Returns the residue view at `index`. If index is negative or bigger than 
    the number of residues minus one, an invalid ResidueView is returned.
    
    :param index: The index
    :type  index: int
    :rtype: :class:`ResidueView`

  .. method:: GetResidueList()

    See :attr:`residues`

  .. method:: InSequence()
  
    See :attr:`in_sequence`

  .. method:: IsValid()

    See :attr:`valid`

  .. method:: RemoveResidue(residue)

    Remove residue from chain. If the residue is not part of the view, the 
    chain is left unchanged.
    
    :param residue: The residue view. May be invalid
    :type  residue: :class:`ResidueView`
    :rtype: None

  .. method:: RemoveResidues()
   
    Remove all residues from this chain view

  .. method:: Select(query, flags=0)

    Perform query on chain view. See :meth:`EntityHandle.Select`.
    
    :param query: The query
    :type  query: :class:`Query` / :class:`str`
    :param flags: An ORed together combination of :class:`QueryFlag`
    :type  flags: :class:`int` / :class:`QueryFlag`
    :rtype: :class:`EntityView`

Residue View
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: ResidueView

  A view representation of a :class:`ResidueHandle`. Mostly, the same
  functionality is provided as for the handle.

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.

  .. attribute:: handle

    The residue handle this view points to. Also available as 
    :meth:`GetHandle`.

    :type: :class:`ResidueHandle`

  .. attribute:: hash_code

    A unique identifier for this residue view. Note, that this is not the same as
    for the residue handle (see :attr:`ResidueHandle.hash_code`).
    Also available as :meth:`GetHashCode`.

    :type: int

  .. attribute:: name
                 qualified_name
                 number
                 one_letter_code
                 bounds
                 mass
                 center_of_atoms
                 center_of_mass
                 geometric_center
                 phi_torsion
                 psi_torsion
                 omega_torsion
                 chem_class
                 chem_type
                 sec_structure
                 is_ligand
                 is_protein
                 peptide_linking
                 nucleotide_linking
                 central_atom
                 central_normal
                 valid
                 next
                 prev

    See the respective attributes in :class:`ResidueHandle`.

  .. attribute:: atoms

     Get list of all atoms of this residue included in the view.
     To access a single atom, use
     :meth:`FindAtom`.

     This property is read-only. Also available as :meth:`GetAtomList`

     :type: :class:`AtomViewList` (list of :class:`AtomView`)

  .. attribute:: atom_count

    Number of atoms included in the view. Read-only. See :meth:`GetAtomCount`.

    :type: :class:`int`

  .. attribute:: chain

    The chain this residue belongs to. Read-only. Also available as 
    :meth:`GetChain`
  
    :type: :class:`ChainView`

  .. attribute:: index

    Residue index (starting at 0) within chain view.

  .. method:: RemoveAtom(atom_view)
  
    Remove atom from residue and all associated bonds. If the atom is not part 
    of the view, the residue view is left unchanged.
    
    :param atom_view: The atom to be removed. May be an invalid handle
    :type  atom_view: :class:`AtomView`
    :rtype: None

  .. method:: ViewForHandle(handle)

    Find atom view of pointing to the given handle.

    :param handle: handle to search for
    :type  residue: :class:`AtomHandle`
    :returns: The view pointing the handle, or an invalid handle if the handle
         is not part of the view
    :rtype: :class:`AtomView`

  .. method:: GetHandle()

    See :attr:`handle`

  .. method:: GetHashCode()
    
    See :attr:`hash_code`

  .. method:: GetName
              GetQualifiedName
              GetNumber
              GetOneLetterCode
              GetBounds
              GetMass
              GetCenterOfAtoms
              GetCenterOfMass
              GetGeometricCenter
              GetPhiTorsion
              GetPsiTorsion
              GetOmegaTorsion
              GetChemClass
              GetChemType
              GetSecStructure
              IsLigand
              SetIsLigand
              IsProtein
              IsPeptideLinking
              IsNucleotideLinking
              GetCentralAtom
              SetCentralAtom
              GetCentralNormal
              IsValid
              GetNext
              GetPrev

    See the respective methods in :class:`ResidueHandle`.

  .. method:: GetChain()

    See :attr:`chain`

  .. method:: FindAtom(atom_name)
    
    Find atom by name. Returns the atom, or an invalid atom if no atoms with 
    the given name is part of the view.
    
    :param atom_name:
    :type  atom_name: str
    :rtype: :class:`AtomView`

  .. method:: IsAtomIncluded(atom_handle)

    Returns true if the given atom is part of the view, false if not.
    
    :param atom_handle:
    :type  atom_handle: :class:`AtomHandle`
    :rtype: bool

  .. method:: AddAtom(atom_handle[, flags])

    Add atom to the view.

    :param atom_handle: Atom handle to be added
    :type  atom_handle: :class:`AtomHandle`
    :param flags: An ORed together combination of :class:`ViewAddFlag`
    :type  flags: :class:`int` / :class:`ViewAddFlag`
    :rtype: :class:`AtomView`

  .. method:: GetAtomList()

    See :attr:`atoms`

  .. method:: GetAtomCount()

    See :attr:`atom_count`

  .. method:: GetIndex()

    See :attr:`index`
    
  .. method:: Select(query, flags=0)
   
    Perform selection on residue view. See :meth:`EntityHandle.Select`.
    
    :param query: The query
    :type  query: :class:`Query` / :class:`str`
    :param flags: An ORed together combination of :class:`QueryFlag`
    :type  flags: :class:`int` / :class:`QueryFlag`
    :rtype: :class:`EntityView`

Atom View
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class:: AtomView

  A view representation of an :class:`AtomHandle`. Mostly, the same
  functionality is provided as for the handle.

  .. attribute:: properties

    All the :class:`generic properties <ost.GenericPropContainer>` are available.

  .. attribute:: handle
  
     The underlying :class:`AtomHandle` of the atom view. Also available as
     :meth:`GetHandle`.
     
     :type: :class:`AtomHandle`

  .. attribute:: hash_code

    A unique identifier for this atom view. Note, that this is not the same as
    for the atom handle (see :attr:`AtomHandle.hash_code`).
    Also available as :meth:`GetHashCode`.

    :type: int

  .. method:: GetHandle()

    See :attr:`handle`

  .. method:: GetHashCode()
    
    See :attr:`hash_code`


Functions
--------------------------------------------------------------------------------


Boolean Operators
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: Intersection(view_a, view_b)

  Calculates and returns the intersection of `view_a` and `view_b`. `view_a`
  and `view_b`  must be views of the same entity.

  :param view_a:    first view
  :type view_a:     EntityView
  :param view_b:    second view
  :type view_b:     EntityView

.. function:: Difference(view_a, view_b)

  Calculates and returns the difference between `view_a` and `view_b`. `view_a`
  and `view_b`  must be views of the same entity.The  returned view will
  contain atoms, residues, chains and bonds that are in `view_a` and not in
  `view_b`.

  :param view_a:    first view
  :type view_a:     EntityView
  :param view_b:    second view
  :type view_b:     EntityView 

.. function:: Union(view_a, view_b)

  Calculates and returns the union of `view_a` and `view_b`. `view_a`
  and `view_b`  must be views of the same entity.The returned view will contain
  all atoms, residues, chains and bonds that are either in `view_a`, `view_b`
  or part of both views.

  :param view_a:    first view
  :type view_a:     EntityView
  :param view_b:    second view
  :type view_b:     EntityView


Other Entity-Related Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: CreateViewFromAtoms(atom_list)
              CreateViewFromAtomList(atom_list)

  Returns a view made up of the atoms in *atom_list*. All atoms are required to
  be atoms of the same entity. Duplicate atoms are only added to the view once.
  
  :param atom_list: the atoms
  :type atom_list: :class:`AtomHandleList` or :class:`AtomViewList`
  :raises: :class:`IntegrityError` if atoms of different entities are
           encountered
  
  :returns: :class:`EntityView`

.. function:: CreateEntityFromView(view, include_exlusive_atoms, \
                                   handle=EntityHandle())
 
  This function behaves exactly like :meth:`EntityHandle.Copy`, except that only
  atoms, residues, chains and bonds that are present in the view will be 
  copied.
   
  :param view: is the view to be converted to a handle
  :param include_exlusive_atoms: if true, atoms that are part of an exclusive
       bond (only one bond partner is included in the view) will also be included
       in the new entity handle.
  :param handle: If invalid a new entity will be created. If valid, the atoms, 
       residues, chains, bonds and torsions will be added to handle. This is 
       useful to combine several entities into one.

  :returns: :class:`EntityHandle`

.. function:: InSequence(res, res_next)

  :return: True, if both *res* and *res_next* are :attr:`~ResidueHandle.valid`,
           *res_next* is the residue following *res* (see 
           :attr:`ResidueHandle.next`), both residues are linking (i.e.
           :attr:`~ChemClass.IsPeptideLinking` or
           :attr:`~ChemClass.IsNucleotideLinking`) and are connected by an
           appropriate bond.
  :rtype:  :class:`bool`
  :param res: First residue to check.
  :type res:  :class:`ResidueHandle`
  :param res_next: Second residue to check.
  :type res_next:  :class:`ResidueHandle`

.. function:: BondExists(atom_a, atom_b)

  :return: True, if *atom_a* and *atom_b* are connected by a bond.
  :rtype:  :class:`bool`
  :param atom_a: First atom to check.
  :type atom_a:  :class:`AtomHandle`
  :param atom_b: Second atom to check.
  :type atom_b:  :class:`AtomHandle`

Residue Numbering
--------------------------------------------------------------------------------

.. class:: ResNum(num, ins_code='\\0')

  Number for a residue. The residue number has a numeric part and an (optional)
  insertion-code. You can work with this object as if it was an integer and
  comparison will look first at the numeric part and then the insertion-code.
  All access to existing objects is read-only. Openstructure supports a range
  of (-8388608 to 8388607) for the numeric part. However, the PDB format only 
  supports a range of (-999, 9999). This becomes relevant when a structure is 
  saved in PDB format where an IOException is raised if the PDB range is not 
  respected.

  :param num: Numeric part of residue number.
  :type num:  :class:`int`
  :param ins_code: Alpha-numeric part of residue number (optional insertion
                   code). Only first character kept.
  :type ins_code:  :class:`str`

  .. attribute:: num

    Numeric part of residue number.

    :type: :class:`int`

  .. attribute:: ins_code

    Alpha-numeric part of residue number (insertion code). Single character.

    :type: :class:`str`

  .. method:: GetNum

    :returns: :attr:`num`

  .. method:: GetInsCode

    :returns: :attr:`ins_code`


ChainType
--------------------------------------------------------------------------------

A ChainType fills the :attr:`ChainHandle.type` attribute. Different types are
described in the :class:`ChainType` enum. The type is set with an editor in
:meth:`EditorBase.SetChainType`. Further convenience functions are described
here.

.. class:: ChainType

  The ChainType enum enumerates all types defined by the PDB for the MMCif file
  format. Following values are supported:

    ``CHAINTYPE_POLY``, ``CHAINTYPE_NON_POLY``, ``CHAINTYPE_WATER``,
    ``CHAINTYPE_POLY_PEPTIDE_D``, ``CHAINTYPE_POLY_PEPTIDE_L``,
    ``CHAINTYPE_POLY_DN``, ``CHAINTYPE_POLY_RN``, ``CHAINTYPE_POLY_SAC_D``,
    ``CHAINTYPE_POLY_SAC_L``, ``CHAINTYPE_POLY_DN_RN``,
    ``CHAINTYPE_UNKNOWN``, ``CHAINTYPE_MACROLIDE``,
    ``CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE``, ``CHAINTYPE_POLY_PEPTIDE_DN_RN``,
    ``CHAINTYPE_BRANCHED``, ``CHAINTYPE_OLIGOSACCHARIDE``,
    ``CHAINTYPE_N_CHAINTYPES``

  Where ``CHAINTYPE_N_CHAINTYPES`` holds the number of different types available.

.. function:: ChainTypeFromString(identifier)

   Create a ChainType item for a given string.

   :param identifier: String to request a type for.
   :type identifier: :class:`str`
   :raises: :class:`runtime_error` if **identifier** is unrecognised.

   :returns: :class:`ChainType`

.. function:: StringFromChainType(type)

   Return the String identifier for a given **type**.

   :param type: To be translated
   :type type: :class:`ChainType`
   :raises: :class:`runtime_error` if **type** is unrecognised.

   :returns: :class:`str`


ViewAddFlag
--------------------------------------------------------------------------------

.. class:: ViewAddFlag

  Defines flags controlling behaviour of routines adding handles to views:

  * ``INCLUDE_ATOMS`` - Include all atoms when adding a residue handle to a view
  * ``INCLUDE_RESIDUES`` - Include all residues when adding a chain to a view
  * ``INCLUDE_CHAINS`` - Include all chains when creating a new entity view
  * ``INCLUDE_ALL`` = ``INCLUDE_ATOMS`` | ``INCLUDE_RESIDUES`` |
    ``INCLUDE_CHAINS`` - Convenience flags to include all substructures
  * ``CHECK_DUPLICATES`` - If set, it will be checked that no duplicates are
    created when adding a new handle

  Flags can be ORed to combine them.
  

SecStructure
--------------------------------------------------------------------------------

.. class:: SecStructure(type)

  Defines a secondary structure type following the types defined by DSSP.

  :param type: Type to be set for this object.
  :type type:  :class:`SecStructure.Type` / :class:`str`

  .. method:: IsHelical

    :return: True, if the set type is any type of helix (i.e. ALPHA_HELIX,
             PI_HELIX or THREE_TEN_HELIX)
  
  .. method:: IsExtended

    :return: True, if the set type is any type of beta sheet (i.e. EXTENDED,
             BETA_BRIDGE)

  .. method:: IsCoil

    :return: True, if the set type is any type of coil (i.e. TURN, BEND or COIL)

  .. method:: __str__

    :return: The character corresponding to the set type (see
             :class:`SecStructure.Type`)

.. class:: SecStructure.Type

  Enumerates all popssible secondary structure types distinguished by DSSP.
  Their values with the corresponding character code are listed here:

  .. hlist::
    :columns: 2

    * ``ALPHA_HELIX``     = 'H'
    * ``PI_HELIX``        = 'I'
    * ``THREE_TEN_HELIX`` = 'G'
    * ``EXTENDED``        = 'E'
    * ``BETA_BRIDGE``     = 'B'
    * ``TURN``            = 'T'
    * ``BEND``            = 'S'
    * ``COIL``            = 'C'


ChemClass
--------------------------------------------------------------------------------

.. class:: ChemClass(chem_class)

  The chemical class is used to broadly categorize residues based on their
  chemical properties. For example, peptides belong to some PEPTIDE_LINKING
  class. Possible values as constant variable names and as characters:

  .. hlist::
    :columns: 2

    * ``PEPTIDE_LINKING``   = 'P'
    * ``D_PEPTIDE_LINKING`` = 'D'
    * ``L_PEPTIDE_LINKING`` = 'L'
    * ``RNA_LINKING``       = 'R'
    * ``DNA_LINKING``       = 'S'
    * ``NON_POLYMER``       = 'N'
    * ``L_SACCHARIDE``      = 'X'
    * ``D_SACCHARIDE``      = 'Y'
    * ``SACCHARIDE``        = 'Z'
    * ``WATER``             = 'W'
    * ``UNKNOWN``           = 'U'

  The constants are defined directly within the :mod:`mol` module.
  Python can implicitly convert characters to objects of this type. Note however
  that only the first character of a :class:`str` object is considered!

  :param chem_class: Chemical class to set.
  :type chem_class:  :class:`str`

  .. method:: IsPeptideLinking

    :return: True, if set class is PEPTIDE_LINKING, D_PEPTIDE_LINKING or
             L_PEPTIDE_LINKING

  .. method:: IsNucleotideLinking

    :return: True, if set class is RNA_LINKING or DNA_LINKING
  

ChemType
--------------------------------------------------------------------------------

.. class:: ChemType

  The chemical type of a residue is a classification of all compounds obtained
  from the PDB component dictionary. For example, ions belong to the class IONS,
  amino acids to AMINOACIDS. Possible values as constant variable names and as
  characters:

  .. hlist::
    :columns: 2

    * ``IONS``             = 'I'
    * ``NONCANONICALMOLS`` = 'M'
    * ``SACCHARIDES``      = 'S'
    * ``NUCLEOTIDES``      = 'N'
    * ``AMINOACIDS``       = 'A'
    * ``COENZYMES``        = 'E'
    * ``WATERCOORDIONS``   = 'C'
    * ``DRUGS``            = 'D'
    * ``WATERS``           = 'W'
    * ``UNKNOWN``          = 'U'

  Python can implicitly convert characters to objects of this type. Note however
  that only the first character of a :class:`str` object is considered!

  :param chem_type: Chemical type to set.
  :type chem_type:  :class:`str`

  .. method:: IsIon

    :return: True, if set type is IONS or WATERCOORDIONS

  .. method:: IsNucleotide

    :return: True, if set type is NUCLEOTIDES

  .. method:: IsSaccharide

    :return: True, if set type is SACCHARIDES

  .. method:: IsAminoAcid

    :return: True, if set type is AMINOACIDS

  .. method:: IsCoenzyme

    :return: True, if set type is COENZYMES

  .. method:: IsDrug

    :return: True, if set type is DRUGS

  .. method:: IsNonCanonical

    :return: True, if set type is NONCANONICALMOLS

  .. method:: IsWater

    :return: True, if set type is WATERS

  .. method:: IsKnown

    :return: True, if set type is not UNKNOWN
