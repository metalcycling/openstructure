Trajectories
================================================================================

.. currentmodule:: ost.mol


.. function:: CreateCoordGroup(atoms)

  :param atoms: List of atoms. All atoms must be from the same entity.
  :type atoms: :class:`AtomHandleList`
  
  :rtype: :class:`CoordGroupHandle`
  :returns: A coord group with zero frames
  
.. class:: CoordGroupHandle

  A collection of coordinate frames, e.g. an MD trajectory. Create with 
  :func:`CreateCoordGroup`.
  
  
  .. attribute:: entity
  
    The attached entity.
    
    :type: :class:`EntityHandle`
    
  .. attribute:: atoms
  
    The atoms of this coord group. The order of atoms is the same as the 
    positions in the coord frames. All atoms are from the same entity.
    
    :type: :class:`AtomHandleList`
  .. method:: AddFrames(frames)

    Combine two trajectories by appending the frames of the second to the first
    :param frames: a valid coord group
    :type  frames: :class:`CoordGroupHandle`

  .. method:: Capture()

    Record the atom positions of the entity attached to this coord group in a 
    new coordinate frame. Note that the atom positions transformed by the entity 
    transform will be stored. Only available for mutable coord groups.
    
    :see: :attr:`AtomHandle.pos`

  .. method:: CaptureInto(frame_index)

    Same as :meth:`Capture`, but doesn't create a new frame and stores the 
    coordinates directly into frame with index *frame_index*. Only available for 
    mutable coord groups.
    
    :param frame_index: index of the frame
    :type  frame_index: int

  .. method:: CopyFrame(frame_index)

    Copies the coordinates of frame with index *frame_index* to the attached 
    entity.
    
    :param frame_index: index of the frame
    :type  frame_index: int

  .. method:: Filter(view)

    Returns a new trajectory containing only coordinates of the atoms in view. 
    Useful to remove water and other solvent molecules from a trajectory to 
    save memory.
    
    :param view: a valid entity view
    :type  view: :class:`EntityView`
    :rtype: :class:`CoordGroupHandle`

  .. method:: GetAtomCount()

    Returns the number of atoms in the coord group
    
    :rtype: int


  .. method:: GetAtomList()

    Returns the atoms of the coord group in the same order they appear in the 
    coordinate frames.
    
    :rtype: :class:`AtomHandleList`

  .. method:: GetAtomPos(frame_index, atom)

    Get position of *atom* in frame with index *frame_index*.
    

    :param frame_index: frame index
    :type  frame_index: int
    :param atom: A valid atom
    :type  atom: :class:`AtomHandle`
    :rtype: :class:`Vec3`

  .. method:: GetEntity()

    Returns the attached entity
    
    :rtype: :class:`EntityHandle`

  .. method:: GetFrameCount()

    Returns the number of frames of this coord group
    
    :rtype: int

  .. method:: IsValid()

    Whether this coord group is valid
    
    :rtype: bool

  .. method:: SetAtomPos(frame_index, atom, pos)

    Set position of *atom* in frame with index *frame_index* to *pos*
    
    :param frame_index: index of the frame
    :type  frame_index: int
    :param atom: a valid atom
    :type  atom: :class:`AtomHandle`
    :param pos: new position of the atom
    :type  pos: :class:`Vec3`

  .. method:: SetFramePositions(frame_index, positions)

    Set the frame positions of frame with index *frame_index*. Order and count 
    of positions must match :attr:`atoms`.
    
    :param frame_index: index of frame
    :type  frame_index: int
    :param positions: list of positions
    :type  positions: :class:`~ost.geom.Vec3List`


.. class:: CoordFrame

  A single frame of coordinates in a :class:`CoordGroupHandle`.
  
  .. method:: GetAngle(atom1, atom2, atom3)

    :param atom1: first atom
    :type  atom1: :class:`AtomHandle`
    :param atom2: second (central) atom
    :type  atom2: :class:`AtomHandle`
    :param atom3: third atom
    :type  atom3: :class:`AtomHandle`
    
    :returns: the angle in radians between the 3 atoms
    :rtype: float

  .. method:: GetAtomPos(atom)

    Returns the position of the atom in the coord frame
    
    :param atom: A valid atom handle
    :type  atom: :class:`AtomHandle`
    

    
    :rtype: :class:`Vec3`

  .. method:: GetCenterOfMassPos(view)

    
    :param view: A valid entity view
    :type  view: :class:`EntityView`
    :rtype: :class:`Vec3`

  .. method:: GetDihedralAngle(atom1, atom2, atom3, atom4)

    Get dihedral angle of the four atoms.
    
    :param atom1: First atom. Must be valid
    :type  atom1: :class:`AtomHandle`
    :param atom2: Second atom. Must be valid
    :type  atom2: :class:`AtomHandle`
    :param atom3: Third atom. Must be valid
    :type  atom3: :class:`AtomHandle`
    :param atom3: Fourth atom. Must be valid
    :type  atom3: :class:`AtomHandle`
    
    :rtype: float

  .. method:: GetDistanceBetwAtoms(atom1, atom2)

    Get distance in (Angstroem) between *atom1* and *atom2* in coord frame.
    
    :param atom1: First atom. Must be valid
    :type  atom1: :class:`AtomHandle`
    :param atom2: Second atom. Must be valid
    :type  atom2: :class:`AtomHandle`
    :rtype: float

  .. method:: GetDistanceBetwCenterOfMass(view1, view2)

    Get distance between center of mass of the first selection and the second.
    
    :param view1: First view. Must be valid
    :type  view1: :class:`EntityView`
    :param view2: Second view. Must be valid
    :type  view2: :class:`EntityView`
    :rtype: float

  .. method:: GetRMSD(view1, view2)

    Get RMSD between two views in the coord frame. The coordinates of the views 
    are taken as is without superposition. The two views must have the same 
    number of atoms. Atoms are matches as they appear in 
    :attr:`EntityView.atoms`.
    
    :param view1: First view. Must be valid
    :type  view1: :class:`EntityView`
    :param view2: Second view. Must be valid
    :type  view2: :class:`EntityView`
    :rtype: float

