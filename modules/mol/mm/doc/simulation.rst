Simulation
================================================================================

.. currentmodule:: ost.mol.mm

The simulation finally connects a :class:`Topology` with an 
:class:`EntityHandle`. While applying minimization or
md tasks, the current positions of the simulation object can be
mapped back to the attached structure at any time.



.. class:: Simulation(entity, settings)

  First constructor that takes an :class:`ost.mol.EntityHandle` as an input
  and automatically constructs an internal :class:`Topology` according
  to the provided **settings**. Be aware, that the **entity** will be altered
  depending on hydrogen/termini definitions in the :class:`Forcefield` attached
  to the **settings**

  :param entity:        Entity from which a default
                        :class:`Topology` will be created according to the
                        **settings**
  :param settings:      Controls the parametrization of
                        the :class:`Simulation` and the buildup of the
                        internal :class:`Topology`

  :type entity:         :class:`EntityHandle` 
  :type settings:       :class:`Settings` 

  :raises:              :class:`RuntimeError` when construction of 
                        :class:`Topology` fails

.. class:: Simulation(topology,entity,settings)

  Second constructor that takes a :class:`Topology`, a consistent
  :class:`ost.mol.EntityHandle` and a :class:`Settings` as input.

  :param topology:      Topology to initialize the :class:`Simulation`
  :param entity:        Entity that is consistent with
                        the provided **topology**
  :param settings:      Controls the parametrization of
                        the :class:`Simulation`

  :type topology:       :class:`Topology`
  :type entity:         :class:`EntityHandle`
  :type settings:       :class:`Settings` 

  :raises:              :class:`RuntimeException` when number of atoms in 
                        **entity** is not equal the number of particles
                        in **topology**


  .. method:: Save(filename, settings)

    Dumps current status of the simulation with its connected **entity** onto 
    disk. All informations from the **entity** except plain chains, residues, 
    atoms and bonds are lost.

    :param filename:    Filename
    :type filename:     :class:`str`

  .. staticmethod:: Load(filename)

    Loads a dumped simulation from disk. You have to make sure, that the provided 
    settings are consistent with those from the saved simulation. Undefined 
    behaviour must be expected otherwise. 

    :param filename: Filename
    :param settings: Settings to reparametrize loaded simulation

    :type filename:     :class:`str`
    :type settings:     :class:`Settings`

  .. staticmethod:: IsPlatformAvailable(settings)

    :return: True, if platform defined in *settings* is available. Otherwise,
             construction of a simulation object will fail with these settings.
    :rtype:  :class:`bool`
    
    :param settings: Controls the parametrization of this class.
                     Only :attr:`Settings.openmm_plugin_directory`,
                     :attr:`Settings.custom_plugin_directory` and
                     :attr:`Settings.platform` are relevant.
    :type settings:  :class:`Settings`

  .. method:: ApplyLBFGS([tolerance=1.0,max_iterations=1000])

    Run minimization using the Limited-memory Broyden–Fletcher–Goldfarb–Shanno 
    (L-BFGS) method.
   
    :param tolerance:   Tolerance describing abort criteria of minimization.
                        The tolerance is interpreted as the root mean square
                        value of all force components.
    :param max_iterations: Maximum number of minimization iterations

    :type tolerance:    :class:`float`
    :type max_iterations: :class:`int`
  

  .. method:: ApplySD([tolerance=100,max_iterations=1000])

    Run minimization using the steepest descend method.
   
    :param tolerance:  Tolerance describing abort criteria of minimization.
                       The tolerance is interpreted as the norm of the maximum 
                       force acting onto the particles.
    :param max_iterations:  Maximum number of minimization iterations

    :type tolerance:    :class:`float`
    :type max_iterations: :class:`int`

    :returns:           :class:`bool` whether tolerance has been reached

  .. method:: Steps(steps)

    Runs the simulation for **steps** number of steps.
   
    :param steps:       Number of steps to perform
    :type steps:        :class:`int`

  .. method:: GetPositions([enforce_periodic_box = false, in_angstrom = true])

    :param enforce_periodic_box: True if periodic box has to be enforced
    :param in_angstrom: True if positions should be in Angstrom, they're in nm
                        otherwise

    :type enforce_periodic_box: :class:`bool`
    :type in_angstrom: :class:`bool`

    :returns: A :class:`ost.geom.Vec3List` of all the positions of the atoms in the Entity.

  .. method:: SetPositions(positions,[,in_angstrom = true])

    :param positions: Positions to be set
    :param in_angstrom: True if **positions** are in Angstrom, they're expected
                        to be in nm otherwise

    :type positions:    :class:`ost.geom.Vec3List`
    :type in_angstrom:  :class:`bool`

  .. method:: UpdatePositions([enforce_periodic_box=false])

    Sets the positions of the atoms in the :class:`EntityHandle` to the current positions of the simulation.

    :param enforce_periodic_box: True if periodic box should be enforced
    :type enforce_periodic_box:  :class:`bool`

  .. method:: GetVelocities()

    :returns: A :class:`~ost.geom.Vec3List` of all the velocities of the atoms in the Entity.
  
  .. method:: GetForces()

    :returns: A :class:`~ost.geom.Vec3List` of all the resulting forces for each atom in the Entity.

  .. method:: SetVelocities(velocities)

    :param velocities:  Velocities to be set
    :type velocities:   :class:`ost.geom.Vec3List`

  .. method:: GetEntity()

    :returns: :class:`EntityHandle` attached to the simulation

  .. method:: GetEnergy()

    :returns: :class:`float` *Potential + Kinetic* energy in kJ/mol of current 
              simulation state

  .. method:: GetPotentialEnergy()

    :returns: :class:`float` *Potential* energy in kJ/mol of current simulation 
              state

  .. method:: GetKineticEnergy()

    :returns: :class:`float` *Kinetic* energy in kJ/mol of current simulation 
              state

  .. method:: GetTopology()

    :returns: :class:`Topology` attached to the simulation

  .. method:: GetPeriodicBoxExtents()

    :returns: :class:`~ost.geom.Vec3` Extents of the periodic box

  .. method:: SetPeriodicBoxExtents(vec)

    :param vec:         Periodic box extents in A
    :type vec:          :class:`ost.geom.Vec3`

  .. method:: Register(observer)

    Adds an observer to the simulation object. This can be used for example to write out the trajectory during the simulation.

    :param observer:    Observer to be attached
    :type observer:     :class:`Observer`

  .. method:: AddPositionConstraint(index)

    Fixes the position of the particle with index given in the argument.
    This requires to reinitialize the internal openmm Context 
    (this is expensive!). Positions, velocities, forces, energies etc. 
    are properly preserved but e.g. states of random number generators etc.
    might be lost.

    :param index:       Particle to be fixed
    :type index:        :class:`int`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of particles

  .. method:: AddPositionConstraints(indices)

    Fixes the position of the atoms with the indices given in in the argument.
    This requires to reinitialize the internal openmm Context 
    (this is expensive!). Positions, velocities, forces, energies etc. 
    are properly preserved but e.g. states of random number generators etc.
    might be lost.

    :param indices:     Particles to be fixed
    :type indices:      :class:`list`

    :raises:            :class:`RuntimeError` when an index in **indices** 
                        exceeds number of particles

  .. method:: ResetPositionConstraints()

    Removes all position constraints.
    This requires to reinitialize the internal openmm Context 
    (this is expensive!). Positions, velocities, forces, energies etc. 
    are properly preserved but e.g. states of random number generators etc.
    might be lost.

  .. method:: ResetHarmonicBond(index, bond_length, force_constant)

    Update of the harmonic bond parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index:       Harmonic Bond to be reset
    :param bond_length: New bond length in nm
    :param force_constant: New force constant in kJ/mol/nm^2

    :type index:        :class:`int`
    :type bond_length:  :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of harmonic bonds

  .. method:: ResetHarmonicAngle()

    Update of the harmonic angle parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index: Harmonic Angle to be reset
    :param angle: New Angle in radians
    :param force_constant: New force constant in kJ/mol/radian^2 

    :type index:        :class:`int`
    :type angle:        :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of harmonic angles

  .. method:: ResetUreyBradleyAngle(index, angle, angle_force_constant, bond_length, bond_force_constant)

    Update of the Urey-Bradley angle parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index: Urey-Bradley angle that has to be reset
    :param angle: New angle in radians
    :param angle_force_constant: New angle force constant kJ/mol/radian^2 
    :param bond_length: New bond length in nm
    :param bond_force_constant: New bond force constant in kJ/mol/nm^2 

    :type index:        :class:`int`
    :type angle:        :class:`float`
    :type angle_force_constant: :class:`float`
    :type bond_length:  :class:`float`
    :type bond_force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of urey-bradley angles

  .. method:: ResetPeriodicDihedral(index, multiplicity, phase, force_constant)

    Update of the periodic dihedral parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index:       Periodic dihedral to be reset
    :param multiplicity: New Multiplicity
    :param phase:       New phase in radians
    :param force_constant: New force constant in kJ/mol/radian^2 

    :type index:        :class:`int`
    :type multiplicity: :class:`int`
    :type phase:        :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of periodic dihedrals


  .. method:: ResetPeriodicImproper(index, multiplicity, phase, force_constant)

    Update of the periodic improper parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index:       Periodic improper to be reset
    :param multiplicity: New multiplicity
    :param phase:       New phase in radians
    :param force_constant: New force constant kJ/mol/radian^2 

    :type index:        :class:`int`
    :type multiplicity: :class:`int`
    :type phase:        :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of periodic impropers

  .. method:: ResetHarmonicImproper(index, phase, force_constant)

    Update of the harmonic improper parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index: Harmonic improper to be reset
    :param phase: New phase in radians
    :param force_constant: New force constant kJ/mol/radian^2 

    :type index:        :class:`int`
    :type phase:        :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of harmonic impropers

  .. method:: ResetLJPair(index,sigma,epsilon)

    Update of the LJ pair parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index:       LJ pair to be reset
    :param sigma:       New sigma in nm
    :param epsilon:     New epsilon in kJ/mol

    :type index:        :class:`int`
    :type sigma:        :class:`float`
    :type epsilon:      :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of LJ pairs


  .. method:: ResetDistanceConstraint(index, constraint_length)

    Update of the distance constraint parameters in the simulation **and**
    in the attached :class:`Topology`. 
    This requires to reinitialize the internal openmm Context 
    (this is expensive!). Positions, velocities, forces, energies etc. 
    are properly preserved but e.g. states of random number generators etc.
    might be lost. 
    
    :param index:       Distance constraint to be reset
    :param constraint_length: New constraint length in nm

    :type index:        :class:`int`
    :type constraint_length: :class:`float` 

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of distance constraints

  .. method:: ResetHarmonicPositionRestraint(index,ref_position, force_constant,x_scale,y_scale,z_scale)

    Update of the harmonic position restraint parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index: Harmonic position restraint to be reset
    :param ref_position: New reference position in nm
    :param force_constant: New force constant in kJ/mol/nm^2 
    :param x_scale:     New x-scale
    :param y_scale:     New y-scale
    :param z_scale:     New z-scale

    :type index:        :class:`int`
    :type ref_position:  :class:`ost.geom.Vec3`
    :type force_constant: :class:`float`
    :type x_scale:      :class:`float`
    :type y_scale:      :class:`float`
    :type z_scale:      :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of harmonic position restraints

  .. method:: ResetHarmonicDistanceRestraint(index, length, force_constant)

    Update of the harmonic distance restraint parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index:       Harmonic distance restraint to be reset
    :param length:      New restraint length in nm
    :param force_constant: New force constant kJ/mol/nm^2 

    :type index:        :class:`int`
    :type length:       :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of harmonic distance restraints

  .. method:: ResetLJ(index,sigma,epsilon)

    Update of the lj parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index:       LJ to be reset
    :param sigma:       New sigma in nm
    :param epsilon:     New epsilon in kJ/mol

    :type index:        :class:`int`
    :type sigma:        :class:`float`
    :type epsilon:        :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of particles

  .. method:: ResetGBSA(index,radius,scaling)

    Update of the GBSA parameters in the simulation **and**
    in the attached :class:`Topology`

    :param index: :     GB to be reset
    :param radius:      New radius in nm
    :param scaling:     New OBC scaling

    :type index:        :class:`int`
    :type radius:       :class:`float`
    :type scaling:      :class:`float`
    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of particles

  .. method:: ResetCharge(index,charge)

    Update of the charge in the simulation **and**
    in the attached :class:`Topology`

    :param index:       Charge to be reset
    :param charge:      New Charge

    :type index:        :class:`int`
    :type charge:       :class:`float`
    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of particles


  .. method:: ResetMass(index,mass)

    Update of the mass in the simulation **and**
    in the attached :class:`Topology`
    This requires to reinitialize the internal openmm Context 
    (this is expensive!). Positions, velocities, forces, energies etc. 
    are properly preserved but e.g. states of random number generators etc.
    might be lost.

    :param index:       Mass to be reset
    :param mass:        New mass

    :type index:        :class:`int`
    :type mass:         :class:`float`

    :raises:            :class:`RuntimeError` when **index** exceeds
                        number of particles



