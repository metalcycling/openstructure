Topology
================================================================================

.. currentmodule:: ost.mol.mm

The :class:`Topology` object is an abstract representation of a protein 
structure or any combination of particles that interact with each other in 
some way. Being initialized with a list of masses representing the particles, 
the topology allows to define arbitrary forces between these particles. The 
:class:`Topology` gets internally translated to a molecular mechanics system 
as soon as a :class:`Simulation` is set up. To make setting up a 
:class:`Topology` more convenient, the :class:`TopologyCreator` class allows 
to create a :class:`Topology` based on an mm :class:`Settings` object, 
that should suffice for most needs.

Automatically create topologies
--------------------------------------------------------------------------------

.. class:: TopologyCreator

  .. method:: Create(entity, settings)

    Creates a topology given an :class:`EntityHandle` and a :class:`Settings`
    object. Note, that the :class:`EntityHandle` is passed by reference
    and gets altered by rebuilding bonds, adding hydrogens and modification of 
    terminal residues based on the information given by the 
    :class:`Settings`/:class:`Forcefield`.
    The Topology Creation process consists of:

    #. Rename to the Gromacs standard
    
    #. Kill all bonds
    
    #. Build difulfide bonds if flag is set in **settings**
    
    #. Rename residues/atoms to the standard defined in the forcefield
    
    #. Extract the building blocks from the forcefield for every residue
    
    #. Apply BlockModifiers on termini and reconstruct hydrogens with HydrogenConstructors
    
    #. Reconnect based on the connectivity definitions in the building blocks
    
    #. Check, whether the residues match with the building blocks
    
    #. Extract bonds/angles etc. based on the flags in the **settings**

    #. Impropers and CMaps get only extracted when they are defined in the :class:`BuildingBlock`
    
    #. Add Parametrization from :class:`Forcefield` to the :class:`Topology`. Note, that you can override the :class:`Forcefield` parametrization if you add this particular :class:`Interaction` in a parameterized form in the :class:`BuildingBlock`
    
    #. Rename entity back to PDB naming if flag is set accordingly in **settings**

    :param entity:      Entity for which to create :class:`Topology`
    :param settings:    Settings to control :class:`Topology` buildup

    :type entity:       :class:`ost.mol.EntityHandle`
    :type settings:     :class:`Settings`


    :raises:            :class:`RuntimeError` if settings has no 
                        :class:`Forcefield` assigned

    :raises:            :class:`RuntimeError` if there is a residue with no
                        appropriate :class:`BuildingBlock` in the 
                        :class:`Forcefield`

    :raises:            :class:`RuntimeError` if there is a 
                        :class:`ResidueHandle` not matching the
                        :class:`BuildingBlock` definition in the
                        :class:`Forcefield`

    :raises:            :class:`RuntimeError` if there is an interaction
                        that has no parameters in the forcefield. This 
                        :class:`RuntimeError` can be supressed by setting the 
                        appropriate flag in the :class:`Settings` object.


The Topology Class
--------------------------------------------------------------------------------

.. class:: Topology(masses)

  :param masses:        :class:`list` of the particles'masses

  .. method:: Save(filename)

    Dump Toplogy into binary file

    :param filename:    Name of file
    :type filename:     :class:`str`


  .. method:: Load(filename)

    Load binary topology from disk

    :param filename:    Name of file
    :type filename:     :class:`str`

    :returns:           loaded :class:`Topology`

    :raises:            :class:`RuntimeError` when file can't be found


  .. method:: Merge(other)

    Merges in another :class:`Topology` and adds all their particles and interactions.
    Assuming the actual :class:`Topology` containing n particles, the indexing of the
    other :class:`Topology` gets altered by adding n everywhere.

    :param other:       Topology to merge in
    :type other:        :class:`Topology`

    :raises:            :class:`RuntimeError` if fudge parameters are not consistent

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains charges but the **other** not

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains sigmas but the **other** not

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains epsilons but the **other** not

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains GBSA Radii but the **other** not

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains OBC scalings but the **other** not

  .. method:: Merge(ent, other, other_ent)

    Merges in another :class:`Topology` and adds all their particles and 
    interactions. Assuming the actual :class:`Topology` contains n particles, 
    the indexing of the other :class:`Topology` gets altered by adding n 
    everywhere. This version of the Merge function assumes the current topology 
    describes the :class:`EntityHandle` ent and the other :class:`Topology`
    describes the :class:`EntityHandle` other_ent. During merging, all
    chains, residues, atoms and bonds from **other_ent** are pulled into **ent** in
    a way, that the current :class:`Topology` is consistent with the
    :class:`EntityHandle` **ent**.

    :param ent:         Entity consistent with current topology
    :param other:       Topology to merge in
    :param other_ent:   Entity consistent with **other**

    :type ent:          :class:`ost.mol.EntityHandle`
    :type other:        :class:`Topology`
    :type other_ent:    :class:`ost.mol.EntityHandle`

    :raises:            :class:`RuntimeError` if fudge parameters are not consistent

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains charges but the other not

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains sigmas but the other not

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains epsilons but the other not

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains GBSA Radii but the other not

    :raises:            :class:`RuntimeError` if current :class:`Topology` 
                        contains OBC scalings but the other not

    :raises:            :class:`RuntimeError` if the number of atoms in one
                        of the entities is not consistent with their related
                        topologies.

    :raises:            :class:`RuntimeError` other_ent contains a chain with
                        a name already present in **ent**.


  .. method:: AddHarmonicBond(index_one, index_two, bond_length, force_constant)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param bond_length: Bond length in nm
    :param force_constant: Force constant in kJ/mol/nm^2

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type bond_length:  :class:`float`
    :type force_constant: :class:`float`


    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles



  .. method:: AddHarmonicAngle(index_one, index_two, index_three, angle, \
                               force_constant)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param angle:       Angle in radians
    :param force_constant: Force constant in kJ/mol/radian^2 

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type angle:        :class:`float`
    :type force_constant: :class:`float`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: AddUreyBradleyAngle(index_one, index_two, index_three, \
                                  angle, angle_force_constant, \
                                  bond_length, bond_force_constant)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param angle:       Angle in radians
    :param angle_force_constant: Angle force constant kJ/mol/radian^2 
    :param bond_length: Bond length in nm
    :param bond_force_constant: Bond force constant kJ/mol/nm^2 

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type angle:        :class:`float`
    :type angle_force_constant: :class:`float`
    :type bond:         :class:`float`
    :type bond_force_constant: :class:`float`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: AddPeriodicDihedral(index_one, index_two, index_three, \
                                  index_four, multiplicity, phase, \
                                  force_constant)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param index_four:  Index of particle 4
    :param multiplicity: Multiplicity
    :param phase:       Phase in radians
    :param force_constant: Force constant in kJ/mol/radian^2 

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type index_four:   :class:`int`
    :type multiplicity: :class:`int`
    :type phase:        :class:`float`
    :type force_constant: :class:`float`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: AddPeriodicImproper(index_one, index_two, index_three, \
                                  index_four, multiplicity, phase, \
                                  force_constant)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param index_four:  Index of particle 4
    :param multiplicity: Multiplicity
    :param phase:       Phase in radians
    :param force_constant: Force constant in kJ/mol/radian^2 

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type index_four:   :class:`int`
    :type multiplicity: :class:`int`
    :type phase:        :class:`float`
    :type force_constant: :class:`float`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: AddHarmonicImproper(index_one, index_two, index_three, \
                                  index_four, angle, force_constant)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param index_four:  Index of particle 4
    :param angle:       Angle in radians
    :param force_constant: Force constant kJ/mol/radian^2 

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type index_four:   :class:`int`
    :type angle:        :class:`float`
    :type force_constant: :class:`float`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: AddCMap(index_one, index_two, index_three, index_four, \
                      index_five, dimension, values)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param index_four:  Index of particle 4
    :param index_five:  Index of particle 5
    :param dimension:   Dimension of CMap
    :param values:      Map of size 
                        **dimension** * **dimension** describing the correction
                        map in kJ/mol. The values have to be ordered like
                        (phi1,psi1),(phi1,psi2)...(phi1,psin),
                        (phi2,psi1),(phi2,psi2)...(phi2,psin),
                        ...(phin,psin).
                        the range for phi and psi are [-pi,pi[

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type index_four:   :class:`int`
    :type index_five:   :class:`int`
    :type dimension:    :class:`int`
    :type values:       :class:`list`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when size of values is not 
                        **dimension** * **dimension** or index exceeds number of 
                        particles.


  .. method:: AddLJPair(index_one, index_two, sigma, epsilon)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param sigma:       Sigma in nm
    :param epsilon:     Epsilon in kJ/mol

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type sigma:        :class:`float`
    :type epsiolon:     :class:`float`


    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles or when particular LJPair is already 
                         parameterized


  .. method:: AddDistanceConstraint(index_one, index_two, distance)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param distance:    Constraint distance in nm

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type distance:     :class:`float`    

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles or when particular distance constraint
                         is already parameterized


  .. method:: AddExclusion(index_one, index_two)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles or when particular exclusion is already set


  .. method:: AddPositionConstraint(index)

    :param index:       Index of particle to be constrained

    :type index:        :class:`int`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: AddHarmonicPositionRestraint(index, reference_pos, \
                                           force_constant, x_scale=1.0, \
                                           y_scale=1.0, z_scale=1.0)

    Adds a harmonic position restraint for a particle in the form:
    force_constant*(x_scale*(x-x0)^2+y_scale*(y-y0)^2+z_scale*(z-z0)^2)

    :param index:       Index of particle to be constrained
    :param reference_pos: Reference position in A
    :param force_constant: Force constant in kJ/mol/nm^2 
    :param x_scale:     X-scale
    :param y_scale:     Y-scale
    :param z_scale:     Z-scale

    :type index:        :class:`int`
    :type reference_pos: :class:`ost.geom.Vec3`
    :type force_constant: :class:`float`
    :type x_scale:      :class:`float`
    :type y_scale:      :class:`float`
    :type z_scale:      :class:`float`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: AddHarmonicDistanceRestraint(index_one, index_two, length, \
                                           force_constant)

    Adds a force between two particles in the form:
    force_constant * (b-b0)^2

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param length:      Restraint length in nm
    :param force_constant: Force constant in kJ/mol/nm^2 

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type length:       :class:`float`
    :type force_constant: :class:`float`

    :returns:           :class:`int` index of added interaction

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: ResetPositionConstraints()

    Removes all position constraints


  .. method:: ResetExclusions()

    Removes all set exclusions


  .. method:: SetSigmas(sigmas)

    :param sigmas:      Sigma parameters for all particles in nm

    :type sigmas:       :class:`list` (:class:`float`)
   
    :raises:            :class:`RuntimeError` when list size is inconsistent
                        with number of particles


  .. method:: SetSigma(index, sigma)

    :param index:       Index of particle in nm
    :param sigma:        Sigma

    :type index:        :class:`int`
    :type sigma:        :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: SetEpsilons(epsilons)

    :param epsilons:    Epsilon parameters for all particles in kJ/mol

    :type epsilons:     :class:`list` (:class:`float`)
   
    :raises:            :class:`RuntimeError` when list size is inconsistent
                        with number of particles


  .. method:: SetEpsilon(index, epsilon)

    :param index:       Index of particle 
    :param epsilon:     Epsilon in kJ/mol

    :type index:        :class:`int`      
    :type epsilon:      :class:`float`      

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: SetGBSARadii(radii)

    :param radii:       GBSA radii for all particles in nm

    :type radii:        :class:`list` (:class:`float`)
   
    :raises:            :class:`RuntimeError` when list size is inconsistent
                        with number of particles


  .. method:: SetGBSARadius(index, radius)

    :param index:       Index of particle
    :param radius:      GBSA radius in nm

    :type index:        :class:`int`
    :type radius:       :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: SetOBCScalings(scalings)

    :param scalings:    OBC scaling factors for all particles
   
    :type scalings:     :class:`list` (:class:`float`)

    :raises:            :class:`RuntimeError` when list size is inconsistent
                        with number of particles


  .. method:: SetOBCScaling(index, scaling)

    :param index:       Index of particle
    :param scaling:     OBC scaling factor

    :type index:        Index of particle
    :type scaling:      :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: SetCharges(charges)

    :param charges:     Charges for all particles

    :type charges:      :class:`list` (:class:`float`)
   
    :raises:            :class:`RuntimeError` when list size is inconsistent
                        with number of particles


  .. method:: SetCharge(index, charge)

    :param index:       Index of particle
    :param charge:      Charge

    :type index:        :class:`int`
    :type charge:       :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles

  .. method:: SetMasses(masses)

    :param masses:      Masses for all particles
    :type masses:       :class:`list` (:class:`float`)
   
    :raises:            :class:`RuntimeError` when list size is inconsistent
                        with number of particles


  .. method:: SetMass(index, mass)

    :param index:       Index of particle
    :param mass:        Mass

    :type index:        :class:`int`
    :type mass:         :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         particles


  .. method:: SetFudgeQQ(fudge)

    :param fudge:       Dampening factor for electrostatic 1,4 interactions
    :type fudge:        :class:`float`

  .. method:: SetFudgeLJ(fudge)

    :param fudge:       Dampening factor for LJ 1,4 interactions       
    :type fudge:        :class:`float`


  .. method:: GetHarmonicBondParameters(index)

    :param index:       Index of harmonic bond
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic bonds


  .. method:: GetHarmonicAngleParameters(index)

    :param index:       Index of harmonic angle
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic angles

  .. method:: GetUreyBradleyAngleParameters(index)

    :param index:       Index of Urey-Bradley angle
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set Urey-Bradley angles

  .. method:: GetPeriodicDiheralParameters(index)

    :param index:       Index of periodic dihedral
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set periodic dihedrals

  .. method:: GetPeriodicImproperParameters(index)

    :param index:       Index of periodic improper
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set periodic impropers

  .. method:: GetHarmonicImproperParameters(index)

    :param index:       Index of harmonic improper
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic impropers


  .. method:: GetCMapParameters(index)

    :param index:       Index of CMap
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set cmaps


  .. method:: GetLJPairParameters(index)

    :param index:       Index of LJ pair
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set lj pairs


  .. method:: GetDistanceConstraintParameters(index)

    :param index:       Index of distance constraint
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set distance constraints


  .. method:: GetHarmonicPositionRestraintParameters(index)

    :param index:       Index of harmonic position restraint
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic position restraints


  .. method:: GetHarmonicDistanceRestraintParameters(index)

    :param index:       Index of harmonic distance restraint
    :type index:        :class:`int`

    :returns:           :class:`tuple` containing the indices of the involved
                        particles and the set parameters

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic distance restraints


  .. method:: SetHarmonicBondParameters(index, bond_length, force_constant)

    Replaces existing parameters for particular interaction

    :param index:       Index of harmonic bond
    :param bond_length: New bond length in nm
    :param force_constant: New force constant in kJ/mol/nm^2

    :type index:        :class:`int`
    :type bond_length:  :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic bonds


  .. method:: SetHarmonicAngleParameters(index, angle, force_constant)

    Replaces existing parameters for particular interaction

    :param index:       Index of harmonic angle
    :param angle:       New angle in radians
    :param force_constant: New force constant kJ/mol/radian^2

    :type index:        :class:`int`
    :type angle:        :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic angles

  .. method:: SetUreyBradleyAngleParameters(index, angle, angle_force_constant,
                                            bond_length, bond_force_constant)

    Replaces existing parameters for particular interaction

    :param index:       Index of Urey-Bradley angle
    :param angle:       New angle in radians 
    :param angle_force_constant: New angle force constant in kJ/mol/radian^2
    :param bond_length: New bond length in nm
    :param bond_force_constant: New bond force constant in kJ/mol/nm^2

    :type index:        :class:`int`
    :type angle:        :class:`float`
    :type angle_force_constant: :class:`float`
    :type bond_length:  :class:`float`
    :type bond_force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set Urey-Bradley angles

  .. method:: SetPeriodicDiheralParameters(index, multiplicity, phase, \
                                           force_constant)

    Replaces existing parameters for a particular interaction

    :param index:       Index of periodic dihedral
    :param multiplicity: New periodicity
    :param phase:       New phase in radians
    :param force_constant: New force constant in kJ/mol/radian^2

    :type index:        :class:`int`
    :type multiplicity: :class:`int`
    :type phase:        :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set periodic dihedrals

  .. method:: SetPeriodicImproperParameters(index, multiplicity, phase, \
                                            force_constant)

    Replaces existing parameters for a particular interaction

    :param index:       Index of periodic improper
    :param multiplicity: New multiplicity
    :param phase:       New phase in radians
    :param force_constant: New force constant kJ/mol/radian^2

    :type  index:       :class:`int`
    :type  multiplicity: :class:`int`
    :type  phase:       :class:`float`
    :type  force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set periodic impropers

  .. method:: SetHarmonicImproperParameters(index, angle, force_constant)

    Replaces existing parameters for a particular interaction

    :param index:       Index of harmonic improper
    :param angle:       New angle in radians
    :param force_constant: New force constant kJ/mol/radian^2

    :type index:        :class:`int`
    :type angle:        :class:`float`
    :type force_constant: :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic impropers


  .. method:: SetCMapParameters(index, dimension, map)

    Replaces existing parameters for a particular interaction

    :param index:       Index of CMap
    :param dimension:   New dimension
    :param map:         New CMap of format described in the AddCMap function.

    :type index:        :class:`int`
    :type dimension:    :class:`int`
    :type map:          :class:`list` 

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set cmaps or size of map is not dimension*dimension


  .. method:: SetLJPairParameters(index, sigma, epsilon)

    Replaces existing parameters for a particular interaction

    :param index:       Index of LJ pair
    :param sigma:       New sigma in nm
    :param epsilon:     New epsilon in kJ/mol

    :type index:        :class:`int`
    :type sigma:        :class:`float`
    :type epsilon:      :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set lj pairs


  .. method:: SetDistanceConstraintParameters(index, distance)

    Replaces existing parameters for a particular interaction

    :param index:       Index of distance constraint
    :param distance:    New constraint distance in nm

    :type index:        :class:`int`
    :type distance:     :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set distance constraints


  .. method:: SetHarmonicPositionRestraintParameters(index, ref_position, k,\
                                                     x_scale=1.0, y_scale=1.0,\
                                                     z_scale=1.0)

    Replaces existing parameters for a particular interaction

    :param index:       Index of harmonic position constraint
    :param ref_position: New reference position in nm
    :param k:           New force constant in kJ/mol/nm^2
    :param x_scale:     New x-scale
    :param y_scale:     New y-scale 
    :param z_scale:     New z-scale

    :type index:        :class:`int`
    :type ref_position: :class:`ost.geom.Vec3`
    :type k:            :class:`float`
    :type x_scale:      :class:`float`
    :type y_scale:      :class:`float`
    :type z_scale:      :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic position restraints


  .. method:: SetHarmonicDistanceRestraintParameters(index, length, k)

    Replaces existing parameters for a particular interaction

    :param index:       Index of harmonic distance restraint
    :param length:      New restraint distance in nm
    :param k:           New force constant in kJ/mol/nm^2

    :type index:       :class:`int`
    :type length:      :class:`float`
    :type k:           :class:`float`

    :raises:            :class:`RuntimeError` when index exceeds number of 
                         set harmonic distance restraints


  .. method:: GetSigmas()

    :returns:           :class:`list` of :class:`float` with previously set 
                        sigmas


  .. method:: GetEpsilons()

    :returns:           :class:`list` of :class:`float` with previously set 
                        epsilons


  .. method:: GetGBSARadii()

    :returns:           :class:`list` of :class:`float` with previously set 
                        GBSA radii


  .. method:: GetOBCScalings()

    :returns:           :class:`list` of :class:`float` with previously set 
                        OBC scalings


  .. method:: GetCharges()

    :returns:           :class:`list` of :class:`float` with previously set 
                        charges


  .. method:: GetMasses()

    :returns:           :class:`list` of :class:`float` with previously set 
                        masses


  .. method:: GetSigma(particle_index)

    :param particle_index: Index of particle
    :type particle_index: :class:`int`

    :returns:           :class:`float`

    :raises:            :class:`RuntimeError` when particle_index exceeds
                        number of particles or sigmas have not been set


  .. method:: GetEpsilon(particle_index)

    :param particle_index: Index of particle
    :type particle_index: :class:`int`

    :returns:           :class:`float`

    :raises:            :class:`RuntimeError` when particle_index exceeds
                        number of particles or epsilons have not been set



  .. method:: GetGBSARadius(particle_index)

    :param particle_index: Index of particle
    :type particle_index: :class:`int`

    :returns:          :class:`float`

    :raises:            :class:`RuntimeError` when particle_index exceeds
                        number of particles or gbsa_radii have not been set


  .. method:: GetOBCScaling(particle_index)

    :param particle_index: Index of particle
    :type particle_index: :class:`int`

    :returns:           :class:`float`

    :raises:            :class:`RuntimeError` when particle_index exceeds
                        number of particles or obc scalings have not been set


  .. method:: GetCharge(particle_index)

    :param particle_index: Index of particle
    :type particle_index: :class:`float`

    :returns:           :class:`float`

    :raises:            :class:`RuntimeError` when particle_index exceeds
                        number of particles or charges have not been set


  .. method:: GetMass(particle_index)

    :param particle_index: Index of particle
    :type particle_index: :class:`int`

    :returns:           :class:`float`

    :raises:            :class:`RuntimeError` when particle_index exceeds
                        number of particles


  .. method:: GetHarmonicBondIndices(index_one, index_two)

    :param index_one:   Index of particle 1
    :param index_one:   Index of particle 2

    :type index_one:    :class:`int`
    :type index_one:    :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of 
                        all harmonic bonds acting on given particle pair

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetHarmonicAngleIndices(index_one, index_two, index_three)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all harmonic
                        angles acting on given particles

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetUreyBradleyAngleIndices(index_one, index_two, index_three)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all 
                        Urey-Bradley angles acting on given particles

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetPeriodicDihedralIndices(index_one, index_two, index_three, \
                                         index_four)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param index_four: Index of particle 4

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type index_four:  :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all periodic
                        dihedrals acting on given particles

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetPeriodicImpropersIndices(index_one, index_two, index_three, \
                                          index_four)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param index_four: Index of particle 4

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type index_four:  :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all periodic
                        impropers acting on given particles

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetHarmonicImpropersIndices(index_one, index_two, index_three, \
                                          index_four)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param index_four: Index of particle 4

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type index_four:  :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all harmonic
                        impropers acting on given particles

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetCMapIndices(index_one, index_two, index_three, \
                             index_four, index_five)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2
    :param index_three: Index of particle 3
    :param index_four: Index of particle 4
    :param index_five: Index of particle 5

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`
    :type index_three:  :class:`int`
    :type index_four:  :class:`int`
    :type index_five:  :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all cmaps
                        acting on given particles

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetLJPairIndex(index_one, index_two)

    :param index_one:   Index of particle 1
    :param index_two:   Index of particle 2

    :type index_one:    :class:`int`
    :type index_two:    :class:`int`

    :returns:           :class:`int` Index of LJPair acting on particles,
                        -1 if there is none

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetDistanceConstraintIndex(index_one, index_two)

    :param index_one:   Index of particle 1
    :param index_one:   Index of particle 2

    :type index_one:    :class:`int`
    :type index_one:    :class:`int`

    :returns:           :class:`int` - Index of distance constraint acting on 
                        particles, -1 if there is none

    :raises:            :class:`RuntimeError` when an index exceeds
                        number of particles in topology


  .. method:: GetHarmonicBondIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all harmonic
                        bonds this particular particle is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetHarmonicAngleIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all harmonic
                        angles this particular particle is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetUreyBradleyAngleIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all 
                        Urey-Bradley angles this particular particle 
                        is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetPeriodicDihedralIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all periodic
                        dihedrals this particular particle is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetPeriodicImproperIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all periodic
                        impropers this particular particle is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology

  .. method:: GetHarmonicImproperIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all harmonic
                        impropers this particular particle is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetCMapIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all cmaps
                        this particular particle is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetLJPairIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all LJ pairs
                        this particular particle is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetDistanceConstraintIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all distance
                        constraints this particular particle is involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetHarmonicPositionRestraintIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all harmonic
                        position restraints this particular particle is 
                        involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology
  

  .. method:: GetHarmonicDistanceRestraintIndices(index)

    :param index:       Index of particle
    :type index:        :class:`int`

    :returns:           :class:`list` of :class:`int` - Indices of all harmonic
                        distance restraints this particular particle is 
                        involved in

    :raises:            :class:`RuntimeError` when index exceeds
                        number of particles in topology


  .. method:: GetNumParticles()

    :returns:           :class:`int`


  .. method:: GetNumHarmonicBonds()

    :returns:           :class:`int`


  .. method:: GetNumHarmonicAngles()

    :returns:           :class:`int`


  .. method:: GetNumUreyBradleyAngles()

    :returns:           :class:`int`


  .. method:: GetNumPeriodicDihedrals()

    :returns:           :class:`int`


  .. method:: GetNumPeriodicImpropers()

    :returns:           :class:`int`


  .. method:: GetNumHarmonicImpropers()

    :returns:           :class:`int`


  .. method:: GetNumCMaps()

    :returns:           :class:`int`


  .. method:: GetNumLJPairs()

    :returns:           :class:`int`


  .. method:: GetNumDistanceConstraints()

    :returns:           :class:`int`


  .. method:: GetNumPositionConstraints()

    :returns:           :class:`int`


  .. method:: GetNumHarmonicPositionRestraints()

    :returns:           :class:`int`


  .. method:: GetNumHarmonicDistanceRestraints()

    :returns:           :class:`int`


  .. method:: GetNumExclusions()

    :returns:           :class:`int`
