The Molecular Mechanics Settings
================================================================================

.. currentmodule:: ost.mol.mm

The :class:`Settings` define all parameters to control the buildup of a 
:class:`Topology` in the :class:`TopologyCreator` and the final setup
of the :class:`Simulation` object.

.. class:: Settings

  .. attribute::        add_bonds  

    Flag, whether bonds should be parametrized
    in :class:`TopologyCreator`. Default: True

  .. attribute::        add_angles 

    Flag, whether angles should be parametrized
    in :class:`TopologyCreator`. Default: True

  .. attribute::        add_dihedrals 

    Flag, whether dihedrals should be parametrized
    in :class:`TopologyCreator`. Default: True

  .. attribute::        add_impropers 

    Flag, whether impropers should be parametrized
    in :class:`TopologyCreator`. Default: True

  .. attribute::        add_cmaps 

    Flag, whether cmaps should be parametrized
    in :class:`TopologyCreator`. Default: True

  .. attribute::        add_exclusions 

    Flag, whether exclusions should be parametrized
    in :class:`TopologyCreator`. Default: True

  .. attribute::        add_nonbonded

    Flag, whether nonbonded interactions should be parametrized
    in :class:`TopologyCreator`. Default: True

  .. attribute::        add_gbsa 

    Flag, whether GBSA interactions should be parametrized
    in :class:`TopologyCreator`. Default: False

  .. attribute::        constrain_hbonds

    Flag, whether bonds involving hydrogens should be constrained
    in :class:`TopologyCreator`. Default: False

  .. attribute::        constrain_bonds

    Flag, whether all bonds should be constrained in 
    :class:`TopologyCreator`. Default: False

  .. attribute::        rigid_water

    Flag, whether water molecules should be made rigid in
    :class:`TopologyCreator`. This is achieved by adding
    a distance constraint on the O-H bonds and an additional
    one between the two hydrogens. Default: False

  .. attribute::        strict_interactions

    Makes sure that an error is thrown when a particular interaction cannot
    be parametrized with the given forcefield in the :class:`TopologyCreator`. 
    By setting it to False, these cases just get ignored. Default: True

  .. attribute::        ideal_bond_length_constraints

    Flag, whether the ideal bond length from the :class:`Forcefield` should be
    taken for distance constraints in the :class:`TopologyCreator`. The actual
    distances from the :class:`~ost.mol.EntityHandle` get taken otherwise.
    Default: True

  .. attribute::        fix_heavy_atoms

    Flag, whether all heavy atom positions (non hydrogens) should be fixed in
    space in the :class:`TopologyCreator`. Default: False

  .. attribute::        kill_electrostatics

    Flag, whether all charges should just be set to zero in the
    :class:`TopologyCreator`. Default: False

  .. attribute::        generate_disulfid_bonds

    Flag, whether disulfid bonds based on a simple geometric criterion 
    (SG-SG dist < 2.5A) should be built in :class:`TopologyCreator`. 
    Default: True

  .. attribute::        nonbonded_method

    Nonbonded method set up at the creation of the :class:`Simulation`.
    Must be one of mm.NoCutoff, mm.CutoffNonPeriodic, mm.CutoffPeriodic
    mm.Ewald or mm.PME. Default: mm.NoCutoff

  .. attribute::        nonbonded_cutoff

    Nonbonded cutoff set when simulation is set up. Default: 10.0

  .. attribute::        remove_cmm_motion

    Whether a center of mass remover object should be attached to the
    :class:`Simulation`. Default: True

  .. attribute::        cmm_frequency

    Frequency regarding simulation steps when the cmm remover should be applied.
    Default: 1

  .. attribute::        periodic_box_extent

    :class:`ost.geom.Vec3` describing the X,Y and Z extents of the rectangular
    simulation_box defined when setting up the :class:`Simulation`. 
    Default: (0,0,0)

  .. attribute::        init_temperature

    Initial simulation velocities get set according to a Boltzman distribution
    controlled by init_temperature(Kelvin). Default: 0.0

  .. attribute::        forcefield

    :class:`Forcefield` that is used in :class:`TopologyCreator`. Default: None

  .. attribute::        termini_exceptions

    Use termini other than the defaults set in the :class:`Forcefield` when
    using :class:`TopologyCreator`. Has to be a :class:`TerminiExceptions`
    object. Default: None

  .. attribute::        platform

    Platform used by OpenMM to do the calculations. Must be one of mm.Reference,
    mm.CPU, mm.CUDA or mm.OpenCL. If anything else than the reference platform 
    is used, the attribute openmm_plugin_directory has to be set accordingly.
    Default: mm.Reference

  .. attribute::        reference_properties

    :class:`dict` of OpenMM specific properties that can be set for the 
    reference platform.
       
  .. attribute::        cpu_properties

    :class:`dict` of OpenMM specific properties that can be set for the 
    cpu platform.

  .. attribute::        opencl_properties

    :class:`dict` of OpenMM specific properties that can be set for the 
    opencl platform.

  .. attribute::        cuda_properties

    :class:`dict` of OpenMM specific properties that can be set for the 
    cuda platform.

  .. attribute::        add_thermostat

    Flag, whether an Andersen thermostat should be attached when settings up
    the :class:`Simulation`. Default: False

  .. attribute::        thermostat_temperature

    Temperature for the Andersen thermostat in K. Default: NaN

  .. attribute::        thermostat_collision_frequency

    Collision frequency of the Andersen thermostat in 1/ps. Default: NaN

  .. attribute::        add_barostat

    Flag, whether an  MonteCarlo barostat should be attached when setting up
    the :class:`Simulation`. Default: False

  .. attribute::        barostat_temperature

    Temperature for the MonteCarlo Barostat in K. Default: NaN

  .. attribute::        barostat_pressure

    Pressure for the MonteCarlo Barostat in bar. Default: NaN

  .. attribute::        barostat_frequency

    Frequency of the MonteCarlo Barostat. Default: 25

  .. attribute::        integrator

    Integrator to move the simulation forward in time. OpenMM offers following
    Integrators: :class:`VerletIntegrator`, :class:`BrownianIntegrator`,
    :class:`LangevinIntegrator`, :class:`VariableVerletIntegrator` and
    :class:`VariableLangevinIntegrator`. Default: None

  .. attribute::        solvent_dielectric

    Solvent dielectric constant. This is used for the GBSA force when the 
    :class:`Simulation` gets set up. Default: 78.3

  ..  attribute::       solute_dielectric

    Solute dielectric constant. This is used for the GBSA force when the
    :class:`Simulation` gets set up. Default: 1.0

  .. attribute::        reaction_field_dielecric

    Sets the reaction field dielectric for the Nonbonded Force when setting
    up the :class:`Simulation`. Default: 78.3

  .. attribute::        use_dispersion_correction

    Flag, whether the dispersion correction should be used when setting
    up the Nonbonded Force in the :class:`Simulation`. Default: True

  .. attribute::        keep_ff_specific_naming

    When running through the :class:`TopologyCreator`, the given entity gets 
    renamed to the forcefield specific naming. If set to true, the naming
    stays like that, if not the entity gets renamed to PDB standard.
    Default: True

  .. attribute::        openmm_plugin_directory

    Path where OpenMM specific plugins are searched. If you want
    to use other platforms than Reference, this has to be set.
    Defaults to *OPEN_MM_PLUGIN_DIR* which was set when
    :ref:`configuring the compilation <cmake-flags>`.

  .. attribute::        custom_plugin_directory

    Path where custom plugins are searched for. Defaults to
    "share/openstructure/openmm_plugins" within the OST installation or to
    :attr:`openmm_plugin_directory` if the OST path could not be determined.




.. class:: TerminiExceptions
  
  Can be used to define exceptions from the standard assignments
  assigned in the :class:`Forcefield`. 

  .. method::           SetException(residue, exception_name)

    :param residue:     Residue for which the exception
                        should be set.
    :param exception_name: Name of the
                        Blockmodifier in the forcefield that should be
                        applied in case of a termini

    :type residue:      :class:`ResidueHandle`
    :type exception_name: :class:`str`

  .. method::           HasException(residue)

    :param residue:     Residue that should be checked for 
                        exceptions
    :type residue:      :class:`ResidueHandle`

    :returns:           :class:`bool` True if **residue** has an exception

  .. method::           GetException(residue)

    :param residue:     Residue containing the exception
    :type residue:      :class:`ResidueHandle`

    :raises:            :class:`RuntimeError` if **residue** has no associated exception

    :returns:           :class:`str` describing the name of the Blockmodifier in
                        the forcefield that should be applied in case of a 
                        terminal residue
