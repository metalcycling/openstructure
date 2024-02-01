The mm Module
================================================================================

.. module:: ost.mol.mm
   :synopsis: Wrapper for OpenMM molecular mechanics simulations

.. currentmodule:: ost.mol.mm

Introduction
--------------------------------------------------------------------------------

The mol.mm module provides a wrapper around the 
`OpenMM <http://openmm.org>`_ molecular mechanics library to provide
basic molecular dynamics (MD) capabilities fully embedded into the OpenStructure universe.
The heart of every simulation is the :class:`Topology` describing how the 
particles of an :class:`EntityHandle` interact. The :class:`Simulation` class connects the :class:`EntityHandle` with a :class:`Topology` and allows you to perform energy minimizations 
or move the simulation through time using an :class:`Integrator`. 
A :class:`Topology` can either be built from scratch by adding
one interaction after the other or automatically using the
:class:`TopologyCreator`. The process of :class:`Topology` building
and setting up a :class:`Simulation` gets controlled with the :class:`Settings`.

Latest Publication of `OpenMM <http://openmm.org>`_:
P. Eastman, J. Swails, J. D. Chodera, R. T. McGibbon, Y. Zhao, K. A. Beauchamp, L.-P. Wang, A. C. Simmonett, M. P. Harrigan, C. D. Stern, R. P. Wiewiora, B. R. Brooks, and V. S. Pande. "OpenMM 7: Rapid development of high performance algorithms for molecular dynamics." PLOS Comp. Biol. 13(7): e1005659. (2017) 


Installation
--------------------------------------------------------------------------------

If you are compiling OpenStructure from source, you need to specifically enable
support for `OpenMM <http://openmm.org>`_. You do this by enabling the
``ENABLE_MM`` flag and setting the ``OPEN_MM_INCLUDE_DIR``, ``OPEN_MM_LIBRARY``
and ``OPEN_MM_PLUGIN_DIR`` flags when calling ``cmake``. See
:ref:`here <cmake-flags>` for details and examples.


Setting up a simple simulation
--------------------------------------------------------------------------------


  .. code-block:: python
    
    from ost.mol import mm

    prot=io.LoadPDB('1crn',remote=True)

    #set up the simulation
    settings = mm.Settings()
    settings.integrator = mm.LangevinIntegrator(310,1,0.002)
    settings.forcefield = mm.LoadCHARMMForcefield()
    sim = mm.Simulation(prot,settings)

    #minimize it
    sim.ApplySD(tolerance = 1.0, max_iterations = 200)

    #create a trajectory observer and register it to the simulation
    #every 10 steps, the actual positions will be written down to disk
    observer = mm.TrajWriter(10,"example_traj.pdb","example_traj.dcd")
    sim.Register(observer)

    #run the simulation
    sim.Steps(10000)

    #Trajectory Observer needs to finalize, otherwise you might get a corrupt dcd file
    observer.Finalize()


Doing more sophisticated stuff
--------------------------------------------------------------------------------

You want to create your own :class:`BuildingBlock` to parameterize custom 
residues? Or even generate your own custom :class:`Forcefield`? 
Check out the mm dir in the examples/code_fragments directory. 

Documentation
--------------------------------------------------------------------------------

.. toctree::
  interaction
  buildingblock
  forcefield
  settings
  topology
  observers
  integrators
  simulation

