from ost.mol.mm import *

"""
MM - Demo

!!!Requires OpenStructure to be compiled with MM-Support!!!

This demo should show the procedure of setting up a simulation and write
the results directly to disk.
"""


prot=io.LoadPDB('1CRN.pdb')

#set up the simulation
settings = Settings()
settings.integrator = LangevinIntegrator(310,1,0.002)
settings.add_gbsa = True
settings.forcefield = LoadCHARMMForcefield()
settings.nonbonded_cutoff = 8.0
settings.nonbonded_method = NonbondedMethod.CutoffNonPeriodic
settings.platform = Platform.CPU
settings.cpu_properties["CpuThreads"] = "2"

sim = Simulation(prot,settings)

#minimize it
sim.ApplySD(tolerance = 1.0, max_iterations = 200)

#create a trajectory observer and register it to the simulation
observer = TrajWriter(10,"gb_example_traj.pdb","gb_example_traj.dcd")
sim.Register(observer)

#run the simulation
sim.Steps(10000)

#Trajectory Observer needs to finalize, otherwise you might get a corrupt dcd file
observer.Finalize()




