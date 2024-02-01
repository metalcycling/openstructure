from ost.gui import trajectory_viewer

"""
MM - Demo

!!!Requires OpenStructure to be compiled with MM-Support!!!

This demo demonstrates how to actually look at the trajectory created
in the gb_example_writing_trajectory.py script.
"""

t = io.LoadCHARMMTraj("gb_example_traj.pdb","gb_example_traj.dcd")
go = gfx.Entity("gb",t.GetEntity())
scene.Add(go)
w = trajectory_viewer.TrajWidget([t],[go])
w.show()
