from ost.mol import mm
from PyQt5 import QtCore

"""
MM - Demo

!!!Requires OpenStructure to be compiled with MM-Support!!!

This demo should show the setup of a simple simulation using implicit solvent.
The script can directly be started from within DNG.
"""

class Anim(QtCore.QTimer):
    def __init__(self,sim,ent,go):
        QtCore.QTimer.__init__(self)
        self.ent=ent
        self.sim=sim
        self.go = go
        self.ed = ent.EditXCS(ost.mol.BUFFERED_EDIT)
        self.timeout.connect(self.OnTimer)

        
    def OnTimer(self):
        self.sim.Steps(5)
        self.sim.UpdatePositions()
        self.go.UpdatePositions()


prot=ost.mol.CreateEntityFromView(io.LoadPDB('1CRN.pdb').Select('peptide=true'),True)

settings = mm.Settings()

settings.integrator = mm.LangevinIntegrator(310,1,0.002)
settings.add_gbsa = True
settings.forcefield = mm.LoadCHARMMForcefield()
settings.platform = mm.Platform.CPU
settings.nonbonded_method = mm.NonbondedMethod.CutoffNonPeriodic

sim = mm.Simulation(prot,settings)
sim.ApplyLBFGS(tolerance = 1.0, max_iterations = 200)
sim.UpdatePositions()

ent = sim.GetEntity()
go = gfx.Entity("crambin",ent)
scene.Add(go)

anim = Anim(sim,ent,go)
anim.start(1)

