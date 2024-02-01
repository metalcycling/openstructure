from ost.mol.mm import  *
from PyQt5 import QtCore
from math import sqrt
from math import sin

"""
MM - Demo

!!!Requires OpenStructure to be compiled with MM-Support!!!

This demo should show the procedure building a topology from scratch and
setting up a simulation. As an additional piece of sugar, it demonstrates
the capabilities of altering force parameters as the simulation runs.
The script can directly be started from within DNG.
"""

class Anim(QtCore.QTimer):
    def __init__(self,sim,go,cc_bond_index):
        QtCore.QTimer.__init__(self)
        self.ent=sim.GetEntity()
        self.sim=sim
        self.go = go
        self.ed = ent.EditXCS()
        self.bond_index = cc_bond_index
        top = sim.GetTopology()
        param = top.GetHarmonicBondParameters(cc_bond_index)
        self.bond_length = param[2]
        self.force_constant = param[3]
        self.steps = 0
        self.counter = 0
        self.timeout.connect(self.OnTimer)

    def OnTimer(self):
        if self.steps % 20 == 0:
          sim.ResetHarmonicBond(self.bond_index,sin(self.counter*3.1416/100)*self.bond_length*0.75 + self.bond_length,self.force_constant)
          self.counter += 1

        self.sim.Steps(1)
        self.sim.UpdatePositions()
        self.go.UpdatePositions()
        self.steps+=1


#create topology by only defining masses

prof = io.IOProfile(dialect='PDB', fault_tolerant=False,
                 quack_mode=False, processor=conop.HeuristicProcessor())

ent = io.LoadPDB('ethanol.pdb', profile=prof)
masses = [12.011,12.011,15.999,1.008,1.008,1.008,1.008,1.008,1.008]
top = Topology(masses)

#all per atom parameters have to be set at once...
charges = [-0.18,0.145,-0.683,0.06,0.06,0.06,0.06,0.06,0.418]
sigmas = [0.37,0.37,0.32,0.24,0.24,0.24,0.24,0.24,0.04]
epsilons = [0.28,0.28,0.64,0.09,0.09,0.09,0.09,0.09,0.2]

top.SetCharges(charges)
top.SetSigmas(sigmas)
top.SetEpsilons(epsilons)

#we can set the fudge parameters for 1,4 interactions also in the topology
top.SetFudgeQQ(0.5)
top.SetFudgeLJ(0.5)

#add the bonds
top.AddHarmonicBond(0,1,0.153,186188.0) #C-C
top.AddHarmonicBond(0,3,0.1111,269449.6) #C-H
top.AddHarmonicBond(0,4,0.1111,269449.6) #C-H
top.AddHarmonicBond(0,5,0.1111,269449.6) #C-H
top.AddHarmonicBond(1,6,0.1111,269449.6) #C-H
top.AddHarmonicBond(1,7,0.1111,269449.6) #C-H
top.AddHarmonicBond(1,2,0.142,358150.4) #C-O
top.AddHarmonicBond(2,8,0.096,456056.0) #O-H

#add the angles
top.AddUreyBradleyAngle(3,0,4,1.892,297.064,0.1802,4518.72) #H-C-H
top.AddUreyBradleyAngle(3,0,5,1.892,297.064,0.1802,4518.72) #H-C-H
top.AddUreyBradleyAngle(4,0,5,1.892,297.064,0.1802,4518.72) #H-C-H
top.AddUreyBradleyAngle(3,0,1,1.922,289.5328,0.2179,18853.104) #H-C-C
top.AddUreyBradleyAngle(4,0,1,1.922,289.5328,0.2179,18853.104) #H-C-C
top.AddUreyBradleyAngle(5,0,1,1.922,289.5328,0.2179,18853.104) #H-C-C
top.AddUreyBradleyAngle(0,1,6,1.922,289.5328,0.2179,18853.104) #H-C-C
top.AddUreyBradleyAngle(0,1,7,1.922,289.5328,0.2179,18853.104) #H-C-C
top.AddUreyBradleyAngle(0,1,2,1.920,633.4576,0.0,0.0) #C-C-O
top.AddUreyBradleyAngle(6,1,2,1.900,384.0912,0.0,0.0) #H-C-O
top.AddUreyBradleyAngle(7,1,2,1.900,384.0912,0.0,0.0) #H-C-O
top.AddUreyBradleyAngle(6,1,7,1.892,297.064,0.1802,4518.72) #H-C-H
top.AddUreyBradleyAngle(1,2,8,1.850,481.16,0.0,0.0) #C-O-H


#add the dihedrals
top.AddPeriodicDihedral(3,0,1,6,3,0.0,0.66944) #H-C-C-H
top.AddPeriodicDihedral(3,0,1,7,3,0.0,0.66944) #H-C-C-H
top.AddPeriodicDihedral(3,0,1,2,3,0.0,0.66944) #H-C-C-O
top.AddPeriodicDihedral(4,0,1,6,3,0.0,0.66944) #H-C-C-H
top.AddPeriodicDihedral(4,0,1,7,3,0.0,0.66944) #H-C-C-H
top.AddPeriodicDihedral(4,0,1,2,3,0.0,0.66944) #H-C-C-O
top.AddPeriodicDihedral(5,0,1,6,3,0.0,0.66944) #H-C-C-H
top.AddPeriodicDihedral(5,0,1,7,3,0.0,0.66944) #H-C-C-H
top.AddPeriodicDihedral(5,0,1,2,3,0.0,0.66944) #H-C-C-O
top.AddPeriodicDihedral(6,1,2,8,3,0.0,0.58576) #H-C-O-H
top.AddPeriodicDihedral(7,1,2,8,3,0.0,0.58576) #H-C-O-H
top.AddPeriodicDihedral(0,1,2,8,1,0.0,5.4392) #C-C-O-H
top.AddPeriodicDihedral(0,1,2,8,2,0.0,1.2552) #C-C-O-H
top.AddPeriodicDihedral(0,1,2,8,3,0.0,1.75728) #C-C-O-H

#we have to manually set the exclusions for the nonbonded interaction...

#1,2 interactions
top.AddExclusion(3,0)
top.AddExclusion(4,0)
top.AddExclusion(5,0)
top.AddExclusion(0,1)
top.AddExclusion(1,6)
top.AddExclusion(1,7)
top.AddExclusion(1,2)
top.AddExclusion(2,8)

#1,3 interactions
top.AddExclusion(3,4)
top.AddExclusion(3,5)
top.AddExclusion(4,5)
top.AddExclusion(3,1)
top.AddExclusion(4,1)
top.AddExclusion(5,1)
top.AddExclusion(0,6)
top.AddExclusion(0,7)
top.AddExclusion(0,2)
top.AddExclusion(6,2)
top.AddExclusion(7,2)
top.AddExclusion(1,8)

#add lj-pairs, note, that the fudge parameters still get applied during system setup
top.AddLJPair(3,6,0.5*(sigmas[3]+sigmas[6]),sqrt(epsilons[3]*epsilons[6]))
top.AddLJPair(3,7,0.5*(sigmas[3]+sigmas[7]),sqrt(epsilons[3]*epsilons[7]))
top.AddLJPair(3,2,0.5*(sigmas[3]+sigmas[2]),sqrt(epsilons[3]*epsilons[2]))
top.AddLJPair(4,6,0.5*(sigmas[4]+sigmas[6]),sqrt(epsilons[4]*epsilons[6]))
top.AddLJPair(4,7,0.5*(sigmas[4]+sigmas[7]),sqrt(epsilons[4]*epsilons[7]))
top.AddLJPair(4,2,0.5*(sigmas[4]+sigmas[2]),sqrt(epsilons[4]*epsilons[2]))
top.AddLJPair(5,6,0.5*(sigmas[5]+sigmas[6]),sqrt(epsilons[5]*epsilons[6]))
top.AddLJPair(5,7,0.5*(sigmas[5]+sigmas[7]),sqrt(epsilons[5]*epsilons[7]))
top.AddLJPair(5,2,0.5*(sigmas[5]+sigmas[2]),sqrt(epsilons[5]*epsilons[2]))
top.AddLJPair(6,8,0.5*(sigmas[6]+sigmas[8]),sqrt(epsilons[6]*epsilons[8]))
top.AddLJPair(7,8,0.5*(sigmas[7]+sigmas[8]),sqrt(epsilons[7]*epsilons[8]))
top.AddLJPair(0,8,0.5*(sigmas[0]+sigmas[8]),sqrt(epsilons[0]*epsilons[8])) 

settings = Settings()
settings.init_temperature =  310
settings.integrator = VerletIntegrator(0.0001)

#let's extract the cc_bond_index to have some fun with it
cc_bond_index = top.GetHarmonicBondIndices(0,1)[0]

#create the simulation
sim = Simulation(top,ent,settings)
go = gfx.Entity("ethanol",sim.GetEntity())
go.SetRenderMode(gfx.CUSTOM)
scene.Add(go)

anim = Anim(sim,go,cc_bond_index)
anim.start(2)


