from ost.mol.mm import *
from PyQt5 import QtCore

"""
MM - Demo

!!!Requires OpenStructure to be compiled with MM-Support!!!

This demo should show the procedure of building a custom forcefield.
The script can directly be started from within DNG.
"""

class Anim(QtCore.QTimer):
    def __init__(self,sim,go):
        QtCore.QTimer.__init__(self)
        self.ent=sim.GetEntity()
        self.sim=sim
        self.go = go
        self.ed = ent.EditXCS()
        self.timeout.connect(self.OnTimer)
        
    def OnTimer(self):
        self.sim.Steps(1)
        self.sim.UpdatePositions()
        self.go.UpdatePositions()


#Let's create an ethanol without its hydrogens
ent = ost.mol.CreateEntity()
ed = ent.EditXCS()
chain = ed.InsertChain('A')
res = ed.AppendResidue(chain,"ETH")
c1 = ed.InsertAtom(res,"C1",geom.Vec3(-0.246,0.0,1.285),"C")
c2 = ed.InsertAtom(res,"C2",geom.Vec3(0.566,0.0,-0.011),"C")
o = ed.InsertAtom(res,"O",geom.Vec3(-0.322,0.0,-1.130),"O")
ed.Connect(c1,c2)
ed.Connect(c2,o)


ff = Forcefield()

#generate hydrogen constructor and add it to the forcefield
h_constructor = GromacsHydrogenConstructor()
h_constructor.AddHydrogenRule(3,4,["H11","H12","H13"],["C1","C2","O"])
h_constructor.AddHydrogenRule(2,6,["H21","H22"],["C2","O","C1"])
h_constructor.AddHydrogenRule(1,2,["OH"],["O","C2","C1"])
ff.AddHydrogenConstructor("ETH",h_constructor)

#generate buildingblock and add it to the forcefield
building_block = BuildingBlock()
building_block.AddAtom("C1","C",-0.18)
building_block.AddAtom("C2","C",0.145)
building_block.AddAtom("O","O",-0.683)
building_block.AddAtom("H11","H",0.06)
building_block.AddAtom("H12","H",0.06)
building_block.AddAtom("H13","H",0.06)
building_block.AddAtom("H21","H",0.06)
building_block.AddAtom("H22","H",0.06)
building_block.AddAtom("OH","OH",0.418)

bonds = list()
for i in range(8):
  bonds.append(Interaction(FuncType.HarmonicBond))

bonds[0].SetNames(["H11","C1"])
bonds[1].SetNames(["H12","C1"])
bonds[2].SetNames(["H13","C1"])
bonds[3].SetNames(["C1","C2"])
bonds[4].SetNames(["C2","H21"])
bonds[5].SetNames(["C2","H22"])
bonds[6].SetNames(["O","C2"])
bonds[7].SetNames(["O","OH"])

for b in bonds:
  building_block.AddBond(b)

ff.AddBuildingBlock("ETH",building_block)

#add masses to the forcefield
ff.AddMass("C",12.011)
ff.AddMass("O",15.999)
ff.AddMass("H",1.008)
ff.AddMass("OH",1.008)

#generate bonds and add them to the forcefield
h_c_bond = Interaction(FuncType.HarmonicBond)
c_c_bond = Interaction(FuncType.HarmonicBond)
c_o_bond = Interaction(FuncType.HarmonicBond)
o_h_bond = Interaction(FuncType.HarmonicBond)

h_c_bond.SetTypes(["H","C"])
c_c_bond.SetTypes(["C","C"])
c_o_bond.SetTypes(["C","O"])
o_h_bond.SetTypes(["O","OH"])

h_c_bond.SetParam([0.1111,269449.6])
c_c_bond.SetParam([0.153,186188.0])
c_o_bond.SetParam([0.142,358150.4])
o_h_bond.SetParam([0.096,456056.0])

ff.AddBond(h_c_bond)
ff.AddBond(c_c_bond)
ff.AddBond(c_o_bond)
ff.AddBond(o_h_bond)

#generate angles and add them to the forcefield
h_c_h_angle = Interaction(FuncType.UreyBradleyAngle)
h_c_c_angle = Interaction(FuncType.UreyBradleyAngle)
h_c_o_angle = Interaction(FuncType.UreyBradleyAngle)
c_c_o_angle = Interaction(FuncType.UreyBradleyAngle)
c_o_h_angle = Interaction(FuncType.UreyBradleyAngle)

h_c_h_angle.SetTypes(["H","C","H"])
h_c_c_angle.SetTypes(["H","C","C"])
h_c_o_angle.SetTypes(["H","C","O"])
c_c_o_angle.SetTypes(["C","C","O"])
c_o_h_angle.SetTypes(["C","O","OH"])

h_c_h_angle.SetParam([1.892,297.064,0.1802,4518.72])
h_c_c_angle.SetParam([1.922,289.5328,0.2179,18853.104])
h_c_o_angle.SetParam([1.900,384.0912,0.0,0.0])
c_c_o_angle.SetParam([1.920,633.4576,0.0,0.0])
c_o_h_angle.SetParam([1.850,481.16,0.0,0.0])

ff.AddAngle(h_c_h_angle)
ff.AddAngle(h_c_c_angle)
ff.AddAngle(h_c_o_angle)
ff.AddAngle(c_c_o_angle)
ff.AddAngle(c_o_h_angle)

#generate dihedrals and add them to the forcefield
x_c_c_x_dihedral = Interaction(FuncType.PeriodicDihedral)
c_c_o_h_dihedral_one = Interaction(FuncType.PeriodicDihedral)
c_c_o_h_dihedral_two = Interaction(FuncType.PeriodicDihedral)
c_c_o_h_dihedral_three = Interaction(FuncType.PeriodicDihedral)
x_c_o_x_dihedral = Interaction(FuncType.PeriodicDihedral)

x_c_c_x_dihedral.SetTypes(["X","C","C","X"])
c_c_o_h_dihedral_one.SetTypes(["C","C","O","OH"])
c_c_o_h_dihedral_two.SetTypes(["C","C","O","OH"])
c_c_o_h_dihedral_three.SetTypes(["C","C","O","OH"])
x_c_o_x_dihedral.SetTypes(["X","C","O","X"])

x_c_c_x_dihedral.SetParam([3,0.0,0.66944])
c_c_o_h_dihedral_one.SetParam([1,0.0,5.4392])
c_c_o_h_dihedral_two.SetParam([2,0.0,1.2552])
c_c_o_h_dihedral_three.SetParam([3,0.0,1.75728])
x_c_o_x_dihedral.SetParam([3,0.0,0.58576])

ff.AddDihedral(x_c_c_x_dihedral)
ff.AddDihedral(c_c_o_h_dihedral_one)
ff.AddDihedral(c_c_o_h_dihedral_two)
ff.AddDihedral(c_c_o_h_dihedral_three)
ff.AddDihedral(x_c_o_x_dihedral)

#add lj parameters
c_lj = Interaction(FuncType.LJ)
o_lj = Interaction(FuncType.LJ)
h_lj = Interaction(FuncType.LJ)
oh_lj = Interaction(FuncType.LJ)

c_lj.SetTypes(["C"])
o_lj.SetTypes(["O"])
h_lj.SetTypes(["H"])
oh_lj.SetTypes(["OH"])

c_lj.SetParam([0.37,0.28])
o_lj.SetParam([0.32,0.64])
h_lj.SetParam([0.24,0.09])
oh_lj.SetParam([0.04,0.2])

ff.AddLJ(c_lj)
ff.AddLJ(o_lj)
ff.AddLJ(h_lj)
ff.AddLJ(oh_lj)

#we do not explicitely set the 1,4 interactions,
#we rather decrease it by a certain factor
ff.SetFudgeLJ(0.5)
ff.SetFudgeQQ(0.5)

#construct settings
settings = Settings()
settings.init_temperature =  310
settings.forcefield = ff

#construct integrator with tiny timesteps
settings.integrator = VerletIntegrator(0.0001)


#create the simulation
sim = Simulation(ent,settings)

io.SavePDB(sim.GetEntity(),'ethanol.pdb')


go = gfx.Entity("ethanol",sim.GetEntity())
go.SetRenderMode(gfx.CUSTOM)
scene.Add(go)

anim = Anim(sim,go)
anim.start(2)


