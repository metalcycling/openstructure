import math

e=io.LoadEntity('../entity/sdh.pdb')
c=e.FindChain('A')

edi=e.EditICS(mol.EditMode.BUFFERED_EDIT)
for r in c.GetResidueList():
  phi=r.GetPhiTorsion()
  psi=r.GetPsiTorsion()
  if phi.IsValid(): 
    edi.SetTorsionAngle(phi, math.radians(-60))    
  if psi.IsValid(): 
    edi.SetTorsionAngle(psi, math.radians(-45))    
  
edi.UpdateXCS()

g=gfx.Entity('helix', e.Select('aname=CA,N,O,C,CB'))
scene.Add(g)
scene.SetCenter(g.GetCenter())



