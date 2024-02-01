sdh=io.LoadPDB('sdh.pdb')
helix=sdh.Select('rnum=75:78 and cname=A and aname=CA,C,N,O')
res_a=helix.residues[0]
res_b=helix.residues[-1]

def MakePlane(name, pt1, pt2, pt3):
  p1=[pt1, pt2, pt3, pt2+(pt1-pt2)+(pt3-pt2)]  
  c1=sum(p1)/4.0
  p1=[c1+(p-c1)*1.5 for p in p1]
  plane=gfx.Quad(name, *p1)
  plane.SetLineWidth(1.0)
  plane.SetFillColor(gfx.Color(0.3, 0.7, 0.3, 0.5))  
  scene.Add(plane)

r1_n=res_a.FindAtom("N").pos
r1_ca=res_a.FindAtom("CA").pos
r1_c=res_a.FindAtom("C").pos
MakePlane("Plane1", r1_n, r1_ca, r1_c)
r2_ca=res_b.FindAtom("CA").pos
r2_c=res_b.FindAtom("C").pos
r2_o=res_b.FindAtom("O").pos
MakePlane("Plane2", r2_ca, r2_c, r2_o)

ca_dist=gfx.PrimList("Calpha dist")
ca_dist.AddLine(r1_ca, r2_ca, gfx.GREEN)
scene.Add(ca_dist)
go=gfx.Entity('helix2', gfx.SIMPLE, helix)
scene.Add(go)
scene.SetBackground(gfx.WHITE)
go.SetRenderMode(gfx.CUSTOM)
scene.center=go.center