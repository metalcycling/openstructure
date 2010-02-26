sdh=io.LoadPDB('sdh.pdb')
helix=sdh.Select('rnum=99:128 and cname=A and aname=CA,C,N,O')
go=gfx.Entity('helix', gfx.SIMPLE, helix)
scene.Add(go)

bbox=mol.BoundingBoxFromEntity(helix)

def RenderBBox(bbox):
  bb=gfx.Cuboid('xxx', bbox)
  bb.SetFillColor(gfx.Color(0.5, 0.8, 0.5, 0.2))
  scene.Add(bb)

print 'C1:',(helix.GetGeometricStart()+helix.GetGeometricEnd())*.5
  
RenderBBox(bbox)
scene.center=go.center