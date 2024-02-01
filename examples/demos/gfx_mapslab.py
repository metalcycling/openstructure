import math,random
from ost import img

# remove all objects from scene, just in case
scene.RemoveAll()

vmax=-10000.0
vmin=+10000.0
mh=img.CreateMap(img.Size(32,32,32))
for p in img.ExtentIterator(mh.GetExtent()):
  val=5*math.sin(0.4*math.sqrt(p[0]*p[0]+p[1]*p[1]))+7*math.cos(0.6*math.sqrt(p[2]*p[2]+p[1]*p[1]))
  mh.SetReal(p,val)
  vmin=min(vmin,val)
  vmax=max(vmax,val)
print(vmin, vmax)
for p in img.ExtentIterator(mh.GetExtent()):
  mh.SetReal(p,(mh.GetReal(p)-vmin)/(vmax-vmin))

pl = gfx.PrimList("box")
edge_list=[[[0,0,0],[0,1,0]],
           [[0,1,0],[1,1,0]],
           [[1,1,0],[1,0,0]],
           [[1,0,0],[0,0,0]],
           [[0,0,0],[0,0,1]],
           [[0,1,0],[0,1,1]],
           [[1,1,0],[1,1,1]],
           [[1,0,0],[1,0,1]],
           [[0,0,1],[0,1,1]],
           [[0,1,1],[1,1,1]],
           [[1,1,1],[1,0,1]],
           [[1,0,1],[0,0,1]]]

corners = [mh.IndexToCoord(mh.GetExtent().GetStart()),mh.IndexToCoord(mh.GetExtent().GetEnd())]

for e in edge_list:
    pl.AddLine(geom.Vec3(corners[e[0][0]][0],
                    corners[e[0][1]][1],
                    corners[e[0][2]][2]),
               geom.Vec3(corners[e[1][0]][0],
                    corners[e[1][1]][1],
                    corners[e[1][2]][2]),
               gfx.Color(1,1,0))

scene.Add(pl)
scene.SetCenter(pl.GetCenter())

go1=gfx.MapIso("iso", mh,0.5)
go1.SetLineWidth(1.5)
scene.Add(go1)
scene.SetCenter(go1.GetCenter())


go2 = gfx.MapSlab("slab",mh,geom.Plane(go1.GetCenter(),geom.Vec3(0.0,0.0,1.0)))
scene.Add(go2)

go2.ColorBy('HEAT_MAP',0.2,0.8)

print('Demo 4: Projecting the density of a map onto a plane...') 
