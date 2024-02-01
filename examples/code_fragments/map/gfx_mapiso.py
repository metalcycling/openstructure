import math,random
from ost import img
mh=img.CreateMap(img.Size(100, 25, 25))
for p in img.ExtentIterator(mh.GetExtent()):
  val=5*math.sin(0.4*math.sqrt(p[0]*p[0]+p[1]*p[1]))+7*math.cos(0.6*math.sqrt(p[2]*p[2]+p[1]*p[1]))
  mh.SetReal(p,val*(1.0+0.0*random.random()))

go=gfx.MapIso("test",mh,4.0)
scene.Add(go)
scene.SetCenter(go.GetCenter())
go.SetLineWidth(1.5)
