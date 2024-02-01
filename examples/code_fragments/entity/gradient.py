lines=gfx.PrimList("lines")

g=gfx.Gradient()
g.SetColorAt(0, gfx.Color(1.0, 0.0, 0.0, 1.0))
g.SetColorAt(0.5, gfx.Color(0.0, 1.0, 0.0, 1.0))
g.SetColorAt(1.0, gfx.Color(0.0, 0.0, 1.0, 1.0))

for i in range(1000):
  c=g.GetColorAt(0.001*i)
  lines.AddLine(geom.Vec3(i*0.05, 0.0, -1), geom.Vec3(i*0.05, 0.0, 1), c)

scene.Add(lines)
scene.center=geom.Vec3(25, 0, 0)
