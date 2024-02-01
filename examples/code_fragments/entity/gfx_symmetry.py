import math
from ost import geom

filename='fragment.pdb'
if len(sys.argv)>1:
  filename=sys.argv[1]
ent=io.LoadEntity(filename)

edi=ent.EditXCS()
m=geom.Mat4()
m.PasteTranslation(-ent.bounds.min)
edi.ApplyTransform(m)
edi.UpdateICS()
frag=gfx.Entity('frag', ent)
sym_ops=gfx.SymmetryOpList()
for i in range(12):
  m=geom.EulerTransformation(0, i*2.0*math.pi/12.0, 0)
  sym_ops.append(gfx.SymmetryOp(m, geom.Vec3(0, 0, 10)))
  
sym=gfx.SymmetryNode('sym', sym_ops)
scene.Add(sym)
sym.Add(frag)

scene.SetCenter(geom.Vec3())
