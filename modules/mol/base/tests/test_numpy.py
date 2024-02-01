import unittest
if __name__== '__main__':
  import sys
  sys.path.insert(0,"../../../../stage/lib64/openstructure/")
  sys.path.insert(0,"../../../../stage/lib/openstructure/")
  
import ost
from ost import geom, mol

if ost.WITH_NUMPY:
  has_numpy=True
  try:
    import numpy
  except ImportError as e:
    has_numpy=False
else:
  has_numpy=False

def v2v(v):
  return geom.Vec3(float(v[0]),float(v[1]),float(v[2]))

def dd(v1,v2):
  return geom.Distance(v1,v2)<1e-8

class TestNumpy(unittest.TestCase):
  def setUp(self):
    pass

  def test_(self):
    if not has_numpy:
      return
    entity=mol.CreateEntity()
    ed=entity.EditXCS()
    ch=ed.InsertChain("X")
    re=ed.AppendResidue(ch,"ALA")
    a0=ed.InsertAtom(re,"A",geom.Vec3(0,0,0))
    self.assertEqual(a0.GetIndex(),0)
    a1=ed.InsertAtom(re,"B",geom.Vec3(1,0,0))
    self.assertEqual(a1.GetIndex(),1)
    a2=ed.InsertAtom(re,"C",geom.Vec3(2,0,0))
    self.assertEqual(a2.GetIndex(),2)
    a3=ed.InsertAtom(re,"D",geom.Vec3(3,0,0))
    self.assertEqual(a3.GetIndex(),3)

    self.assertTrue(dd(a0.pos,geom.Vec3(0,0,0)))
    self.assertTrue(dd(a1.pos,geom.Vec3(1,0,0)))
    self.assertTrue(dd(a2.pos,geom.Vec3(2,0,0)))
    self.assertTrue(dd(a3.pos,geom.Vec3(3,0,0)))

    ed.SetAtomTransformedPos(entity.GetAtomList(),
                             numpy.array([[0,1,0],[0,2,0],[0,3,0],[0,4,0]], dtype=numpy.float32))

    self.assertTrue(dd(a0.pos,geom.Vec3(0,1,0)))
    self.assertTrue(dd(a1.pos,geom.Vec3(0,2,0)))
    self.assertTrue(dd(a2.pos,geom.Vec3(0,3,0)))
    self.assertTrue(dd(a3.pos,geom.Vec3(0,4,0)))

    na=entity.positions

    self.assertTrue(dd(v2v(na[0]),geom.Vec3(0,1,0)))
    self.assertTrue(dd(v2v(na[1]),geom.Vec3(0,2,0)))
    self.assertTrue(dd(v2v(na[2]),geom.Vec3(0,3,0)))
    self.assertTrue(dd(v2v(na[3]),geom.Vec3(0,4,0)))

    ed.SetAtomTransformedPos([3,99,2],
                             numpy.array([[0,0,-3],[-1,-1,-1],[0,0,-2]], dtype=numpy.float32))

    self.assertTrue(dd(a0.pos,geom.Vec3(0,1,0)))
    self.assertTrue(dd(a1.pos,geom.Vec3(0,2,0)))
    self.assertTrue(dd(a2.pos,geom.Vec3(0,0,-2)))
    self.assertTrue(dd(a3.pos,geom.Vec3(0,0,-3)))

if __name__== '__main__':
    unittest.main()

