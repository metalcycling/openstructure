#------------------------------------------------------------------------------
# This file is part of the OpenStructure project <www.openstructure.org>
#
# Copyright (C) 2008-2020 by the OpenStructure authors
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#------------------------------------------------------------------------------

import unittest
if __name__== '__main__':
  import sys
  sys.path.insert(0,"../../../stage/lib64/openstructure/")
  sys.path.insert(0,"../../../stage/lib/openstructure/")

import ost
import ost.mol as mol
import ost.gfx as gfx
import ost.geom as geom

if ost.WITH_NUMPY:
  has_numpy=True
  try:
    import numpy
  except ImportError as e:
    has_numpy=False
else:
  has_numpy=False

def col_delta(c1,c2):
  return geom.Distance(geom.Vec3(c1[0],c1[1],c1[2]),geom.Vec3(c2[0],c2[1],c2[2]))

class MyGfxObj(gfx.GfxObj):
  def __init__(self,name):
    gfx.GfxObj.__init__(self,name)
    self.rendered=False
    
  def CustomRenderGL(self,render_pass):
    self.rendered=True

class TestGfx(unittest.TestCase):
  def runTest(self):
    self.test_gradient()
    self.test_color()
    self.test_primlist()
    self.test_entity_reset()
    self.test_custom_gfx_obj()
    self.test_gfxobj_conv()

  def test_gfxobj_conv(self):
    e=mol.CreateEntity()
    gfx.Scene().Add(gfx.Entity("foo2",e))
    gfx.Scene()["foo2"].SetColor(gfx.YELLOW)
    
  def test_custom_gfx_obj(self):
    myobj=MyGfxObj("foo")
    gfx.Scene().Add(myobj)
    #self.assertTrue(myobj.rendered)

  def test_entity_reset(self):
    eh=mol.CreateEntity()
    eh2=mol.CreateEntity()
    go=gfx.Entity("foo",eh)
    go.Reset(eh)
    go.Reset("rnum=1")
    go.Reset(entity=eh2)
    go.Reset(eh.Select("rnum=2"))
    go.Reset(eh,"rnum=3",mol.MATCH_RESIDUES)
    go.Reset("rnum=4",eh2)

  def test_gradient(self):
    gs=[gfx.Gradient(),
        gfx.Gradient({0.0: [1,0,0], 1.0: gfx.Color(0,1,0)}),
        gfx.Gradient([[1,0,0], gfx.Color(0,1,0)])]
    gs[0].SetColorAt(0.0,gfx.Color(1.0,0.0,0.0))
    gs[0].SetColorAt(1.0,gfx.Color(0.0,1.0,0.0))
    for g in gs:
      self.assertAlmostEqual(col_delta(g.GetColorAt(0.0),gfx.Color(1.0,0,0)),0.0)
      self.assertAlmostEqual(col_delta(g.GetColorAt(0.5),gfx.Color(0.5,0.5,0)),0.0)
      self.assertAlmostEqual(col_delta(g.GetColorAt(1.0),gfx.Color(0,1.0,0)),0.0)

  def test_color(self):
    c=gfx.Color(0.5,0.3,0.2)
    self.assertAlmostEqual(c.r,0.5)
    self.assertAlmostEqual(c.g,0.3)
    self.assertAlmostEqual(c.b,0.2)
    self.assertAlmostEqual(c.a,1.0)
    self.assertAlmostEqual(c.red,0.5)
    self.assertAlmostEqual(c.green,0.3)
    self.assertAlmostEqual(c.blue,0.2)
    self.assertAlmostEqual(c.alpha,1.0)
    self.assertAlmostEqual(c[0],0.5)
    self.assertAlmostEqual(c[1],0.3)
    self.assertAlmostEqual(c[2],0.2)
    self.assertAlmostEqual(c[3],1.0)
    c.r=0.9
    self.assertAlmostEqual(c.r,0.9)
    self.assertAlmostEqual(c.red,0.9)
    self.assertAlmostEqual(c[0],0.9)

  def test_primlist(self):
    pl=gfx.PrimList("foo")
    pl.AddPoint([0,0,0])
    pl.AddPoint(geom.Vec3(1,2,3),color=gfx.RED)
    pl.AddLine([0,0,0],[1,2,3])
    pl.AddLine(geom.Vec3(0,0,0),geom.Vec3(1,2,3),color=gfx.BLUE)
    pl.AddSphere([0,0,0],radius=2.0)
    pl.AddSphere(geom.Vec3(1,2,3),color=gfx.RED,radius=3.0)
    pl.AddCyl([0,0,0],[1,2,3],radius=0.5,color=gfx.YELLOW)
    pl.AddCyl(geom.Vec3(0,0,0),geom.Vec3(1,2,3),radius1=0.5,radius2=0.1,color1=gfx.BLUE,color2=gfx.GREEN)
    pl.AddText("foo",[0,2,3])
    pl.AddText("bar",[-2,0,0],color=gfx.WHITE,point_size=8)
    if has_numpy:
      pl.AddMesh(numpy.zeros((5,3),dtype=numpy.float32),
                 numpy.zeros((5,3),dtype=numpy.float32),
                 numpy.zeros((5,4),dtype=numpy.float32),
                 numpy.zeros((2,3),dtype=numpy.uint32))
      pl.AddMesh(numpy.zeros((7,3),dtype=numpy.float32),
                 None,
                 None,
                 numpy.zeros((4,3),dtype=numpy.uint32))

      # Passing wrong data type should fail
      with self.assertRaises(Exception):
        pl.AddMesh(numpy.zeros((5, 3), dtype=numpy.uint32),
                   numpy.zeros((5, 3), dtype=numpy.float32),
                   numpy.zeros((5, 4), dtype=numpy.float32),
                   numpy.zeros((2, 3), dtype=numpy.uint32))
      with self.assertRaises(Exception):
        pl.AddMesh(numpy.zeros((5, 3), dtype=numpy.float32),
                   numpy.zeros((5, 3), dtype=numpy.uint32),
                   numpy.zeros((5, 4), dtype=numpy.float32),
                   numpy.zeros((2, 3), dtype=numpy.uint32))
      with self.assertRaises(Exception):
        pl.AddMesh(numpy.zeros((5, 3), dtype=numpy.float32),
                   numpy.zeros((5, 3), dtype=numpy.float32),
                   numpy.zeros((5, 4), dtype=numpy.uint32),
                   numpy.zeros((2, 3), dtype=numpy.uint32))
      with self.assertRaises(Exception):
        pl.AddMesh(numpy.zeros((5, 3), dtype=numpy.float32),
                   numpy.zeros((5, 3), dtype=numpy.float32),
                   numpy.zeros((5, 4), dtype=numpy.float32),
                   numpy.zeros((2, 3), dtype=numpy.float32))
                 

if __name__== '__main__':
  from ost import testutils
  testutils.RunTests()

