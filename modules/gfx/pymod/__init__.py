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
from ._ost_gfx import *
from .py_gfx_obj import PyGfxObj
import functools

WHITE=RGB(1.0,1.0,1.0)
BLACK=RGB(0.0,0.0,0.0)
GREY=RGB(0.5,0.5,0.5)
RED=RGB(1.0,0.0,0.0)
DARKRED=RGB(0.5,0.0,0.0)
LIGHTRED=RGB(1.0,0.5,0.5)
GREEN=RGB(0.0,1.0,0.0)
DARKGREEN=RGB(0.0,0.5,0.0)
LIGHTGREEN=RGB(0.5,1.0,0.5)
BLUE=RGB(0.0,0.0,1.0)
DARKBLUE=RGB(0.0,0.0,0.5)
LIGHTBLUE=RGB(0.5,0.5,1.0)
YELLOW=RGB(1.0,1.0,0.0)
DARKYELLOW=RGB(0.5,0.5,0.0)
LIGHTYELLOW=RGB(1.0,1.0,0.5)
CYAN=RGB(0.0,1.0,1.0)
DARKCYAN=RGB(0.0,0.5,0.5)
LIGHTCYAN=RGB(0.5,1.0,1.0)
MAGENTA=RGB(1.0,0.0,1.0)
DARKMAGENTA=RGB(0.5,0.0,0.5)
LIGHTMAGENTA=RGB(1.0,0.5,1.0)
PURPLE=MAGENTA
DARKPURPLE=DARKMAGENTA
LIGHTPURPLE=LIGHTMAGENTA
ORANGE=RGB(1.0,0.5,0.0)
DARKORANGE=RGB(0.5,0.25,0.0)
LIGHTORANGE=RGB(1.0,0.75,0.5)

def Stereo(mode,flip=None,alg=None):
  """
  Stereo control

  :param mode: 0=off, 1=quad-buffered, 2=interlaced
  :type mode: int
  :param flip: invert order of left/right display
  :type flip: bool
  :param alg: stereo algorithm (0 or 1)
  :type param: int
  """
  if(flip):
    Scene().SetStereoFlip(flip)
  if(alg):
    Scene().SetStereoAlg(alg)

  Scene().SetStereoMode(mode)

def FitToScreen(gfx_ent, width=None, height=None, margin=0.05):
  """
  Setup camera such that it is centered on the graphical entity and the entity 
  fits the entire viewport. The longest axes of the entity are aligned along 
  the x- and y- axes of the screen.
  
  :param gfx_ent: The graphical entity
  :type gfx_ent: str or :class:`Entity`
  
  
  """
  from ost import geom
  import math
  def _XYZ(view, axes):
    """
    returns the vectors in x, y and z direction respectively. The smallest 
    vector is in z, then y, and the largest along x.
    """
    rows=[axes.GetRow(i) for i in range(3)]
    lengths=[]
    for axe in rows:
      min_proj=geom.Dot(axe, view.atoms[0].pos)
      max_proj=min_proj
      for atom in view.atoms[1:]:
        proj=geom.Dot(axe, atom.pos)
        min_proj=min(proj, min_proj)
        max_proj=max(proj, max_proj)
      lengths.append(max_proj-min_proj)
    def cmp_x(rhs, lhs):
      # replaced cmp when porting to Python 3
      #return cmp(lhs[1], rhs[1])
      return (lhs[1] > rhs[1]) - (lhs[1] < rhs[1])
    sorted_axes=sorted(zip(rows, lengths), key=functools.cmp_to_key(cmp_x))
    return [r*l for r,l in sorted_axes]
  scene=Scene()
  if not isinstance(gfx_ent, Entity):
    gfx_ent=scene[str(gfx_ent)]
  width=width and width or scene.viewport.width
  height=height and height or scene.viewport.height
  atom_positions=geom.Vec3List([atom.pos for atom in gfx_ent.view.atoms])
  axes=atom_positions.principal_axes
  sorted_axes=_XYZ(gfx_ent.view, axes)
  x_bigger_than_y=geom.Length(sorted_axes[0])>geom.Length(sorted_axes[1])
  if x_bigger_than_y:
    if width>height:    
      x_axes=geom.Normalize(sorted_axes[0])
      y_axes=geom.Normalize(sorted_axes[1])
    else:
      x_axes=geom.Normalize(sorted_axes[1])
      y_axes=geom.Normalize(sorted_axes[0])
  else:
    if width>height:    
      x_axes=geom.Normalize(sorted_axes[1])
      y_axes=geom.Normalize(sorted_axes[0])
    else:
      x_axes=geom.Normalize(sorted_axes[0])
      y_axes=geom.Normalize(sorted_axes[1])
  z_axes=geom.Normalize(geom.Cross(x_axes, y_axes))
  rotation=geom.Mat3(x_axes[0], x_axes[1], x_axes[2],
                     y_axes[0], y_axes[1], y_axes[2],
                     z_axes[0], z_axes[1], z_axes[2])
  rtc=geom.Mat4(rotation)
  
  center=gfx_ent.center
  aspect=float(width)/float(height)
  factor_y=1.0/math.tan(math.radians(scene.fov))
  factor_x=factor_y/aspect
  z_off=geom.Length(sorted_axes[2])*0.5
  rtc[0,3]=center[0]
  rtc[1,3]=center[1]
  rtc[2,3]=center[2]
  rtc[3,0]=0
  rtc[3,1]=0
  rtc[3,2]=-(max(factor_x*(1+margin)*geom.Length(sorted_axes[0]), 
                factor_y*(1+margin)*geom.Length(sorted_axes[1]))+z_off)
  scene.SetRTC(rtc)


class GfxNodeListAttrProxy:
  def __init__(self, node_list, name):
    self._node_list=node_list
    self._name=name

  def __iter__(self):
    for node in self._node_list:
      yield getattr(node, self._name)

  def __call__(self, *args, **kwargs):
    for node in self._node_list:
      bound_method=getattr(node, self._name)
      bound_method(*args, **kwargs)

class GfxNodeListProxy(object):
  def __init__(self, node_list):
    self._nodes=node_list

  def __getattr__(self, name):
    if name.startswith('_'):
      return super(GfxNodeListProxy, self).__getattr__(name)
    return GfxNodeListAttrProxy(self._nodes, name)

  def __setattr__(self, name, value):
    if name.startswith('_'):
      super(GfxNodeListProxy, self).__setattr__(name, value)
    for node in self._nodes:
      setattr(node, name, value)

def _Match(scene, pattern="*"):
  import os
  import fnmatch
  def _Recurse(path, node, pattern):
    matches=[]
    for child in node.children:
      full_name=os.path.join(path, child.name)
      if fnmatch.fnmatchcase(full_name, pattern):
        matches.append(child)
      matches.extend(_Recurse(full_name, child, pattern))
    return matches
  return GfxNodeListProxy(_Recurse("", Scene().root_node, pattern))

SceneSingleton.Match=_Match

def _to_vec3(p):
  import ost.geom
  if isinstance(p,ost.geom.Vec3):
    return p
  else:
    try:
      return ost.geom.Vec3(p[0],p[1],p[2])
    except:
      raise TypeError("expected either a sequence or a geom.Vec3 object")
  

def _primlist_add_point(self,pos,color=None):
  pos=_to_vec3(pos)
  if not color:
    color=WHITE
  self._add_point(pos,color)
  
def _primlist_add_line(self,pos1,pos2,color1=None,color2=None,color=None):
  pos1=_to_vec3(pos1)
  pos2=_to_vec3(pos2)
  if not color:
    color=WHITE
  if not color1:
    color1=color
  if not color2:
    color2=color1
  self._add_line(pos1,pos2,color1,color2)

def _primlist_add_sphere(self,cen,radius=1.0,color=None):
  pos=_to_vec3(cen)
  if not color:
    color=WHITE
  self._add_sphere(pos,radius,color)
  
def _primlist_add_cyl(self,pos1,pos2,radius1=None,radius2=None,radius=None,color1=None,color2=None,color=None,):
  pos1=_to_vec3(pos1)
  pos2=_to_vec3(pos2)
  if radius is None:
    radius=1.0
  if radius1 is None:
    radius1=radius
  if radius2 is None:
    radius2=radius1
  if not color:
    color=WHITE
  if not color1:
    color1=color
  if not color2:
    color2=color1
  self._add_cyl(pos1,pos2,radius1,radius2,color1,color2)

def _primlist_add_text(self,text,pos,color=None,point_size=None):
  pos=_to_vec3(pos)
  if not color:
    color=WHITE
  if not point_size:
    point_size=1.0
  self._add_text(text,pos,color,point_size)

PrimList.AddPoint=_primlist_add_point
PrimList.AddLine=_primlist_add_line
PrimList.AddSphere=_primlist_add_sphere
PrimList.AddCyl=_primlist_add_cyl
PrimList.AddText=_primlist_add_text

# entity reset

def _entity_reset(self,*args,**kwargs):
  import ost.mol as mol
  eh=None
  ev=None
  qr=None
  qf=None
  for a in args:
    if isinstance(a,mol.Query):
      if qr:
        raise TypeError("Reset: more than one query string given")
      qr=a
    elif isinstance(a,mol.EntityHandle):
      if eh:
        raise TypeError("Reset: more than one entity handle given")
      eh=a
    elif isinstance(a,mol.EntityView):
      if ev:
        raise TypeError("Reset: more than one entity view given")
      ev=a
    elif isinstance(a,str):
      if qr:
        raise TypeError("Reset: more than one query string given")
      qr=mol.Query(a)
    elif isinstance(a,int):
      if qf:
        raise TypeError("Reset: more than one QueryFlags given")
      qf=a
    else:
      raise TypeError("Reset: unknown option of type '%s' given"%type(a))

  for key,val in kwargs.items():
    if key=="entity":
      if not isinstance(val,mol.EntityHandle):
        raise TypeError("Reset: expected mol.EntityHandle for 'entity' option")
      if eh:
        raise TypeError("Reset: more than one entity handle given")
      eh=val
    elif key=="view":
      if not isinstance(val,mol.EntityView):
        raise TypeError("Reset: expected mol.EntityView for 'view' option")
      if ev:
        raise TypeError("Reset: more than one entity view given")
      ev=val
    elif key=="query":
      if isinstance(val,mol.Query):
        pass
      elif isinstance(val,str):
        val=mol.Query(val)
      else:
        raise TypeError("Reset: expected mol.Query or string for 'query' option")
      if qr:
        raise TypeError("Reset: more than one query string given")
      qr=val
    elif key=="query_flags":
      if not isinstance(val,int):
        raise TypeError("Reset: expected integer for 'query_flags' option")
      if qf:
        raise TypeError("Reset: more than one query flags given")
      qf=val
    else:
      raise TypeError("Reset: unknown key '%s'"%key)
        
  if eh and ev:
    raise TypeError("Reset: entity and view are mutually exclusive options")

  if ev:
    self._reset4(ev)
  else:
    if not eh:
      eh = self.query_view.entity
    if not qr:
      qr = self.query_view.query
    if not qf:
      qf = self.query_view.GetFlags()
    self._reset3(eh,qr,qf)

Entity.Reset=_entity_reset

def _scene_export(self,*args,**kwargs):
  """
  scene.Export(Exporter)
  scene.Export("file.png")
  scene.Export("file.png",(width,height),samples=0,transparency=False)
  deprecated:
  scene.Export("file.png",width,height,samples=0,transparency=False)
  scene.Export("file.png",width,height,transparency) 
  """
  scene=Scene()
  tp=False
  sa=0
  if "tp" in kwargs:
    tp=int(kwargs["tp"])
  if "transparency" in kwargs:
    tp=int(kwargs["transparency"])
  if "samples" in kwargs:
    sa=int(kwargs["samples"])

  if len(args)==1:
    if isinstance(args[0],Exporter):
      scene._export_via_exporter(args[0])
      return
    elif type(args[0])==type(""):
      scene._export_screen(args[0],tp)
      return
  elif len(args)==2:
    if type(args[0]==type("")):
      # assume second argument is a dimension
      width=int(args[1][0])
      height=int(args[1][1])
      scene._export_buffer(args[0],width,height,sa,tp)
      return
  elif len(args)==3:
    if type(args[0]==type("")):
      width=int(args[1])
      height=int(args[2])
      scene._export_buffer(args[0],width,height,sa,tp)
      return
  elif len(args)==4:
    if type(args[0]==type("")):
      width=int(args[1])
      height=int(args[2])
      tp=int(args[3])
      scene._export_buffer(args[0],width,height,sa,tp)
      return
  # getting here indicates an error
  raise RuntimeError("""invalid arguments to scene.Export; expected one of
  Export(gfx.Exporter)
  Export('file.png')
  Export('file.png',(width,height),samples=0, transparency=False)
  Export('file.png',width,height,samples=0, transparency=False) -> deprecated
  Export('file.png',width,height,transparency) -> deprecated
  """)

SceneSingleton.Export=_scene_export

import __main__ as main_mod
main_mod.scene=Scene()
main_mod.scene.Stereo=Stereo

import ost as ost_mod
ost_mod.scene=Scene()
ost_mod.scene.Stereo=Stereo

def GostExporter(file,scale=1.0,to_origin=True):
  e=GostExporter_(file)
  e.scale=scale
  e.to_origin=to_origin
  return e

def ColladaExporter(file,scale=1.0,to_origin=True):
  e=ColladaExporter_(file)
  e.scale=scale
  e.to_origin=to_origin
  return e

def _go_get_vis(go):
  return go.IsVisible()

def _go_set_vis(go,flag):
  if flag:
    go.Show()
  else:
    go.Hide()

GfxObj.visible=property(_go_get_vis,_go_set_vis)


