//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------

#include <ost/mol/entity_property_mapper.hh>

#include "mapped_property.hh"

namespace ost { namespace gfx { namespace impl {

namespace {

float InterpolatedValue(const img::ImageHandle& mh, const geom::Vec3& xyz)
{
  geom::Vec3 uvw=mh.CoordToIndex(xyz);
  float x=fmod(float(uvw[0]), float(1.0));
  float y=fmod(float(uvw[1]), float(1.0));
  float z=fmod(float(uvw[2]), float(1.0));
  uint u=static_cast<uint>(uvw[0]);
  uint v=static_cast<uint>(uvw[1]);
  uint w=static_cast<uint>(uvw[2]);

  float value=0.0;
  if(u<mh.GetSize()[0]-1 && v<mh.GetSize()[1]-1 && w<mh.GetSize()[2]-1) {
   value=mh.GetReal(ost::img::Point(u,v,w))*(1.0-x)*(1.0-y)*(1.0-z)+
     mh.GetReal(ost::img::Point(u+1,v,w))*(x)*(1.0-y)*(1.0-z)+
     mh.GetReal(ost::img::Point(u,v+1,w))*(1.0-x)*(y)*(1.0-z)+
     mh.GetReal(ost::img::Point(u,v,w+1))*(1.0-x)*(1.0-y)*(z)+
     mh.GetReal(ost::img::Point(u+1,v,w+1))*(x)*(1.0-y)*(z)+
     mh.GetReal(ost::img::Point(u,v+1,w+1))*(1.0-x)*(y)*(z)+
     mh.GetReal(ost::img::Point(u+1,v+1,w))*(x)*(y)*(1.0-z)+
     mh.GetReal(ost::img::Point(u+1,v+1,w+1))*(x)*(y)*(z);
  }
  return static_cast<float>(value);
}  
}

float DLLEXPORT_OST_GFX MappedProperty(const mol::EntityView& ev, 
                                       const String& prop,
                                       const geom::Vec3& pos)
{
  mol::EntityPropertyMapper epm(prop, mol::Prop::UNSPECIFIED);
  mol::AtomViewList avl = ev.FindWithin(pos,5.0);
  if(avl.empty()) return 0.0;

  mol::AtomView cla = avl[0];
  float dist2 = Length2(cla.GetPos()-pos);
  for(uint i=1;i<avl.size();++i) {
    float d2 = geom::Length2(avl[i].GetPos()-pos);
    if(d2<dist2) {
      cla=avl[i];
      dist2=d2;
    }
  }
  return epm.Get(cla,0.0);
}

Color DLLEXPORT_OST_GFX MappedProperty(const mol::EntityView& ev, 
                                       const String& prop,
                                       const Gradient& g, float minv, float maxv,
                                       const geom::Vec3& pos)
{
  Color nrvo;
  nrvo=g.GetColorAt(Normalize(MappedProperty(ev,prop,pos), minv, maxv));
  return nrvo;
}
                        
float DLLEXPORT_OST_GFX MappedProperty(const img::MapHandle& mh, 
                                       const String& prop,
                                       const geom::Vec3& pos)
{
  return InterpolatedValue(mh,pos);
}

Color DLLEXPORT_OST_GFX MappedProperty(const img::MapHandle& mh, 
                                       const String& prop,
                                       const Gradient& g, float minv, float maxv,
                                       const geom::Vec3& pos)
{
  Color nrvo;
  nrvo=g.GetColorAt(Normalize(MappedProperty(mh,prop,pos), minv, maxv));
  return nrvo;
}

}}}
