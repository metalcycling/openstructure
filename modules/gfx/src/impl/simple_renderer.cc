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

/*
  Authors: Ansgar Philippsen, Stefan Scheuber, Marco Biasini
 */

#include <ost/gfx/gl_helper.hh>
#include <ost/gfx/scene.hh>
#include <ost/gfx/entity.hh>
#include <ost/gfx/povray.hh>
#include "simple_renderer.hh"

namespace ost { namespace gfx { namespace impl {


  
SimpleRenderer::SimpleRenderer(): options_(new SimpleRenderOptions()) 
{
  this->SetName("Fast Bonds");
  this->SetFixedPickRadius(0.2);
}

bool SimpleRenderer::CanSetOptions(RenderOptionsPtr& render_options)
{
 return render_options.get()->GetRenderMode()==RenderMode::SIMPLE;
}

void SimpleRenderer::SetOptions(RenderOptionsPtr& render_options)
{
  options_ = boost::static_pointer_cast<SimpleRenderOptions>(render_options);
}

RenderOptionsPtr SimpleRenderer::GetOptions()
{
 return options_;
}

void SimpleRenderer::PrepareRendering()
{
  ConnectRendererBase::PrepareRendering();
  va_.Clear();
  this->PrepareRendering(view_, va_);
  sel_va_.Clear();  
  if (this->HasSelection()) {
    this->PrepareRendering(sel_view_, sel_va_);
    sel_va_.SetLineWidth(options_->GetLineWidth()+5.0);
    sel_va_.SetColorMaterial(false);
  }
  sel_state_=0;
  state_=0;
}

geom::Vec3 SimpleRenderer::GetDoubleBondPlane(mol::BondHandle b)
{
  int n=0;
  geom::Vec3 vec;
  mol::AtomHandle a1=b.GetFirst();
  mol::AtomHandle a2=b.GetSecond();
  mol::AtomHandleList partner1=a1.GetBondPartners();
  mol::AtomHandleList partner2=a2.GetBondPartners();
  geom::Vec3 bond_vec=a2.GetPos()-a1.GetPos();
  GetBondPartnerNormal(vec, n, bond_vec, partner1, a1, a2);
  GetBondPartnerNormal(vec, n, bond_vec, partner2, a1, a2);
  if(n==0) {
    //set vec into x direction if no partners are present
    vec.x=1;
  }
  vec=geom::Normalize(geom::Cross(vec,bond_vec));
  return vec;
}

void SimpleRenderer::GetBondPartnerNormal(geom::Vec3& vec,
                                           int& n,
                                           geom::Vec3& bond_vec,
                                           mol::AtomHandleList& ahl,
                                           mol::AtomHandle& a1,
                                           mol::AtomHandle& a2)
{
  // find plane by calculating average cross product between double
  // bond vector and all neighboring bond vectors
  for (mol::AtomHandleList::const_iterator i=ahl.begin();i!=ahl.end();++i) {
    if((*i)!=a1 && (*i)!=a2) {
      geom::Vec3 vec_=geom::Normalize(geom::Cross(
                                              i->GetPos()-a1.GetPos(),
                                              bond_vec));
      if(geom::Dot(vec,vec_)<0) vec_*=-1;
      vec+=vec_;
      n++;
    }
  }
}

void SimpleRenderer::PrepareRendering(GfxView& view, IndexedVertexArray& va)
{
  if (options_!=NULL) {
    static const float cross_size=0.25;
    va.SetLighting(false);
    va.SetCullFace(false);
    va.SetColorMaterial(false);
    va.SetLineWidth(options_->GetLineWidth());
    va.SetPointSize(options_->GetLineWidth());
    va.SetMode(0x2);
    va.SetTwoSided(true);
    va.SetAALines(options_->GetAALines());

    for (BondEntryList::const_iterator it=view.bond_list.begin();
         it!=view.bond_list.end();++it) {
      const geom::Vec3& p0=it->bond.GetFirst().GetPos();
      const geom::Vec3& p2=it->bond.GetSecond().GetPos();
      geom::Vec3 p1=(p0+p2)*0.5;

      if(it->bond.GetBondOrder()!=2 || !options_->GetBondOrderFlag()) {
        VertexID id0 = va.Add(p0,geom::Vec3(),it->atom1->color);
        VertexID id1 = va.Add(p1,geom::Vec3(),it->atom1->color);
        VertexID id2 = va.Add(p1,geom::Vec3(),it->atom2->color);
        VertexID id3 = va.Add(p2,geom::Vec3(),it->atom2->color);
        va.AddLine(id0,id1);
        va.AddLine(id2,id3);
      }
      if(it->bond.GetBondOrder()>=2 && options_->GetBondOrderFlag()) {
        //draw double bond lines for all higher order bonds
        geom::Vec3 vec;
        if(it->bond.GetBondOrder()==2) {
          vec=GetDoubleBondPlane(it->bond);
        } else {
          vec.x=1;
        }
        Real dbd=options_->GetBondOrderDistance();
        VertexID id00 = va.Add(p0+dbd*vec,geom::Vec3(),it->atom1->color);
        VertexID id01 = va.Add(p0-dbd*vec,geom::Vec3(),it->atom1->color);
        VertexID id10 = va.Add(p1+dbd*vec,geom::Vec3(),it->atom1->color);
        VertexID id11 = va.Add(p1-dbd*vec,geom::Vec3(),it->atom1->color);
        VertexID id20 = va.Add(p1+dbd*vec,geom::Vec3(),it->atom2->color);
        VertexID id21 = va.Add(p1-dbd*vec,geom::Vec3(),it->atom2->color);
        VertexID id30 = va.Add(p2+dbd*vec,geom::Vec3(),it->atom2->color);
        VertexID id31 = va.Add(p2-dbd*vec,geom::Vec3(),it->atom2->color);
        va.AddLine(id00,id10);
        va.AddLine(id20,id30);
        va.AddLine(id01,id11);
        va.AddLine(id21,id31);
      }
    }
    for (std::vector<long>::const_iterator oit=view.orphan_atom_list.begin();
         oit!=view.orphan_atom_list.end();++oit) {
      AtomEntryMap::const_iterator it = view.atom_map.find(*oit);
      if(it==view.atom_map.end()) continue;
      geom::Vec3 p0 = it->second.atom.GetPos();
      VertexID id0 = va.Add(p0-geom::Vec3(cross_size,0.0,0.0),
                             geom::Vec3(),it->second.color);
      VertexID id1 = va.Add(p0+geom::Vec3(cross_size,0.0,0.0),
                             geom::Vec3(),it->second.color);
      va.AddLine(id0,id1);
      id0 = va.Add(p0-geom::Vec3(0.0,cross_size,0.0),
                    geom::Vec3(),it->second.color);
      id1 = va.Add(p0+geom::Vec3(0.0,cross_size,0.0),
                    geom::Vec3(),it->second.color);
      va.AddLine(id0,id1);
      id0 = va.Add(p0-geom::Vec3(0.0,0.0,cross_size),
                    geom::Vec3(),it->second.color);
      id1 = va.Add(p0+geom::Vec3(0.0,0.0,cross_size),
                    geom::Vec3(),it->second.color);
      va.AddLine(id0,id1);
    }
  }  
}

BondEntryList& SimpleRenderer::GetBondEntryList()
{
  return view_.bond_list;
}

void SimpleRenderer::RenderPov(PovState& pov, const std::string& name)
{
  if(view_.atom_map.empty() && view_.bond_list.empty()) return;
  pov.write_merge_or_union(name);
  
  for (AtomEntryMap::const_iterator it=view_.atom_map.begin();
       it!=view_.atom_map.end();++it) {
    pov.write_sphere(it->second.atom.GetPos(),
                     options_->GetLineWidth()*0.05,
                     it->second.color,
                     name);
  }

  for(BondEntryList::const_iterator it=view_.bond_list.begin();it!=view_.bond_list.end();++it) {
    const geom::Vec3& p0=it->bond.GetFirst().GetPos();
    const geom::Vec3& p2=it->bond.GetSecond().GetPos();
    geom::Vec3 p1=(p0+p2)*0.5;
    pov.write_cyl(p0,p1,options_->GetLineWidth()*0.05,
                  it->atom1->color,name,true);
    pov.write_cyl(p1,p2,options_->GetLineWidth()*0.05,
                  it->atom2->color,name,true);
  }
  pov.inc() << " }\n";
}

}}}
