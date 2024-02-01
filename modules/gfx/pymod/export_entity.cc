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
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/gfx/entity.hh>
using namespace ost;
using namespace ost::gfx;
#include <ost/export_helper/pair_to_tuple_conv.hh>
#include "color_by_def.hh"

namespace {

void color_by_chain_01(Entity* e)
{
  e->ColorByChain();
}

void color_by_chain_02(Entity* e, const String& selection)
{
  e->ColorByChain(selection);
}

void color_by_element_01(Entity* e)
{
  e->ColorByElement();
}

void color_by_element_02(Entity* e, const String& selection)
{
  e->ColorByElement(selection);
}

void color_by_01(Entity* e,
                 const String& prop, 
                 const Gradient& gradient,
                 float minv,float maxv,
                 mol::Prop::Level hint)
{
  e->ColorBy(prop,gradient,minv,maxv,hint);
}

void color_by_02(Entity* e,
                 const String& prop, 
                 const Gradient& gradient,
                 float minv,float maxv)
{
  e->ColorBy(prop,gradient,minv,maxv);
}

void color_by_03(Entity* e,
                 const String& prop, 
                 const Gradient& gradient,
                 mol::Prop::Level hint)
{
  e->ColorBy(prop,gradient,hint);
}

void color_by_04(Entity* e,
                 const String& prop, 
                 const Gradient& gradient)
{
  e->ColorBy(prop,gradient);
}

void color_by_05(Entity* e,
                 const String& prop, 
                 const Color& c1, const Color& c2,
                 float minv,float maxv,
                 mol::Prop::Level hint)
{
  e->ColorBy(prop,c1,c2,minv,maxv,hint);
}

void color_by_06(Entity* e,
                 const String& prop, 
                 const Color& c1, const Color& c2,
                 float minv,float maxv)
{
  e->ColorBy(prop,c1,c2,minv,maxv);
}

void color_by_07(Entity* e,
                 const String& prop, 
                 const Color& c1, const Color& c2,
                 mol::Prop::Level hint)
{
  e->ColorBy(prop,c1,c2,hint);
}

void color_by_08(Entity* e,
                 const String& prop, 
                 const Color& c1, const Color& c2)
{
  e->ColorBy(prop,c1,c2);
}

void color_by_09(Entity* e,
                 const String& prop,
                 const Gradient& gradient,
                 const String& selection)
{
  e->ColorBy(prop,gradient,selection);
}

void color_by_10(Entity* e,
                 const String& prop,
                 const Gradient& gradient,
                 float minv, float maxv,
                 bool clamp)
{
  e->ColorBy(prop,gradient,minv,maxv,clamp);
}

// temporary, see comment in gfx/entity.hh
void detail_color_by_02(Entity* e,
                        const String& prop, 
                        const Gradient& gradient,
                        float minv,float maxv)
{
  e->DetailColorBy(prop,gradient,minv,maxv);
}


void radius_by_01(Entity* e,
                  const String& prop, 
                  float rmin,float rmax,
                  float minv,float maxv,
                  mol::Prop::Level hint)
{
  e->RadiusBy(prop,rmin,rmax,minv,maxv,hint);
}

void radius_by_02(Entity* e,
                  const String& prop, 
                  float rmin,float rmax,
                  float minv,float maxv)
{
  e->RadiusBy(prop,rmin,rmax,minv,maxv);
}

void radius_by_03(Entity* e,
                  const String& prop, 
                  float rmin,float rmax,
                  mol::Prop::Level hint)
{
  e->RadiusBy(prop,rmin,rmax,hint);
}

void radius_by_04(Entity* e,
                  const String& prop, 
                  float rmin,float rmax)
{
  e->RadiusBy(prop,rmin,rmax);
}


void ent_set_color1(Entity* e, const Color& c) {
  e->SetColor(c);
}
void ent_set_color2(Entity* e, const Color& c, const String& s) {
  e->SetColor(c,s);
}


void ent_apply_11(Entity* e, UniformColorOp& uco, bool store){
  e->Apply(uco,store);
}
void ent_apply_12(Entity* e, UniformColorOp& uco){
  e->Apply(uco);
}

void ent_apply_21(Entity* e, ByElementColorOp& beco, bool store){
  e->Apply(beco,store);
}
void ent_apply_22(Entity* e, ByElementColorOp& beco){
  e->Apply(beco);
}

void ent_apply_31(Entity* e, ByChainColorOp& beco, bool store){
  e->Apply(beco,store);
}
void ent_apply_32(Entity* e, ByChainColorOp& beco){
  e->Apply(beco);
}

void ent_apply_41(Entity* e, EntityViewColorOp& evco, bool store){
  e->Apply(evco,store);
}
void ent_apply_42(Entity* e, EntityViewColorOp& evco){
  e->Apply(evco);
}

void ent_apply_51(Entity* e, GradientLevelColorOp& glco, bool store){
  e->Apply(glco,store);
}
void ent_apply_52(Entity* e, GradientLevelColorOp& glco){
  e->Apply(glco);
}

void ent_apply_61(Entity* e, MapHandleColorOp& mhco, bool store){
  e->Apply(mhco,store);
}
void ent_apply_62(Entity* e, MapHandleColorOp& mhco){
  e->Apply(mhco);
}

RenderOptionsPtr ent_sline_opts(Entity* ent)
{
  return ent->GetOptions(RenderMode::SLINE);
}

void (Entity::*set_rm1)(RenderMode::Type, const mol::EntityView&, bool)=&Entity::SetRenderMode;
void (Entity::*set_rm3)(RenderMode::Type, const String&, bool)=&Entity::SetRenderMode;
void (Entity::*set_rm2)(RenderMode::Type)=&Entity::SetRenderMode;

void (Entity::*set_vis1)(const mol::EntityView&, bool)=&Entity::SetVisible;
void (Entity::*set_vis2)(const String&, bool)=&Entity::SetVisible;

RenderOptionsPtr ent_trace_opts(Entity* ent)
{
  return ent->GetOptions(RenderMode::TRACE);
}

RenderOptionsPtr ent_simple_opts(Entity* ent)
{
  return ent->GetOptions(RenderMode::SIMPLE);
}

RenderOptionsPtr ent_custom_opts(Entity* ent)
{
  return ent->GetOptions(RenderMode::CUSTOM);
}

RenderOptionsPtr ent_tube_opts(Entity* ent)
{
  return ent->GetOptions(RenderMode::TUBE);
}

RenderOptionsPtr ent_hsc_opts(Entity* ent)
{
  return ent->GetOptions(RenderMode::HSC);
}

RenderOptionsPtr ent_cpk_opts(Entity* ent)
{
  return ent->GetOptions(RenderMode::CPK);
}

void set_query(Entity* e, object o)
{
  LOG_WARNING("SetQuery is deprecated, use source property instead");

  if(o==object()) {
    e->SetQuery(mol::Query());
    return;
  }

  extract<String> str(o);
  if(str.check()) {
    e->SetQuery(mol::Query(str()));
    return;
  }

  extract<mol::Query> qry(o);
  if(qry.check()) {
    e->SetQuery(qry());
    return;
  }

  throw Error("expected string or mol::Query as parameter");
}

RenderOptionsPtr ent_ltrace_opts(Entity* ent)
{
  return ent->GetOptions(RenderMode::LINE_TRACE);
}

void set_selection(Entity* ent, object sel)
{
  object none;
  if (sel==none) {
    ent->SetSelection(ent->GetView().CreateEmptyView());
    return;
  }
  try {
    String sel_string=extract<String>(sel);
    ent->SetSelection(ent->GetView().Select(sel_string));
  } catch (error_already_set& e) {
    PyErr_Clear();
    mol::EntityView view=extract<mol::EntityView>(sel);
    ent->SetSelection(view);
  }
}

}

void export_Entity()
{
  void (Entity::*reset1)(const mol::EntityHandle&) = &Entity::Reset;
  void (Entity::*reset2)(const mol::EntityHandle&, const mol::Query&) = &Entity::Reset;
  void (Entity::*reset3)(const mol::EntityHandle&, const mol::Query&, mol::QueryFlags) = &Entity::Reset;
  void (Entity::*reset4)(const mol::EntityView&) = &Entity::Reset;

  class_<Entity, boost::shared_ptr<Entity>, bases<GfxObj>, boost::noncopyable>("Entity", init<const String&, const mol:: EntityHandle&, optional<const mol:: Query&, mol::QueryFlags> >())
    .def(init<const String&, RenderMode::Type, const mol::EntityHandle&, optional<const mol::Query&, mol::QueryFlags> >())
    .def(init<const String&, const mol::EntityView&>())
    .def(init<const String&, RenderMode::Type, const mol::EntityView&>())
    .def("_reset1",reset1)
    .def("_reset2",reset2)
    .def("_reset3",reset3)
    .def("_reset4",reset4)
    .def("SetColor",ent_set_color1)
    .def("SetColor",ent_set_color2)
    .def("SetDetailColor", &Entity::SetDetailColor, arg("sel")=String(""))
    .def("SetColorForAtom", &Entity::SetColorForAtom)
    .def("Rebuild", &Entity::Rebuild)
    .def("UpdatePositions",&Entity::UpdatePositions)
    .def("BlurSnapshot", &Entity::BlurSnapshot)
    .def("SetBlurFactors",&Entity::SetBlurFactors)
    .def("SetBlur",&Entity::SetBlur)
    .def("GetBoundingBox",&Entity::GetBoundingBox)
    .def("SetSelection",&Entity::SetSelection)
    .def("GetSelection",&Entity::GetSelection)    
    .add_property("selection", &Entity::GetSelection, 
                  &set_selection)
    .def("GetView", &Entity::GetView)
    .add_property("view",&Entity::GetView)
    .def("UpdateView", &Entity::UpdateView)
    .def("SetQuery", set_query)
    .def("SetQueryView",&Entity::SetQueryView)
    .def("GetQueryView",&Entity::GetQueryView)
    .add_property("query_view",&Entity::GetQueryView,&Entity::SetQueryView)
    .def("GetRenderModeName", &Entity::GetRenderModeName)
    .def("GetNotEmptyRenderModes", &Entity::GetNotEmptyRenderModes)
    .def("SetRenderMode", set_rm1, (arg("mode"), arg("view"), arg("keep")=false))
    .def("SetRenderMode", set_rm2)
    .def("SetRenderMode", set_rm3, (arg("mode"), arg("sel"), arg("keep")=false))    
    .def("SetEnableRenderMode", &Entity::SetEnableRenderMode)
    .def("IsRenderModeEnabled", &Entity::IsRenderModeEnabled)
    .def("SetVisible", set_vis1, (arg("view"), arg("flag")=true))
    .def("SetVisible", set_vis2, (arg("sel"), arg("flag")=true))    
    .def("ColorBy", color_by_10) // this line must be before color_by_01 because
                                 // of boost python overload resolution
    .def("ColorBy", color_by_01)
    .def("ColorBy", color_by_02)
    .def("ColorBy", color_by_03)
    .def("ColorBy", color_by_04)
    .def("ColorBy", color_by_05)
    .def("ColorBy", color_by_06)
    .def("ColorBy", color_by_07)
    .def("ColorBy", color_by_08)
    .def("ColorBy", color_by_09)
    .def("DetailColorBy", detail_color_by_02)
    COLOR_BY_DEF()
    .def("RadiusBy", radius_by_01)
    .def("RadiusBy", radius_by_02)
    .def("RadiusBy", radius_by_03)
    .def("RadiusBy", radius_by_04)
    .def("ResetRadiusBy", &Entity::ResetRadiusBy)
    .def("PickAtom", &Entity::PickAtom)
    .def("PickBond", &Entity::PickBond)
    .def("ColorByElement", color_by_element_01)
    .def("ColorByElement", color_by_element_02)
    .def("ColorByChain", color_by_chain_01)
    .def("ColorByChain", color_by_chain_02)
    .def("CleanColorOps", &Entity::CleanColorOps)
    .def("ReapplyColorOps", &Entity::ReapplyColorOps)
    .def("GetOptions", &Entity::GetOptions)
    .add_property("sline_options", &ent_sline_opts)
    .add_property("simple_options", &ent_simple_opts)    
    .add_property("tube_options", &ent_tube_opts)
    .add_property("custom_options", &ent_custom_opts)
    .add_property("cartoon_options", &ent_hsc_opts)    
    .add_property("cpk_options", &ent_cpk_opts)
    .add_property("trace_options", &ent_trace_opts)
    .add_property("line_trace_options", &ent_ltrace_opts)
    .def("ApplyOptions", &Entity::ApplyOptions)
    .def("SetOptions", &Entity::SetOptions)
    .def("Apply",&ent_apply_11)
    .def("Apply",&ent_apply_12)
    .def("Apply",&ent_apply_21)
    .def("Apply",&ent_apply_22)
    .def("Apply",&ent_apply_31)
    .def("Apply",&ent_apply_32)
    .def("Apply",&ent_apply_41)
    .def("Apply",&ent_apply_42)
    .def("Apply",&ent_apply_51)
    .def("Apply",&ent_apply_52)
    .def("Apply",&ent_apply_61)
    .def("Apply",&ent_apply_62)
    .add_property("seq_hack",&Entity::GetSeqHack,&Entity::SetSeqHack)
  ;
  //register_ptr_to_python<EntityP>();
  
  to_python_converter<std::pair<GfxObjP, mol::AtomHandle>, 
                      PairToTupleConverter<GfxObjP, mol::AtomHandle> >();

  class_<RenderModeTypes>("RenderModeTypes", init<>())
    .def(vector_indexing_suite<RenderModeTypes, true >());
}
