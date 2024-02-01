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
  Author: Stefan Scheuber, Marco Biasini
*/

#include <ost/gfx/gl_helper.hh>

#include <ost/mol/view_op.hh>

#include <ost/gfx/input.hh>
#include <ost/gfx/gfx_object.hh>
#include <ost/gfx/scene.hh>

#include "entity_renderer.hh"

#if OST_SHADER_SUPPORT_ENABLED
#include <ost/gfx/shader.hh>
#endif //OST_SHADER_SUPPORT_ENABLED

namespace ost { namespace gfx { namespace impl {

EntityRenderer::EntityRenderer():
  name_(""),
  enabled_(true)
{}

void EntityRenderer::FlagPositionsDirty()
{
  state_|=DIRTY_VA;
  sel_state_|=DIRTY_VA;
}

 void EntityRenderer::Debug(unsigned int flags)
 {
   debug_flags_=flags;
   va_.DrawNormals(debug_flags_&0x1);
 }

const String& EntityRenderer::GetName() const{
  return name_;
}

void EntityRenderer::SetName(const String& name)
{
  name_=name;
}

void EntityRenderer::SetEnabled(bool enabled){
  enabled_=enabled;
}

bool EntityRenderer::IsEnabled() const
{
  return enabled_;
}

const Color& EntityRenderer::GetSelectionColor() const
{
  static Color selection_color(0.3,0.9,0.1,0.6);  
  return selection_color;
}

const Color& EntityRenderer::GetSelectionColorOutline() const
{
#if OST_SHADER_SUPPORT_ENABLED
  static Color selection_color_outline(0.3,0.9,0.1,1.0);
#else
  static Color selection_color_outline=this->GetSelectionColor();
#endif
  return selection_color_outline;
}

bool EntityRenderer::IsDirty() const
{
  return sel_state_ || state_;
}

void EntityRenderer::AddView(const mol::EntityView& view)
{
  full_view_=full_view_.IsValid() ? Union(full_view_, view) : view.Copy();
  effective_view_=hidden_view_.IsValid() ? mol::Difference(full_view_,hidden_view_): full_view_;

  state_=DIRTY_VIEW|DIRTY_VA;
  sel_state_=DIRTY_VIEW|DIRTY_VA;
}

void EntityRenderer::SubstractView(const mol::EntityView& view)
{
  if (full_view_.IsValid()) {
    full_view_=mol::Difference(full_view_, view);
    effective_view_=hidden_view_.IsValid() ? mol::Difference(full_view_,hidden_view_): full_view_;
  }
  state_=DIRTY_VIEW|DIRTY_VA;
  sel_state_=DIRTY_VIEW|DIRTY_VA;
}

void EntityRenderer::ClearViews()
{
  if (full_view_.IsValid() && full_view_.HasAtoms()) {
    full_view_=full_view_.CreateEmptyView();
    if(effective_view_.IsValid()){
      effective_view_=effective_view_.CreateEmptyView();
    }
    full_sel_=full_view_.CreateEmptyView();
    sel_=full_view_.CreateEmptyView();
    state_=DIRTY_VIEW|DIRTY_VA;
    sel_state_=DIRTY_VIEW|DIRTY_VA;    
  }
}

bool EntityRenderer::HasDataToRender() const
{
  return effective_view_ && effective_view_.HasAtoms();
}

mol::EntityView EntityRenderer::GetFullView()
{
  return full_view_;
}

mol::EntityView EntityRenderer::GetEffectiveView()
{
  return effective_view_;
}

void EntityRenderer::Render(RenderPass pass)
{
  assert(sel_state_==0 && state_==0);
  if (pass==STANDARD_RENDER_PASS) {
    va_.DrawNormals(debug_flags_&0x1);
    va_.RenderGL();
  } else if (pass==GLOW_RENDER_PASS && this->HasSelection()) {
    sel_va_.RenderGL();
  }
}

void EntityRenderer::RenderPov(PovState& pov, const std::string& name)
{
  va_.RenderPov(pov,name);
}

void EntityRenderer::Export(Exporter* ex)
{
  va_.Export(ex);
}

bool EntityRenderer::HasSelection() const
{
  return (sel_.IsValid() && sel_.HasAtoms());
}

void EntityRenderer::SetSelection(const mol::EntityView& sel)
{
  if(sel.IsValid()){
    full_sel_=sel;
    sel_=(hidden_view_.IsValid())? mol::Difference(full_sel_,hidden_view_) : full_sel_;
    sel_state_=DIRTY_VIEW|DIRTY_VA;
  }
}

void EntityRenderer::SetVisible(const mol::EntityView& view, bool visible)
{
  if(visible){
    if(hidden_view_.IsValid()){
      hidden_view_=mol::Difference(hidden_view_, view);
    }
  }
  else{
    hidden_view_=hidden_view_.IsValid() ? Union(hidden_view_, view) : view.Copy();
  }

  if(sel_.IsValid()){
    sel_= hidden_view_.IsValid() ? mol::Difference(full_sel_,hidden_view_) : sel_;
  }

  if(full_view_.IsValid()){
    effective_view_=hidden_view_.IsValid() ? mol::Difference(full_view_,hidden_view_): full_view_;
  }
  state_=DIRTY_VIEW|DIRTY_VA;
  sel_state_=DIRTY_VIEW|DIRTY_VA;
}

void EntityRenderer::RenderOptionsChanged()
{
  state_=DIRTY_VA;
  sel_state_=DIRTY_VA;
}

}}}
