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

#include "by_element_color_op.hh"

#include <ost/gfx/entity.hh>
#include <ost/dyn_cast.hh>


namespace ost { namespace gfx {

ByElementColorOp::ByElementColorOp() : ColorOp()
{
  this->Init();
}

ByElementColorOp::ByElementColorOp(const String& selection, int mask) : ColorOp(selection,mask)
{
  this->Init();
}

ByElementColorOp::ByElementColorOp(const mol::QueryViewWrapper& query_view, int mask) : ColorOp(query_view,mask)
{
  this->Init();
}

void ByElementColorOp::Init()
{
  this->SetName("By element");
}

bool ByElementColorOp::CanApplyTo(const GfxObjP& obj) const{
  if(dynamic_cast<Entity*>(obj.get()))
    return true;
  return false;
}

void ByElementColorOp::ApplyTo(GfxObjP& objP) const{
  GfxObj* obj = objP.get();
  if(Entity* ent=dynamic_cast<Entity*>(obj)){
    ent->Apply(*this, false);
  }
}

gfx::ByElementColorOp ByElementColorOp::FromInfo(info::InfoGroup& group){
  gfx::ColorOp op = ColorOp::FromInfo(group);
  mol::QueryViewWrapper wrapper(op.GetSelection(),op.GetSelectionFlags());
  int mask = op.GetMask();
  return gfx::ByElementColorOp(wrapper, mask);
}

}}

