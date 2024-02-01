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
#include <math.h>

#include "by_chain_color_op.hh"

#include <ost/gfx/entity.hh>
#include <ost/dyn_cast.hh>


namespace ost { namespace gfx {

ByChainColorOp::ByChainColorOp() : ColorOp(),chain_count_(0),color_grad_("RAINBOW")
{
  this->Init();
}

ByChainColorOp::ByChainColorOp(const String& selection, int mask) : ColorOp(selection,mask),chain_count_(0),color_grad_("RAINBOW")
{
  this->Init();
}

ByChainColorOp::ByChainColorOp(const mol::QueryViewWrapper& query_view, int mask) : ColorOp(query_view,mask),chain_count_(0),color_grad_("RAINBOW")
{
  this->Init();
}

void ByChainColorOp::Init()
{
  this->SetName("By chain");
}

bool ByChainColorOp::CanApplyTo(const GfxObjP& obj) const{
  if(dynamic_cast<Entity*>(obj.get()))
    return true;
  return false;
}

void ByChainColorOp::ApplyTo(GfxObjP& objP) const{
  GfxObj* obj = objP.get();
  if(Entity* ent=dynamic_cast<Entity*>(obj)){
    ent->Apply(*this, false);
  }
}

/*
  The first chain will be assigned a gradient stop of 0.
  The next chain will be assigned a gradient stop of n/(N-1), where
  n is the number of chains assigned so far, and N is the total
  number of chains; thus for 3 total chains, the second gets a stop
  of 0.5 (1/(3-1)), and the third gets a stop of 1.0 (2/(3-1))
*/
gfx::Color ByChainColorOp::GetColor(const String& ident) const
{
  std::map<String,Color>::const_iterator cit = colors_.find(ident);
  if(cit == colors_.end()) {
    float f = static_cast<float>(colors_.size())*cm_;
    Color c = color_grad_.GetColorAt(f);
    colors_[ident] = c;
    return c;
  } else {
    return cit->second;
  }
}

unsigned int ByChainColorOp::GetChainCount() const
{
  return chain_count_;
}

void ByChainColorOp::SetChainCount(unsigned int chain_count)
{
  chain_count_ = chain_count;
  cm_=chain_count>1 ? 1.0/static_cast<float>(chain_count_-1) : 0.0;
  // invalidate all assigned colors so far
  colors_.clear();
}

gfx::ByChainColorOp ByChainColorOp::FromInfo(info::InfoGroup& group){
  gfx::ColorOp op = ColorOp::FromInfo(group);
  mol::QueryViewWrapper wrapper(op.GetSelection(),op.GetSelectionFlags());
  int mask = op.GetMask();
  return gfx::ByChainColorOp(wrapper, mask);
}

}}

