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
  Author: Stefan Scheuber
*/
#include <ost/gfx/input.hh>
#include <ost/gfx/entity.hh>

#include "render_options.hh"

namespace ost { namespace gfx {

RenderOptions::RenderOptions(){}

bool RenderOptions::AddObserver(EntityP entity)
{
  for (EntityWObservers::const_iterator
       i = observers_.begin(), e = observers_.end(); i != e; ++i) {
    if (i->lock() == entity)
      return false;
  }
  observers_.push_back(entity);
  return true;
}

bool RenderOptions::RemoveObserver(EntityP entity)
{
  for (EntityWObservers::iterator i = observers_.begin(), e = observers_.end(); 
       i != e; ++i) {
    if (i->lock() == entity) {
      observers_.erase(i);
      return true;
    }
    
  }
  return false;
}

void RenderOptions::NotifyStateChange()
{
  for (EntityWObservers::iterator
       i = observers_.begin(), e = observers_.end(); i != e; ++i) {
    i->lock()->OptionsChanged(this->GetRenderMode());
  }
}

}} // ns
