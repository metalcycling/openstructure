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

#ifndef OST_GFX_RENDER_OPTIONS_HH
#define OST_GFX_RENDER_OPTIONS_HH
#include <vector>

#include <boost/shared_ptr.hpp>

#include <ost/gfx/module_config.hh>
#include <ost/gfx/render_mode.hh>
#include <ost/gfx/entity_fw.hh>



namespace ost { namespace gfx {

class RenderOptions;
typedef boost::shared_ptr<RenderOptions> RenderOptionsPtr;

class DLLEXPORT_OST_GFX RenderOptions {
public:
  RenderOptions();
  virtual ~RenderOptions(){}

  virtual RenderMode::Type GetRenderMode()=0;
  virtual bool CanApplyRenderOptions(RenderOptionsPtr render_options)=0;
  virtual void ApplyRenderOptions(RenderOptionsPtr render_options)=0;

  bool AddObserver( EntityP entity);
  bool RemoveObserver( EntityP entity);
  void NotifyStateChange();

private:
  typedef std::vector< EntityW > EntityWObservers;
  EntityWObservers observers_;

};

}} //ns

#endif /* OST_GFX_RENDER_OPTIONS_HH */
