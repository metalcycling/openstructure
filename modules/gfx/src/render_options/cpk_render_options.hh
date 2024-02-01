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

#ifndef OST_GFX_CPK_RENDER_OPTIONS_HH
#define OST_GFX_CPK_RENDER_OPTIONS_HH

#include <boost/shared_ptr.hpp>
#include <ost/base.hh>

#include <ost/gfx/module_config.hh>

#include "render_options.hh"

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX CPKRenderOptions: public RenderOptions {
public:
  CPKRenderOptions();

  // RenderOptions interface
  virtual RenderMode::Type GetRenderMode();
  virtual bool CanApplyRenderOptions(RenderOptionsPtr render_options);
  virtual void ApplyRenderOptions(RenderOptionsPtr render_options);

  // own interface
  /// number of arc subdivisions per pi/2
  void SetSphereDetail(uint detail);
  uint GetSphereDetail();
  /// 0 = triangulated spheres, 1 = fast 3D sprites (default)
  void SetSphereMode(uint mode);
  uint GetSphereMode();
  /// Radius multiplier, default=1.0
  void SetRadiusMult(float m);
  float GetRadiusMult() const;


private:
  uint sphere_detail_;
  uint cpk_mode_;
  float rad_mult_;
};

typedef boost::shared_ptr<CPKRenderOptions> CPKRenderOptionsPtr;

}}

#endif /* OST_GFX_CPK_RENDER_OPTIONS_HH */
