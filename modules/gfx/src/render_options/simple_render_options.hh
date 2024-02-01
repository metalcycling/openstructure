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

#ifndef OST_GFX_SIMPLE_RENDER_OPTIONS_HH
#define OST_GFX_SIMPLE_RENDER_OPTIONS_HH

#include <boost/shared_ptr.hpp>

#include <ost/gfx/module_config.hh>

#include "line_render_options.hh"

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX SimpleRenderOptions: public LineRenderOptions {
public:
  SimpleRenderOptions();

  // remaining RenderOptions interface not define in LineRenderOptions
  virtual RenderMode::Type GetRenderMode();
  virtual bool CanApplyRenderOptions(RenderOptionsPtr render_options);

  // own interface
  bool GetBlurFlag() const;
  bool GetBondOrderFlag() const;
  Real GetBondOrderDistance() const;
  void SetBlurFlag(bool flag);
  const std::pair<Real, Real>& GetBlurFactors() const;
  void SetBlurFactors(Real bf1, Real bf2);
  void SetBondOrderFlag(bool flag);
  void SetBondOrderDistance(Real bod);

private:
  bool                     blur_flag_;
  std::pair<Real, Real>    blur_factors_;
  bool                     bond_order_flag_;
  Real                     bond_order_distance_;
};

typedef boost::shared_ptr<SimpleRenderOptions> SimpleRenderOptionsPtr;

}

}

#endif /* OST_GFX_SIMPLE_RENDER_OPTIONS_HH */
