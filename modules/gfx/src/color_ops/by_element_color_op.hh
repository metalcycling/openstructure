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
#ifndef OST_COLOR_OPS_BY_ELEMENT_COLOR_OP_HH
#define OST_COLOR_OPS_BY_ELEMENT_COLOR_OP_HH

#include <boost/shared_ptr.hpp>

#include <ost/info/info.hh>
#include <ost/info/info_fw.hh>

#include <ost/gfx/color_ops/color_op.hh>

/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX ByElementColorOp: public ColorOp {
public:
  ByElementColorOp();

  ByElementColorOp(const mol::QueryViewWrapper& query_view, int mask=DETAIL_COLOR|MAIN_COLOR);
  ByElementColorOp(const String& selection, int mask=DETAIL_COLOR|MAIN_COLOR);


  virtual bool CanApplyTo(const GfxObjP& obj) const;
  virtual void ApplyTo(GfxObjP& obj) const;

  static ost::gfx::ByElementColorOp FromInfo(info::InfoGroup& group);
private:
  void Init();
};

}}

#endif

