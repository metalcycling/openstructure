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
#ifndef OST_COLOR_OPS_MAP_HANDLE_COLOR_OP_HH
#define OST_COLOR_OPS_MAP_HANDLE_COLOR_OP_HH

#include <ost/gfx/gradient.hh>

#include <ost/info/info_fw.hh>

#include <ost/img/map.hh>

#include <ost/gfx/color_ops/gradient_color_op.hh>


/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX MapHandleColorOp: public GradientColorOp {
public:
  MapHandleColorOp();
  MapHandleColorOp(const String& selection, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const img::MapHandle& mh);
  MapHandleColorOp(const String& selection, int mask, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const img::MapHandle& mh);
  MapHandleColorOp(const mol::QueryViewWrapper& query_view, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const img::MapHandle& mh);
  MapHandleColorOp(const mol::QueryViewWrapper& query_view, int mask, const String& property, const gfx::Gradient& gradient, float minv, float maxv, const img::MapHandle& mh);

  virtual bool CanApplyTo(const GfxObjP& obj) const;
  virtual void ApplyTo(GfxObjP& obj) const;

  virtual void SetMapHandle(const img::MapHandle& mh);
  virtual const img::MapHandle& GetMapHandle() const;

  virtual void ToInfo(info::InfoGroup& group) const;
  static gfx::MapHandleColorOp FromInfo(info::InfoGroup& group);

private:
  void Init();
  img::MapHandle mh_;
};

}}

#endif

