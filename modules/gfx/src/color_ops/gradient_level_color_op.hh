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
#ifndef OST_COLOR_OPS_GRADIENT_LEVEL_COLOR_OP_HH
#define OST_COLOR_OPS_GRADIENT_LEVEL_COLOR_OP_HH

#include <ost/mol/property_id.hh>



#include <ost/info/info.hh>
#include <ost/info/info_fw.hh>

#include <ost/gfx/gradient.hh>
#include <ost/gfx/color_ops/gradient_color_op.hh>

/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX GradientLevelColorOp: public GradientColorOp {
public:
  GradientLevelColorOp();
  GradientLevelColorOp(const String& selection, const String& property, const gfx::Gradient& gradient,
                       float minv, float maxv, 
                       mol::Prop::Level level=mol::Prop::UNSPECIFIED);

  GradientLevelColorOp(const String& selection, int mask, const String& property, const gfx::Gradient& gradient,
                       float minv, float maxv,
                       mol::Prop::Level level=mol::Prop::UNSPECIFIED);

  GradientLevelColorOp(const String& selection, const String& property, const gfx::Gradient& gradient,
                       mol::Prop::Level level=mol::Prop::UNSPECIFIED);

  GradientLevelColorOp(const String& selection, int mask, const String& property, const gfx::Gradient& gradient,
                       mol::Prop::Level level=mol::Prop::UNSPECIFIED);

  GradientLevelColorOp(const mol::QueryViewWrapper& query_view, const String& property, const gfx::Gradient& gradient,
                       float minv, float maxv,
                       mol::Prop::Level level=mol::Prop::UNSPECIFIED);

  GradientLevelColorOp(const mol::QueryViewWrapper& query_view, int mask, const String& property, const gfx::Gradient& gradient,
                       float minv, float maxv,
                       mol::Prop::Level level=mol::Prop::UNSPECIFIED);

  GradientLevelColorOp(const mol::QueryViewWrapper& query_view, const String& property, const gfx::Gradient& gradient,
                       mol::Prop::Level level=mol::Prop::UNSPECIFIED);

  GradientLevelColorOp(const mol::QueryViewWrapper& query_view, int mask, const String& property, const gfx::Gradient& gradient,
                       mol::Prop::Level level=mol::Prop::UNSPECIFIED);

  virtual bool CanApplyTo(const GfxObjP& obj) const;
  virtual void ApplyTo(GfxObjP& obj) const;

  virtual void SetLevel(mol::Prop::Level level);
  virtual mol::Prop::Level GetLevel() const;

  virtual void ToInfo(info::InfoGroup& group) const;
  static gfx::GradientLevelColorOp FromInfo(info::InfoGroup& group);

private:
  void Init();
  mol::Prop::Level level_;
};

}}

#endif

