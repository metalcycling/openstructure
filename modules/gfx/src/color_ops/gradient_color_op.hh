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
#ifndef OST_COLOR_OPS_GRADIENT_COLOR_OP_HH
#define OST_COLOR_OPS_GRADIENT_COLOR_OP_HH

#include <ost/info/info.hh>
#include <ost/info/info_fw.hh>

#include <ost/gfx/gradient.hh>

#include <ost/gfx/color_ops/color_op.hh>

/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX GradientColorOp: public ColorOp {
public:
  GradientColorOp();
  GradientColorOp(const String& selection, const String& property, const gfx::Gradient& gradient,
                  float minv, float maxv);
  GradientColorOp(const String& selection, int mask, const String& property, const gfx::Gradient& gradient,
                    float minv, float maxv);
  GradientColorOp(const String& selection, const String& property, const gfx::Gradient& gradient);
  GradientColorOp(const String& selection, int mask, const String& property, const gfx::Gradient& gradient);
  GradientColorOp(const mol::QueryViewWrapper& query_view, const String& property, const gfx::Gradient& gradient,
                  float minv, float maxv);
  GradientColorOp(const mol::QueryViewWrapper& query_view, int mask, const String& property, const gfx::Gradient& gradient,
                    float minv, float maxv);
  GradientColorOp(const mol::QueryViewWrapper& query_view, const String& property, const gfx::Gradient& gradient);
  GradientColorOp(const mol::QueryViewWrapper& query_view, int mask, const String& property, const gfx::Gradient& gradient);

  virtual void SetProperty(const String& property);
  virtual String GetProperty() const;

  virtual void SetGradient(const gfx::Gradient& gradient);
  virtual gfx::Gradient GetGradient() const;

  virtual bool GetCalculateMinMax() const;
  virtual void SetMaxV(float maxv) const;
  virtual float GetMaxV() const;

  virtual void SetMinV(float minv) const;
  virtual float GetMinV() const;

  virtual void ToInfo(info::InfoGroup& group) const;
  static gfx::GradientColorOp FromInfo(info::InfoGroup& group);

  // if true then color with clamped colors outside given range as well
  void SetClamp(bool f) {clamp_=f;}
  bool GetClamp() const {return clamp_;}

private:
  void Init();
  String property_;
  gfx::Gradient gradient_;
  bool clamp_;
  bool calculate_;
  mutable float minv_;
  mutable float maxv_;
};

}}

#endif

