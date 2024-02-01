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
#ifndef OST_COLOR_OPS_BY_CHAIN_COLOR_OP_HH
#define OST_COLOR_OPS_BY_CHAIN_COLOR_OP_HH
#include <map>

#include <ost/info/info.hh>
#include <ost/info/info_fw.hh>

#include <ost/gfx/gradient.hh>
#include <ost/gfx/color_ops/color_op.hh>

/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX ByChainColorOp: public ColorOp {
public:
  ByChainColorOp();
  ByChainColorOp(const String& selection, int mask=DETAIL_COLOR|MAIN_COLOR);
  ByChainColorOp(const mol::QueryViewWrapper& query_view, int mask=DETAIL_COLOR|MAIN_COLOR);

  // Color Op interface
  virtual bool CanApplyTo(const GfxObjP& obj) const;
  virtual void ApplyTo(GfxObjP& obj) const;

  // this interface
  Color GetColor(const String& ident) const;
  unsigned int GetChainCount() const;
  void SetChainCount(unsigned int chain_count);

  //virtual void ToInfo(info::InfoGroup& group) const;
  static gfx::ByChainColorOp FromInfo(info::InfoGroup& group);

private:
  void Init();
  gfx::Color GenerateColor(String& ident) const;

  unsigned int chain_count_;
  float cm_; // 1 over chain_count
  mutable std::map<String,Color> colors_;

  gfx::Gradient color_grad_;
};

}}

#endif

