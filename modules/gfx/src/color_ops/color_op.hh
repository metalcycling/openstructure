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
#ifndef OST_COLOR_OPS_COLOR_OP_HH
#define OST_COLOR_OPS_COLOR_OP_HH

#include <string>

#include <ost/mol/query_view_wrapper.hh>

#include <ost/gfx/gfx_object_fw.hh>
#include <ost/gfx/module_config.hh>

#include <ost/info/info.hh>
#include <ost/info/info_fw.hh>

/*
  Author: Stefan Scheuber
*/

namespace ost { namespace gfx {

typedef enum {
  DETAIL_COLOR = 1,
  MAIN_COLOR = 2
} ColorBit;

typedef unsigned char ColorMask;
  
class DLLEXPORT_OST_GFX ColorOp {
public:
  ColorOp();
  virtual ~ColorOp() {}
  ColorOp(const String& selection, int mask=DETAIL_COLOR|MAIN_COLOR);
  ColorOp(const mol::QueryViewWrapper& query_view, int mask=DETAIL_COLOR|MAIN_COLOR);
  virtual bool CanApplyTo(const GfxObjP& obj) const;
  virtual void ApplyTo(GfxObjP& obj) const;

  virtual const String& GetName() const;

  ColorMask GetMask() const;
  
  void SetMask(ColorMask mask);
  
  virtual void SetSelection(const String& selection);
  virtual String GetSelection() const;

  virtual void SetSelectionFlags(mol::QueryFlags flags);
  virtual mol::QueryFlags GetSelectionFlags() const;

  virtual bool IsSelectionOnly() const;
  virtual void SetView(const mol::EntityView& view);
  virtual mol::EntityView GetView() const;

  virtual void ToInfo(info::InfoGroup& group) const;
  static gfx::ColorOp FromInfo(info::InfoGroup& group);

protected:
  virtual void SetName(const String& name);

private:
  String name_;
  mol::QueryViewWrapper query_view_;
  ColorMask mask_;
};

}}

#endif

