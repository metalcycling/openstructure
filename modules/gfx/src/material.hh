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
#ifndef OST_GFX_MATERIAL_HH
#define OST_GFX_MATERIAL_HH

#include "color.hh"

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX Material
{
public:
  Material();
  Material(const Color& amb, const Color& diff, const Color& spec, float shin, const Color& emm);
  Material(float amb, float diff, float spec, float shin, float emm);
  void SetAmb(const Color& c);
  void SetAmb(float c);
  void SetDiff(const Color& c);
  void SetDiff(float c);
  void SetSpec(const Color& c);
  void SetSpec(float c);
  void SetShin(float s);
  void SetEmm(const Color& c);
  void SetEmm(float c);

  void RenderGL();

private:
  Color amb_;
  Color diff_;
  Color spec_;
  Color emm_;
  float shin_;
};

}} // ns


#endif
