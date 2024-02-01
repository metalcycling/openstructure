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
  Author: Ansgar Philippsen
*/

#include "gfx_prim.hh"

namespace ost { namespace gfx {

void CylinderPrim::calc_rotmat()
{
  geom::Vec3 ref=geom::Normalize(end-start);
  
  float xi = std::acos(ref[2]);
  float om = ref[1]!=0.0 ? std::atan2(ref[1],ref[0]) : (ref[0]>=0.0 ? 0.0 : -M_PI);
  
  rotmat_t = geom::Mat3(cos(xi),0.0,-sin(xi),0.0,1.0,0.0,sin(xi),0.0,cos(xi))*
    geom::Mat3(cos(om),sin(om),0.0,-sin(om),cos(om),0.0,0.0,0.0,1.0);
  rotmat = geom::Transpose(rotmat_t);
}

}} // ns
