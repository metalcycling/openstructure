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
#ifndef OST_GFX_PRIM_HH
#define OST_GFX_PRIM_HH

/*
  Author: Ansgar Philippsen
*/

#include <vector>

#include <ost/geom/geom.hh>

#include "color.hh"

namespace ost { namespace gfx {

struct SpherePrim {
  SpherePrim():
    position(),radius(1.0),color() 
  {}
  SpherePrim(const geom::Vec3& pos, float rad, const Color& col):
    position(pos), radius(rad), color(col) 
  {}

  geom::Vec3 position;
  float radius;
  Color color;
};

typedef std::vector<SpherePrim> SpherePrimList;

struct CylinderPrim {
  CylinderPrim():
    start(), end(),
    radius1(1.0), radius2(1.0),
    color1(), color2(),
    length(1.0), rotmat(), rotmat_t()
  {
    calc_rotmat();
  }

  CylinderPrim(const geom::Vec3& st, const geom::Vec3& en, float rad, const Color& col):
    start(st), end(en),
    radius1(rad), radius2(rad),
    color1(col), color2(col),
    length(geom::Length(end-start)), rotmat(), rotmat_t() 
  {
    calc_rotmat();
  }

  CylinderPrim(const geom::Vec3& st, const geom::Vec3& en, float rad, const Color& col1, const Color& col2):
    start(st), end(en),
    radius1(rad), radius2(rad),
    color1(col1), color2(col2),
    length(geom::Length(end-start)), rotmat(), rotmat_t() 
  {
    calc_rotmat();
  }

  CylinderPrim(const geom::Vec3& st, const geom::Vec3& en, float radius_1, float radius_2, const Color& col):
    start(st), end(en),
    radius1(radius_1), radius2(radius_2),
    color1(col), color2(col),
    length(geom::Length(end-start)), rotmat(), rotmat_t() 
  {
    calc_rotmat();
  }

  CylinderPrim(const geom::Vec3& st, const geom::Vec3& en, float radius_1, float radius_2, const Color& col1, const Color& col2):
    start(st), end(en),
    radius1(radius_1), radius2(radius_2),
    color1(col1), color2(col2),
    length(geom::Length(end-start)), rotmat(), rotmat_t() 
  {
    calc_rotmat();
  }

  void calc_rotmat();

  geom::Vec3 start,end;
  float radius1,radius2;
  Color color1, color2;
  float length;
  geom::Mat3 rotmat;
  geom::Mat3 rotmat_t;
};

typedef std::vector<CylinderPrim> CylinderPrimList;

struct TextPrim {
  TextPrim(): str(""), position(),color(),points(1.0) {}
  TextPrim(const String& s, const geom::Vec3& p, const Color& c, float ps):
    str(s), position(p), color(c), points(ps) {}
  String str;
  geom::Vec3 position;
  Color color;
  float points;
};

typedef std::vector<TextPrim> TextPrimList;


}} // ns

#endif
