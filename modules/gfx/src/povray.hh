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
#ifndef OST_GFX_POVRAY_HH
#define OST_GFX_POVRAY_HH

/*
  Author: Ansgar Philippsen
*/

#include <string>
#include <fstream>

#include <ost/geom/transform.hh>

#include "color.hh"

namespace ost { namespace gfx {

class PovState {
public:
  PovState(const std::string& pov_file, const std::string& inc_file, const std::string& wdir);

  void write_preamble();
  void write_postamble();
  void write_background(const Color& c);
  void write_camera(float fov, float zdist);
  void write_default_light();
  void write_fog(float zdist, float znear, float zfar, const Color& c, bool f);

  void start_obj(const std::string& name, float tp, float lw, float ps);
  void write_obj_tex(const std::string& name);
  void end_obj();

  void write_sphere(const geom::Vec3& p, float r, const Color& c, const std::string& name);
  void write_cyl(const geom::Vec3& p1, const geom::Vec3& p2, float r, const Color& c, const std::string& name, bool open);

  void write_merge_or_union(const std::string& name);

  std::string write_coord(const geom::Vec3& v);
  std::string write_coord(float* v);

  std::string write_norm(const geom::Vec3& n);
  std::string write_norm(float* n);

  std::ostream& pov() {return pov_;}
  std::ostream& inc() {return inc_;}

public:
  bool use_tf;
  geom::Transform tf;

private:
  std::string pov_file_;
  std::string inc_file_;
  std::string wdir_;
  std::ofstream pov_;
  std::ofstream inc_;
  std::vector<std::string> obj_list_;
};


}} // ns

#endif
