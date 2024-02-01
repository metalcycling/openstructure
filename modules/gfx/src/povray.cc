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

#include <sstream>
#include <boost/format.hpp>

#include "povray.hh"

using namespace ost;
using namespace ost::gfx;

PovState::PovState(const std::string& pfile, const std::string& ifile, const std::string& wdir):
  use_tf(false),
  tf(),
  pov_file_(pfile),
  inc_file_(ifile),
  wdir_(wdir),
  pov_((wdir_+pov_file_).c_str()),
  inc_((wdir_+inc_file_).c_str()),
  obj_list_()
{
}

void PovState::write_preamble()
{
  pov_ << "#version 3.6;\n\n";
}

void PovState::write_postamble()
{
  pov_ << boost::format("#include \"%s\"\n") % inc_file_;
  pov_ << "union {\n";
  for(std::vector<std::string>::const_iterator it=obj_list_.begin();it!=obj_list_.end();++it) {
    pov_ << " object {_" << (*it) << "}\n";
  }
  if(!use_tf) {
    geom::Vec3 cen=tf.GetCenter();
    geom::Vec3 tra=tf.GetTrans();
    geom::Mat3 rot=tf.GetRot();
    pov_ << " /*\n"
      "  When uncommented, this applies a clipped\n"
      "  object in POV camera coordinates, but then\n" 
      "  transforms it into the original object\n"
      "  coordinates. This slightly complicated\n"
      "  procedure is necessary for clipping to work\n"
      "  properly\n"
      " */\n";
    pov_ << " /*\n";
    pov_ << " clipped_by {\n";
    pov_ << "  // in camera coordinates\n";
    pov_ << "  plane {z 0}\n";
    pov_ << "  // object transformation\n";
    pov_ << boost::format("  //translate <%.4f,%.4f,%.4f>\n") % -tra[0] % -tra[1] % 0.0;
    pov_ << boost::format("  matrix <%.4f,%.4f,%.4f,\n"
                          "          %.4f,%.4f,%.4f,\n"
                          "          %.4f,%.4f,%.4f,\n"
                          "          0,0,0>\n") % rot(0,0) % rot(0,1) % rot(0,2) % rot(1,0) % rot(1,1) % rot(1,2) % rot(2,0) % rot(2,1) % rot(2,2);
    pov_ << boost::format(" translate <%.4f,%.4f,%.4f>\n") % cen[0] % cen[1] % cen[2];
    pov_ << " }\n */\n\n";
    
    
    pov_ << " /*\n";
    pov_ << "  This transformation, together with\n";
    pov_ << "  the camera translation above,\n";
    pov_ << "  mimics the original scene view:\n";
    pov_ << "   1. translate to center of rotation\n";
    pov_ << "   2. apply rotation matrix\n";
    pov_ << "   3. apply xy translation (z component is in camera)\n";
    pov_ << " */\n";
    pov_ << boost::format(" translate <%.4f,%.4f,%.4f>\n") % -cen[0] % -cen[1] % -cen[2];
    pov_ << boost::format(" matrix <%.4f,%.4f,%.4f,\n"
                          "         %.4f,%.4f,%.4f,\n"
                          "         %.4f,%.4f,%.4f,\n"
                          "         0,0,0>\n") % rot(0,0) % rot(1,0) % rot(2,0) % rot(0,1) % rot(1,1) % rot(2,1) % rot(0,2) % rot(1,2) % rot(2,2);
    pov_ << boost::format(" //translate <%.4f,%.4f,%.4f>\n") % tra[0] % tra[1] % 0.0;
  }
  pov_ << "}\n\n";
}

void PovState::write_background(const Color& c)
{
  pov_ << boost::format("background {color rgb <%.4f,%.4f,%.4f>}\n\n") % c[0] % c[1] % c[2];
}

void PovState::write_camera(float fov, float zdist)
{
  pov_ << boost::format(
   "/*\n"
   "  This camera setup is the same as in GL; the view is\n"
   "  along the negative z-axis, and the camera is translated\n"
   "  towards the positive z, looking at the origin. The\n"
   "  x-axis points to the right, and the y-axis up.\n"
   "  ori_dist is the distance along the positive z (out from\n"
   "  the screen towards the viewer) from the origin to the\n"
   "  camera\n"
   "*/\n\n"
   "#declare ori_dist = %.4f;\n\n"
   "camera {\n"
   " perspective\n"
   " location <0,0,0>\n"
   " direction <0,0,-1>\n"
   " up <0,1,0>\n"
   " right <1,0,0>\n"
   " angle %.4f\n"
   " translate <0,0,ori_dist>\n"
   "}\n\n"
   ) % zdist % fov;
}

void PovState::write_default_light()
{
  pov_ << "light_source {<0,0,ori_dist> color rgb 1}\n\n";
}

void PovState::write_fog(float zdist, float znear, float zfar, const Color& c,bool f)
{
  pov_ << 
    "/*\n"
    "  Rotated ground fog is used to start the depth\n"
    "  effect at the precise starting offset (fog_offset).\n"
    "  In contrast to OpenGL, the fog equation here has\n"
    "  an exponential decay, and hence the end is\n"
    "  approximated (distance). Decrease distance to get\n"
    "  a faster falloff, and hence a stronger fog.\n"
    "*/\n";
  pov_ << "#if (" << (f ? "1" : "0") << ") // 1=on, 0=off\n";
  pov_ << boost::format(
   "fog {\n"
   " fog_type 2\n"
   " distance %.4f\n"
   " color rgb <%.4f,%.4f,%.4f>\n"
   " fog_offset %.4f\n"
   " fog_alt 0.01\n"
   " up <0,0,1>\n"
   "}\n"
   "#end\n\n") % (zfar-znear) % c[0] % c[1] % c[2] % (zdist-znear);
}

void PovState::start_obj(const std::string& name, float tp, float lw, float ps)
{
  pov_ << boost::format(
   "/* Definitions for %s */\n"
   "#declare _%s_tex = texture {\n"
   " pigment {color rgbft <0,0,0,1,1>}\n"
   " normal {granite 0.0 scale 1.0}\n"
   " finish {\n"
   "  ambient 0.10\n"
   "  diffuse 0.80 brilliance 1.00\n"
   "  specular 0.10 roughness 0.02\n"
   "  reflection 0.0\n"
   " }\n"
   "}\n"
   "#declare _%s_mat = material {texture {_%s_tex}}\n"
   "#declare _%s_tp = %.4f;  // transparency in pigment color\n"
   "#declare _%s_fi = %.4f;  // filter in pigment color\n"
   "#declare _%s_rf = 1.0; // multiplier for sphere and cyl radius\n"
   "#declare _%s_merge = 0; // use merge for objects, union otherwise (if applicable)\n"
   "\n\n"
   ) % name % name % name % name % name % (1.0f-tp) % name % 0.0f % name % name;

  inc_ << boost::format(
   "#declare _%s = object {\n"
   ) % name;

  obj_list_.push_back(name);
}

void PovState::end_obj()
{
  inc_ << "}\n";
}

void PovState::write_obj_tex(const std::string& name)
{
  for(unsigned int r=0;r<8;++r) {
    float fr=static_cast<float>(r)/7.0;
    for(unsigned int g=0;g<8;++g) {
      float fg=static_cast<float>(g)/7.0;
      for(unsigned int b=0;b<8;++b) {
        float fb=static_cast<float>(b)/7.0;
        inc_ << boost::format("#declare _%s_tex_%1d%1d%1d = texture {_%s_tex pigment {color rgbft <%.4f,%.4f,%.4f,_%s_fi,_%s_tp>}}\n") % name % r % g % b % name % fr % fg % fb % name % name ;
      }
    }
  }
}
void PovState::write_sphere(const geom::Vec3& p, float r, const Color& c, const std::string& name)
{
  inc_ << boost::format(" sphere {%s, _%s_rf*%.4f") % write_coord(p) % name % r;
  inc_ << boost::format(" texture { _%s_tex pigment {color rgbft <%.4f,%.4f,%.4f,_%s_fi,_%s_tp>}}}\n") % name % c[0] % c[1] % c[2] % name % name;
}

void PovState::write_cyl(const geom::Vec3& p1, const geom::Vec3& p2, float r, const Color& c, const std::string& name, bool open)
{
  inc_ << boost::format(" cylinder {%s, %s, _%s_rf*%.4f") % write_coord(p1) % write_coord(p2) % name % r;
  if(open) inc_ << " open";
  inc_ << boost::format(" texture { _%s_tex pigment {color rgbft <%.4f,%.4f,%.4f,_%s_fi,_%s_tp>}}}\n") % name % c[0] % c[1] % c[2] % name % name;
}

void PovState::write_merge_or_union(const std::string& name)
{
  inc_ << boost::format("#if (_%s_merge)\n") % name;
  inc_ << " merge {\n"
          "#else\n"
          " union {\n"
          "#end\n";
}

std::string PovState::write_coord(const geom::Vec3& c)
{
  geom::Vec3 t = (use_tf) ? tf.Apply(c) : c;
  std::ostringstream s;
  s << boost::format("<%.4f,%.4f,%.4f>") % t[0] % t[1] % t[2];
  return s.str();
}

std::string PovState::write_coord(float* v)
{
  return write_coord(geom::Vec3(v));
}

std::string PovState::write_norm(const geom::Vec3& c)
{
  geom::Vec3 t = (use_tf) ? tf.GetRot()*c : c;
  std::ostringstream s;
  s << boost::format("<%.4f,%.4f,%.4f>") % t[0] % t[1] % t[2];
  return s.str();
}

std::string PovState::write_norm(float* v)
{
  return write_norm(geom::Vec3(v));
}
