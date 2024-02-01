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
Author: Juergen Haas
*/

#ifndef OST_GFX_GL_HELPER_HH
#define OST_GFX_GL_HELPER_HH

/*
  must be included before any other gfx headers for gl and glext to
  be correctly 
*/
#ifdef OST_GFX_GLEXT_INCLUDE_HH
#error gl_helper.hh must be included before any occurence of glext_include.hh
#endif
#ifdef OST_GFX_GL_INCLUDE_HH
#error gl_helper.hh must be included before any occurence of gl_include.hh
#endif
#include <ost/geom/vec3.hh>

#include "glext_include.hh"

#include <ost/log.hh>

inline void check_gl_error(const std::string& m="")
{
  #ifndef NDEBUG
  GLenum error_code;
  if((error_code=glGetError())!=GL_NO_ERROR) {
    if(!m.empty()) {
      LOG_VERBOSE("GL error in [" << m << "]: " << gluErrorString(error_code));
    }
  }
  #endif
}

inline void glVertex3v(double* v){
  glVertex3dv(v);
}

inline void glVertex3v(const double* v){
  glVertex3dv(v);
}

inline void glVertex3v(float* v){
  glVertex3fv(v);
}


inline void glVertex3v(const float* v){
  glVertex3fv(v);
}

inline void glVertex3(const geom::Vec3& v)
{
  glVertex3v(&v[0]);
}

inline void glMultMatrix(float* v) {
  glMultMatrixf(v);
}

inline void glMultMatrix(double* v) {
  glMultMatrixd(v);
}

inline void glNormal3v(double* v){
  glNormal3dv(v);
}

inline void glNormal3v(float* v){
  glNormal3fv(v);
}

inline void glTexCoord2v(float* v){
  glTexCoord2fv(v);
}

inline void glTexCoord2v(double* v){
  glTexCoord2dv(v);
}


inline void glGetv(GLenum pname, double* v){
  glGetDoublev(pname, v);
}

inline void glGetv(GLenum pname, float* v){
  glGetFloatv(pname, v);
}


inline void glLoadMatrix(float* arr) {
  glLoadMatrixf(arr);
}

inline void glLoadMatrix(double* arr) {
  glLoadMatrixd(arr);
}

#if OST_SHADER_SUPPORT_ENABLED

inline void glLoadTransposeMatrix(float* arr) {
  glLoadTransposeMatrixf(arr);
}

inline void glLoadTransposeMatrix(double* arr) {
  glLoadTransposeMatrixd(arr);
}

#endif

#endif
