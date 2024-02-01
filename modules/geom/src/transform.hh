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
#ifndef OST_GEOM_TRANSFORM_HH
#define OST_GEOM_TRANSFORM_HH

/*
  Author: Ansgar Philippsen
*/

#include <ost/config.hh>

#include "module_config.hh"
#include "vec3.hh"
#include "vec4.hh"
#include "mat3.hh"
#include "mat.hh"
#include "aligned_cuboid.hh"

namespace geom {

/// \brief basic and essential transformation class, including translation,
///    rotation and center of rotation
class DLLEXPORT_OST_GEOM Transform {
public:
  Transform();
  
  /// \brief reset to identity
  void Reset() {*this=Transform();}
  
  /// \brief retrieve transformation matrix
  Mat4 GetMatrix() const {return tm_;}
  /// \brief retrieve transposed transformation matrix
  Mat4 GetTransposedMatrix() const {return ttm_;}
  /// \brief retrieve inverted transformation matrix
  Mat4 GetInvertedMatrix() const {return itm_;}

  /// \brief set center of rotation
  void SetCenter(const Vec3& c);
  /// \brief retrieve center of rotation
  Vec3 GetCenter() const;

  // directly set tm, messing up rot/trans/cen !
  void SetMatrix(const Mat4& m);

  /// \name rotation
  //@{
  /// \brief rotate around x-axis
  void ApplyXAxisRotation(float delta);
  void ApplyYAxisRotation(float delta);
  void ApplyZAxisRotation(float delta);
  void ApplyAxisRotation(float delta, const Vec3& axis);  
  
  Mat3 GetXAxisRotation(float delta);
  Mat3 GetYAxisRotation(float delta);
  Mat3 GetZAxisRotation(float delta);
  
  Mat3 GetRot() const;
  void SetRot(const Mat3& r);
  //@}
  /// \brief translation
  //@{
  void ApplyXAxisTranslation(float delta);
  void ApplyYAxisTranslation(float delta);
  void ApplyZAxisTranslation(float delta);
  void SetTrans(const Vec3& t);
  Vec3 GetTrans() const;  
  //@}

  // apply to a vec3 and return result
  Vec3 Apply(const Vec3& v) const;
  // apply to a vec4 and return result
  Vec4 Apply(const Vec4& v) const;
  // apply inverse to a vec3 and return result
  Vec3 ApplyInverse(const Vec3& v) const;
  // apply inverse to a vec4 and return result
  Vec4 ApplyInverse(const Vec4& v) const;
  // apply to an aligned cuboid and return result
  AlignedCuboid Apply(const AlignedCuboid& c) const;
  Transform Apply(const Transform& tf) const;

private:
  Mat3 rot_;
  Vec3 trans_;
  Vec3 cen_;
  Mat4 tm_;
  Mat4 ttm_;
  Mat4 itm_;

  void update_tm();
  void update_components();
};

} // ns

#endif
