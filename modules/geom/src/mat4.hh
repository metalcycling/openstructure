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
#ifndef GEOM_MAT4_HH
#define GEOM_MAT4_HH

#include <cassert>
#include <cstddef> // for size_t
#include <ostream>
#include <vector>
#include <stdexcept>

#include <boost/operators.hpp>

#include <ost/geom/module_config.hh>

namespace geom {

// fw decl
class Mat2;
class Mat3;
class Vec3;

class DLLEXPORT_OST_GEOM Mat4:
    private boost::equality_comparable<Mat4>,
    private boost::additive<Mat4>,
    private boost::multiplicative<Mat4, Real>
{
public:
  static Mat4 Identity();

  //! Default initialization, identity matrix
  Mat4();

  //! In with 16 values in row-major order
  /*!
    row-major order means that the matrix

    | a b c d |
    | e f g h |
    | i j k l |
    | m n o p |

    is initialized with (a,b,c,d, e,f,g,h, i,j,k,l, m,n,o,p)
  */

  Mat4(Real i00, Real i01, Real i02, Real i03, Real i10, Real i11, Real i12, Real i13,
       Real i20, Real i21, Real i22, Real i23, Real i30, Real i31, Real i32, Real i33);
  //! Copy ctor
  Mat4(const Mat4& m);

  //! implicit initialization from mat2
  explicit Mat4(const Mat2& m);

  //! implicit initialization from mat3
  explicit Mat4(const Mat3& m);

  //! initialization from array
  explicit Mat4(const float[16]);
  
  explicit Mat4(const double[16]);
  //! assignement op
  Mat4& operator=(const Mat4& m);

  //! comparable
  bool operator==(const Mat4& rhs) const;

  const Real& At(std::size_t r, std::size_t c) const
  {
    if (r>3 || c>3) {
      throw std::out_of_range("indices must be smaller than 4");
    }
    return data_[r][c];
  }

  Real& At(std::size_t r, std::size_t c)
  {
    if (r>3 || c>3) {
      throw std::out_of_range("indices must be smaller than 4");
    }
    return data_[r][c];
  }

  Real& operator()(std::size_t r, std::size_t c)
  {
    assert(r<4 && c < 4);
    return data_[r][c];
  }

  const Real& operator()(std::size_t r, std::size_t c) const
  {
    assert(r<4 && c < 4);
    return data_[r][c];
  }

  //! addable op
  Mat4& operator+=(const Mat4& rhs);
  //! subtractable op
  Mat4& operator-=(const Mat4& rhs);

  Mat4& operator*=(const Real d);
  Mat4& operator/=(const Real d);

  Mat4& operator*=(const Mat4& m);

  Mat3 ExtractRotation() const;
  void PasteRotation(const Mat3& m);
  Vec3 ExtractTranslation() const;
  void PasteTranslation(const Vec3& v);

  Real* Data() {return &data_[0][0];}
  const Real* Data() const {return &data_[0][0];}

private:
  Real data_[4][4];

  void set(float i00, float i01, float i02, float i03, float i10, float i11, float i12,
float i13, float i20, float i21, float i22, float i23, float i30, float i31, float i32,
float i33);
  void set(double i00, double i01, double i02, double i03, double i10, double i11, double
i12, double i13, double i20, double i21, double i22, double i23, double i30, double i31,
double i32, double i33);

};

typedef std::vector<Mat4> Mat4List;

DLLEXPORT_OST_GEOM std::ostream& operator<<(std::ostream& os, const Mat4& m);

} // ns geom


#endif
