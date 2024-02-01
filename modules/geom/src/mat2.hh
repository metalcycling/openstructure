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
#ifndef GEOM_MAT2_HH
#define GEOM_MAT2_HH

#include <cstddef> // for size_t
#include <ostream>
#include <cassert>
#include <stdexcept>

#include <boost/operators.hpp>

#include <ost/geom/module_config.hh>

namespace geom {

class DLLEXPORT_OST_GEOM Mat2:
    private boost::equality_comparable<Mat2>,
    private boost::additive<Mat2>,
    private boost::multiplicative2<Mat2, Real>
{
public:
  static Mat2 Identity();

  //! Default initialization, identity matrix
  Mat2();

  //! In with 4 values in row-major order
  /*!
    row-major order means that the matrix

    | a b |
    | c d |

    is initialized with (a,b, c,d)
  */
  Mat2(Real i00, Real i01, Real i10, Real i11);

  //! Copy ctor
  Mat2(const Mat2& m);

  //! initialization from array
  explicit Mat2(const Real[4]);

  //! assignement op
  Mat2& operator=(const Mat2& m);

  //! comparable
  bool operator==(const Mat2& rhs) const;


  //! const element access
  const Real& At(std::size_t r, std::size_t c) const
  {
    if (r>1 || c>1) {
      throw std::out_of_range("indices must be smaller than 2");
    }
    return data_[r][c];
  }

  //! element access
  Real& At(std::size_t r, std::size_t c)
  {
    if (r>1 || c>1) {
      throw std::out_of_range("indices must be smaller than 2");
    }
    return data_[r][c];
  }

  //! element access
  Real& operator()(std::size_t r, std::size_t c)
  {
    assert(r<2 && c<2);
    return data_[r][c];
  }
  //! const element access
  const Real& operator()(std::size_t r, std::size_t c) const
  {
    assert(r<2 && c<2);
    return data_[r][c];
  }

  //! addable op
  Mat2& operator+=(const Mat2& rhs);
  //! subtractable op
  Mat2& operator-=(const Mat2& rhs);

  Mat2& operator*=(const Real d);
  Mat2& operator/=(const Real d);

  Mat2& operator*=(const Mat2& m);

  Real* Data() {return &data_[0][0];}
  const Real* Data() const {return &data_[0][0];}

private:
  Real data_[2][2];

  void set(Real i00, Real i01, Real i10, Real i11);
};

DLLEXPORT_OST_GEOM std::ostream& operator<<(std::ostream& os, const Mat2& m);

} // ns geom


#endif
