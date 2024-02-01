//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  transformation interface

  Author: Ansgar Philippsen
*/

#ifndef _IMG_TRANSFORMATION_H
#define _IMG_TRANSFORMATION_H

#include <cmath>

#include <ost/img/vecmat.hh>
#include <ost/img/extent.hh>
#include <ost/img/point.hh>

#include <ost/img/alg/module_config.hh>
namespace ost { namespace img { namespace alg {

//! transformation base class
/*!
  A transformation is defined by its matrix
  and the corresponding inverse. The transformation routine
  of the image class will use both to construct the result.

  The extent of the newly created images can be defined in 
  several ways: (1) size of the original image (2) maximal
  size to fit the transform

  Usage of this class allows an arbitrary transformation, in form
  of a 4x4 matrix, to be specified. If only the forward matrix is
  given, the inverse will be calculated automatically.
*/

class DLLEXPORT_IMG_ALG Transformation {
public:
  enum SizeRule {ORIGINAL, FIT};

  Transformation();
  Transformation(const Mat4& tmat);
  Transformation(const Mat4& tmat,const Mat4& imat);

  Mat4& Matrix() {return tmat_;}
  const Mat4& Matrix() const {return tmat_;}

  Mat4& InverseMatrix() {return imat_;}
  const Mat4& InverseMatrix() const {return imat_;}

  Point Apply(const Point& p);
  Extent Apply(const Extent& e);


protected:
  Mat4 tmat_;
  Mat4 imat_;
};

}}} // ns

#endif

