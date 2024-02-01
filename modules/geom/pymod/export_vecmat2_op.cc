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
#include <boost/python.hpp>
using namespace boost::python;

#include <ost/geom/geom.hh>
using namespace geom;
  
Real (*Vec2Length)(const Vec2& )                                      = &Length;
Real (*Vec2Length2)(const Vec2& )                                     = &Length2;
bool   (*Vec2Equal)(const Vec2& v1, const Vec2& v2, Real epsilon)     = &Equal;
bool   (*Mat2Equal)(const Mat2& m1, const Mat2& m2, Real epsilon)     = &Equal;
Real (*Mat2Dot)(const Vec2& v1, const Vec2& v2)                       = &Dot;
Vec2   (*Vec2CompMultiply)(const Vec2& v1, const Vec2& v2)              = &CompMultiply;
Vec2   (*Vec2CompDivide)(const Vec2& v1, const Vec2& v2)                = &CompDivide;
Real (*Mat2Det)(const Mat2& m)                                        = &Det;
Mat2   (*Mat2Invert)(const Mat2& m)                                     = &Invert;
Mat2   (*Mat2Transpose)(const Mat2& m)                                  = &Transpose;
Real (*Vec2Angle)(const Vec2& v1, const Vec2& v2)                     = &Angle;
Real (*Vec2SignedAngle)(const Vec2& v1, const Vec2& v2)                     = &SignedAngle;
Vec2   (*Vec2Normalize)(const Vec2& v1)                                 = &Normalize;
Vec2   (*Vec2Rotate)(const Vec2& v1,Real ang)                         = &Rotate;
Vec2 (*Vec2Min)(const Vec2&, const Vec2&) = &Min;
Vec2 (*Vec2Max)(const Vec2&, const Vec2&) = &Max;

void export_VecMat2_op()
{
  using namespace geom;
  
  def("Length",Vec2Length);
  def("Length2",Vec2Length2);
  def("Equal",Vec2Equal, (arg("v1"), arg("v2"), arg("epsilon")=EPSILON));
  def("CompMultiply",Vec2CompMultiply);
  def("CompDivide",Vec2CompDivide);
  def("Equal",Mat2Equal, (arg("m1"), arg("m2"), arg("epsilon")=EPSILON));
  def("Dot",Mat2Dot);
  def("Det",Mat2Det);
  def("Invert",Mat2Invert);
  def("Transpose",Mat2Transpose);
  def("Angle",Vec2Angle);
  def("SignedAngle",Vec2SignedAngle);
  def("Normalize",Vec2Normalize);
  def("Rotate",Vec2Rotate);
  def("Min",Vec2Min);
  def("Max",Vec2Max);
}
