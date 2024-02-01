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

Real (*Vec4Length)(const Vec4& )                                      = &Length;
Real (*Vec4Length2)(const Vec4& )                                     = &Length2;
bool   (*Vec4Equal)(const Vec4& v1, const Vec4& v2, Real epsilon)     = &Equal;
bool   (*Mat4Equal)(const Mat4& m1, const Mat4& m2, Real epsilon)     = &Equal;
Real (*Mat4Dot)(const Vec4& v1, const Vec4& v2)                       = &Dot;
Vec4   (*Vec4CompMultiply)(const Vec4& v1, const Vec4& v2)              = &CompMultiply;
Vec4   (*Vec4CompDivide)(const Vec4& v1, const Vec4& v2)                = &CompDivide;
Real (*Mat4Det)(const Mat4& m)                                        = &Det;
Mat4   (*Mat4Invert)(const Mat4& m)                                     = &Invert;
Mat4   (*Mat4Transpose)(const Mat4& m)                                  = &Transpose;
Real (*Mat4Comp)(const Mat4& m, unsigned int i, unsigned int j)       = &Comp;
Real (*Mat4Minor)(const Mat4& m, unsigned int i, unsigned int j)      = &Minor;
Real (*Vec4Angle)(const Vec4& v1, const Vec4& v2)                     = &Angle;
Vec4 (*Vec4Min)(const Vec4&, const Vec4&) = &Min;
Vec4 (*Vec4Max)(const Vec4&, const Vec4&) = &Max;


void export_VecMat4_op()
{
  using namespace geom;
  
  def("Equal",Vec4Equal, (arg("v1"), arg("v2"), arg("epsilon")=EPSILON));
  def("Length2",Vec4Length2);
  def("Equal",Vec4Equal);
  def("CompMultiply",Vec4CompMultiply);
  def("CompDivide",Vec4CompDivide);
  def("Equal",Mat4Equal, (arg("m1"), arg("m2"), arg("epsilon")=EPSILON));
  def("Dot",Mat4Dot);
  def("Det",Mat4Det);
  def("Transpose",Mat4Transpose);
  def("Invert",Mat4Invert);
  def("Comp",Mat4Comp);
  def("Minor",Mat4Minor);
  def("Angle",Vec4Angle);
  def("Min",Vec4Min);
  def("Max",Vec4Max);
}
