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

Real (*Vec3Length)(const Vec3& )                                      = &Length;
Real (*Vec3Length2)(const Vec3& )                                     = &Length2;
bool   (*Vec3Equal)(const Vec3& v1, const Vec3& v2, Real epsilon)     = &Equal;
bool   (*Mat3Equal)(const Mat3& m1, const Mat3& m2, Real epsilon)     = &Equal;
Real (*Mat3Dot)(const Vec3& v1, const Vec3& v2)                       = &Dot;
Vec3   (*Vec3CompMultiply)(const Vec3& v1, const Vec3& v2)              = &CompMultiply;
Real (*Vec3Distance)(const Vec3& v1, const Vec3& v2)                  = &Distance;
Vec3   (*Vec3CompDivide)(const Vec3& v1, const Vec3& v2)                = &CompDivide;
Real (*Mat3Det)(const Mat3& m)                                        = &Det;
Vec3   (*Vec3Cross)(const Vec3& v1, const Vec3& v2)                     = &Cross;
Mat3   (*Mat3Invert)(const Mat3& m)                                     = &Invert;
Mat3   (*Mat3Transpose)(const Mat3& m)                                  = &Transpose;
Real (*Mat3Comp)(const Mat3& m, unsigned int i, unsigned int j)       = &Comp;
Real (*Mat3Minor)(const Mat3& m, unsigned int i, unsigned int j)      = &Minor;
Vec3 (*Vec3Min)(const Vec3&, const Vec3&) = &Min;
Vec3 (*Vec3Max)(const Vec3&, const Vec3&) = &Max;
Real (*Vec3SignedAngle)(const Vec3& v1, const Vec3& v2, const Vec3& v3) = &SignedAngle;
Real (*Vec3Distance2WithPBC)(const Vec3&, const Vec3&, const Vec3&)   = &Distance2WithPBC;
Real (*Vec3DistanceWithPBC)(const Vec3&, const Vec3&, const Vec3&)    = &DistanceWithPBC;
Vec3 (*wrap_vec3_1)(const Vec3&, const Vec3&, const Vec3&)               = &WrapVec3;
Vec3 (*wrap_vec3_2)(const Vec3&, const Vec3&, const Vec3&, const Vec3&)  = &WrapVec3;
Vec3List (*wrap_vec3list1)(const Vec3List&, const Vec3&, const Vec3&)               = &WrapVec3List;
Vec3List (*wrap_vec3list2)(const Vec3List&, const Vec3&, const Vec3&, const Vec3&)  = &WrapVec3List;
void export_VecMat3_op()
{
  using namespace geom;
  
  def("Length",Vec3Length);
  def("Length2",Vec3Length2);
  def("Equal",Vec3Equal, (arg("v1"), arg("v2"), arg("epsilon")=EPSILON));
  def("CompMultiply",Vec3CompMultiply);
  def("CompDivide",Vec3CompDivide);
  def("Distance",Vec3Distance);  
  def("Equal",Mat3Equal, (arg("m1"), arg("m2"), arg("epsilon")=EPSILON));
  def("DihedralAngle", &DihedralAngle);
  def("SignedAngle", Vec3SignedAngle, (arg("v1"), arg("v2"),arg("ref_normal")));
  def("Dot",Mat3Dot);
  def("Det",Mat3Det);
  def("Cross",Vec3Cross);
  def("Transpose",Mat3Transpose);
  def("Invert",Mat3Invert);
  def("Comp",Mat3Comp);
  def("Minor",Mat3Minor);
  def("EulerTransformation",EulerTransformation);
  def("AxisRotation",AxisRotation);
  def("OrthogonalVector",OrthogonalVector);
  def("Min",Vec3Min);
  def("Max",Vec3Max);
  def("Distance2WithPBC",Vec3Distance2WithPBC);
  def("DistanceWithPBC",Vec3DistanceWithPBC);
  def("MinDistance",MinDistance);
  def("MinDistanceIndices",MinDistanceIndices);
  def("MinDistanceWithPBC",MinDistanceWithPBC);
  def("WrapVec3",wrap_vec3_1,(arg("vector"),arg("cell_center"),arg("cell_size")));
  def("WrapVec3",wrap_vec3_2,(arg("vector"),arg("cell_center"),arg("cell_size"),arg("cell_angles")));
  def("WrapVec3List",wrap_vec3list1,(arg("vector_list"),arg("cell_center"),arg("cell_size")));
  def("WrapVec3List",wrap_vec3list2,(arg("vector_list"),arg("cell_center"),arg("cell_size"),arg("cell_angles")));  
}
