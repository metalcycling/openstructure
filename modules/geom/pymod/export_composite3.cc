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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <ost/geom/geom.hh>

using namespace boost::python;
void export_Composite3()
{
  using namespace geom;

  class_<Line3>("Line3",init<>())
    .def(init<const Vec3&,const Vec3&>())
    .def(self == self)
    .def("At",&Line3::At)
    .def("GetOrigin",&Line3::GetOrigin)
    .def("GetDirection",&Line3::GetDirection)
    .add_property("origin", &Line3::GetOrigin)
    .add_property("direction", &Line3::GetDirection)
    .def(self_ns::str(self))
  ;

{ // scope

scope PlaneScope =
  class_<Plane>("Plane",init<>())
    .def(init<const Vec3&,const Vec3&>())
    .def(init<const Vec3&, const Vec3&, const Vec3&>())
    .def(init<const Line3&, const Vec3&>())
    .def(init< Real, Real, Real, Real>())
    .def(init< Real, Real, Real>())
    .def(self == self)
    .def("At",&Plane::At)
    .def("GetOrigin",&Plane::GetOrigin)
    .def("GetNormal",&Plane::GetNormal)
    .add_property("normal", &Plane::GetNormal)
    .add_property("origin", &Plane::GetOrigin)
    .add_property("p", &Plane::GetP)
    .def("GetP",&Plane::GetP)
  ;

    enum_<Plane::PLANE_TYPE>("PLANE_TYPE")
      .value("YZ", Plane::YX)
      .value("ZY", Plane::ZY)
      .value("XZ", Plane::XZ)
      .value("ZX", Plane::ZX)
      .value("XY", Plane::XY)
      .value("YX", Plane::YX)
     . export_values()
    ;



} // scope

   class_<Rotation3>("Rotation3",init<>())
    .def(init<Real,Real,Real,optional<const Vec3&> >())
    .def(init<const Vec3&, Real,optional<const Vec3&> >())
    .def(init<const Line3&, Real>())
    .def(init<const Mat3&,optional<const Vec3&> >())
    .def(init<const Quat&,optional<const Vec3&> >())
    .def("GetOrigin",&Rotation3::GetOrigin)
    .def("SetOrigin",&Rotation3::SetOrigin)
    .def("GetPhi",&Rotation3::GetPhi)
    .def("SetPhi",&Rotation3::SetPhi)
    .def("GetTheta",&Rotation3::GetTheta)
    .def("SetTheta",&Rotation3::SetTheta)
    .def("GetPsi",&Rotation3::GetPsi)
    .def("SetPsi",&Rotation3::SetPsi)
    .def("GetRotationAxis",&Rotation3::GetRotationAxis)
    .def("SetRotationAxis",&Rotation3::SetRotationAxis)
    .def("GetRotationAngle",&Rotation3::GetRotationAngle)
    .def("SetRotationAngle",&Rotation3::SetRotationAngle)
    .def("GetRotationMatrix",&Rotation3::GetRotationMatrix)
    .def("SetRotationMatrix",&Rotation3::SetRotationMatrix)
    .def("GetQuaternion",&Rotation3::GetQuat)
    .def("SetQuaternion",&Rotation3::SetQuat)
    .def("Apply",&Rotation3::Apply)
    .def(self == self)
    ;

  class_<Rotation3List>("Rotation3List")
    .def(vector_indexing_suite<Rotation3List>())
    ;

  def("ExportEulerAngles",&ExportEulerAngles);
  def("ImportEulerAngles",&ImportEulerAngles);

  implicitly_convertible<std::string,boost::filesystem::path>();

  class_<Sphere>("Sphere",init<>())
    .def(init<const Vec3&,Real>())
    .def("GetOrigin",&Sphere::GetOrigin)
    .def("GetRadius",&Sphere::GetRadius)
    .add_property("origin", &Sphere::GetOrigin, &Sphere::SetOrigin)
    .add_property("radius", &Sphere::GetRadius, &Sphere::SetRadius)
  ;
  
  class_<CuboidAxis>("CuboidAxis", init<>())
    .def(init<const geom::Vec3&, Real>())
    .def("GetHalfExtent", &CuboidAxis::GetHalfExtent)
    .def("GetVector", &CuboidAxis::GetVector,
         return_value_policy<copy_const_reference>())
    .def("GetExtent", &CuboidAxis::GetExtent)
    
    .add_property("vector", make_function(&CuboidAxis::GetVector,
                  return_value_policy<copy_const_reference>()))
    .add_property("half_extent", &CuboidAxis::GetHalfExtent)
    .add_property("extent", &CuboidAxis::GetExtent)
  ;
  
  class_<Cuboid>("Cuboid", init<>())
    .def(init<const geom::Vec3&, const CuboidAxis&,const CuboidAxis&, 
              const CuboidAxis&>())
    .def("GetCenter",&Cuboid::GetCenter)
    .add_property("center", &Cuboid::GetCenter)
    .def("GetVecA", &Cuboid::GetVecA,
         return_value_policy<copy_const_reference>())
    .add_property("vec_a", make_function(&Cuboid::GetVecA,
              return_value_policy<copy_const_reference>()))
    .def("GetVecB", &Cuboid::GetVecB,
         return_value_policy<copy_const_reference>())
   .add_property("vec_b", make_function(&Cuboid::GetVecB,
             return_value_policy<copy_const_reference>()))         
    .def("GetVecC", &Cuboid::GetVecC,
         return_value_policy<copy_const_reference>())
    .add_property("vec_c", make_function(&Cuboid::GetVecC,
                  return_value_policy<copy_const_reference>()))         
    .def("GetAxisA", &Cuboid::GetAxisA, 
         return_value_policy<copy_const_reference>())
    .add_property("axis_a", make_function(&Cuboid::GetAxisA,
                 return_value_policy<copy_const_reference>()))
    .def("GetAxisB", &Cuboid::GetAxisB, 
         return_value_policy<copy_const_reference>())
    .add_property("axis_b", make_function(&Cuboid::GetAxisB,
                  return_value_policy<copy_const_reference>()))
    .def("GetAxisC", &Cuboid::GetAxisC, 
         return_value_policy<copy_const_reference>())              
    .add_property("axis_c", make_function(&Cuboid::GetAxisC,
                  return_value_policy<copy_const_reference>()))             
    .def("GetHalfExtents", &Cuboid::GetHalfExtents)
    .add_property("half_extents", &Cuboid::GetHalfExtents)
  ;
  class_<AlignedCuboid>("AlignedCuboid", init<geom::Vec3, geom::Vec3>())
    .def("GetMin", &AlignedCuboid::GetMin, 
         return_value_policy<copy_const_reference>())
    .def("GetMax", &AlignedCuboid::GetMax, 
         return_value_policy<copy_const_reference>())
    .def("GetCenter", &AlignedCuboid::GetCenter)
    .def("GetSize", &AlignedCuboid::GetSize)
    .add_property("max", make_function(&AlignedCuboid::GetMax,
         return_value_policy<copy_const_reference>()))
    .add_property("min", make_function(&AlignedCuboid::GetMin,
         return_value_policy<copy_const_reference>()))
    .add_property("center", &AlignedCuboid::GetCenter)
    .add_property("size", &AlignedCuboid::GetSize)
    .def(self_ns::str(self))
  ;
}

