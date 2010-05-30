//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
#include <ost/mol/transform.hh>
#include <ost/info/info.hh>

using namespace boost::python;
using namespace ost::mol;

void export_Entity();
void export_EntityView();
void export_ChainView();
void export_Surface();
void export_Query();
void export_Chain();
void export_Residue();
void export_Atom();
void export_Bond();
void export_Torsion();
void export_Visitor();
void export_AtomView();
void export_ResidueView();
void export_Editors();
void export_CoordGroup();
void export_PropertyID();
void export_BoundingBox();
void export_QueryViewWrapper();
BOOST_PYTHON_MODULE(_mol)
{
  export_Entity();
  export_Surface();
  export_Query();
  export_Chain();
  export_Residue();
  export_Atom();
  export_Bond();
  export_Visitor();
  export_Torsion();
  export_AtomView();
  export_ResidueView();
  export_ChainView();  
  export_EntityView(); 
  export_Editors();
  export_CoordGroup();
  export_PropertyID();  
  export_BoundingBox();
  export_QueryViewWrapper();

  class_<Transform>("Transform", init<>())
    .def(init<const Transform&>()) // shouldn't this be there automatically ?
    .def("GetMatrix",&Transform::GetMatrix)
    .def("GetTransposedMatrix",&Transform::GetTransposedMatrix)
    .def("SetTrans",&Transform::SetTrans)
    .def("GetTrans",&Transform::GetTrans)
    .def("SetCenter",&Transform::SetCenter)
    .def("GetCenter",&Transform::GetCenter)
    .def("SetRot",&Transform::SetRot)
    .def("GetRot",&Transform::GetRot)
    .def("ApplyXAxisRotation",&Transform::ApplyXAxisRotation)
    .def("ApplyYAxisRotation",&Transform::ApplyYAxisRotation)
    .def("ApplyZAxisRotation",&Transform::ApplyZAxisRotation)
    .def("ApplyXAxisTranslation",&Transform::ApplyXAxisTranslation)
    .def("ApplyYAxisTranslation",&Transform::ApplyYAxisTranslation)
    .def("ApplyZAxisTranslation",&Transform::ApplyZAxisTranslation)
    .def("ApplyAxisRotation",&Transform::ApplyAxisRotation)
    ;
  def("TransformToInfo", &TransformToInfo);
  def("TransformFromInfo", &TransformFromInfo);

}
