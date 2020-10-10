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

#include <ost/mol/alg/find_membrane.hh>

using namespace boost::python;

namespace{
  ost::mol::alg::FindMemParam FindMembraneView(ost::mol::EntityView& v,
                                               bool assign_membrane_representation, 
                                               bool fast) {
    return ost::mol::alg::FindMembrane(v, assign_membrane_representation, fast);
  }

  ost::mol::alg::FindMemParam FindMembraneHandle(ost::mol::EntityHandle& h,
                                                 bool assign_membrane_representation, 
                                                 bool fast) {
    return ost::mol::alg::FindMembrane(h, assign_membrane_representation, fast);
  }
}

void export_find_membrane() {

  class_<ost::mol::alg::FindMemParam>("FindMemParam", no_init)
    .def_readonly("tilt", &ost::mol::alg::FindMemParam::tilt)
    .def_readonly("angle", &ost::mol::alg::FindMemParam::angle)
    .def_readonly("width", &ost::mol::alg::FindMemParam::width)
    .def_readonly("pos", &ost::mol::alg::FindMemParam::pos)
    .def_readonly("energy", &ost::mol::alg::FindMemParam::energy)
    .def_readonly("membrane_asa", &ost::mol::alg::FindMemParam::membrane_asa)
    .def_readonly("axis", &ost::mol::alg::FindMemParam::axis)
    .def_readonly("tilt_axis", &ost::mol::alg::FindMemParam::tilt_axis)
    .def_readonly("membrane_axis", &ost::mol::alg::FindMemParam::GetMembraneAxis)
    .def_readonly("membrane_representation", &ost::mol::alg::FindMemParam::membrane_representation)
  ;
  
  def("FindMembrane", FindMembraneView, (arg("ent"), arg("assign_membrane_representation")=true, 
                                         arg("fast")=false));
  def("FindMembrane", FindMembraneHandle, (arg("ent"), arg("assign_membrane_representation")=true,
                                           arg("fast")=false));
}
