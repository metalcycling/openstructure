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
  Author: Ansgar Philippsen
*/

#include <boost/python.hpp>
#include <ost/units.hh>

using namespace boost::python;
using namespace ost;

#ifdef OST_STATIC_PROPERTY_WORKAROUND

struct Units_ {
    Real get_m() { return Units::m; }
    Real get_mm() { return Units::mm; }
    Real get_nm() { return Units::nm; }
    Real get_um() { return Units::um; }
    Real get_pm() { return Units::pm; }
    Real get_A() { return Units::A; }
    Real get_deg() { return Units::deg; }
    Real get_rad() { return Units::rad; }
    Real get_mrad() { return Units::mrad; }
    Real get_J() { return Units::J; }
    Real get_eV() { return Units::eV; }
    Real get_V() { return Units::V; }
    Real get_kV() { return Units::kV; }
    Real get_s() { return Units::s; }
    Real get_h() { return Units::h; }
    Real get_kg() { return Units::kg; }
    Real get_g() { return Units::g; }
    Real get_mg() { return Units::mg; }
    Real get_kDa() { return Units::kDa; }
    Real get_Da() { return Units::Da; }
    Real get_C() { return Units::C; }
};


#endif

void export_Units()
{

#ifdef OST_STATIC_PROPERTY_WORKAROUND
  object units=class_<Units_>("Units_")
    .add_property("m",&Units_::get_m)
    .add_property("mm",&Units_::get_mm)
    .add_property("nm",&Units_::get_nm)
    .add_property("um",&Units_::get_um)
    .add_property("pm",&Units_::get_pm)
    .add_property("A",&Units_::get_A)
    .add_property("deg",&Units_::get_deg)
    .add_property("rad",&Units_::get_rad)
    .add_property("mrad",&Units_::get_mrad)
    .add_property("J",&Units_::get_J)
    .add_property("eV",&Units_::get_eV)
    .add_property("V",&Units_::get_V)
    .add_property("kV",&Units_::get_kV)
    .add_property("s",&Units_::get_s)
    .add_property("h",&Units_::get_h)
    .add_property("kg",&Units_::get_kg)
    .add_property("g",&Units_::get_g)
    .add_property("mg",&Units_::get_mg)
    .add_property("kDa",&Units_::get_kDa)
    .add_property("Da",&Units_::get_Da)
    .add_property("C",&Units_::get_C)
  ;
  scope().attr("Units")=units();
#else
  class_<Units>("Units")
    .def_readonly("m",Units::m)
    .def_readonly("mm",Units::mm)
    .def_readonly("nm",Units::nm)
    .def_readonly("um",Units::um)
    .def_readonly("pm",Units::pm)
    .def_readonly("A",Units::A)
    .def_readonly("deg",Units::deg)
    .def_readonly("rad",Units::rad)
    .def_readonly("mrad",Units::mrad)
    .def_readonly("J",Units::J)
    .def_readonly("eV",Units::eV)
    .def_readonly("V",Units::V)
    .def_readonly("kV",Units::kV)
    .def_readonly("s",Units::s)
    .def_readonly("h",Units::h)
    .def_readonly("kg",Units::kg)
    .def_readonly("g",Units::g)
    .def_readonly("mg",Units::mg)
    .def_readonly("kDa",Units::kDa)
    .def_readonly("Da",Units::Da)
    .def_readonly("C",Units::C)
  ;
#endif
}
