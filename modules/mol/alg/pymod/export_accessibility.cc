//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2017 by the OpenStructure authors
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

#include <ost/mol/alg/accessibility.hh>

using namespace boost::python;

namespace{

Real WrapAccessibilityHandle(ost::mol::EntityHandle& ent, 
                             Real probe_radius, 
                             bool include_hydrogens, 
                             bool include_hetatm,
                             bool include_water,
                             const String& selection,
                             const String& asa_abs, 
                             const String& asa_rel,
                             const String& asa_atom) {
                                 
  return ost::mol::alg::Accessibility(ent, probe_radius, include_hydrogens, 
                                      include_hetatm, include_water,
                                      selection, asa_abs, asa_rel, asa_atom);
}

Real WrapAccessibilityView(ost::mol::EntityView& ent, 
                           Real probe_radius, 
                           bool include_hydrogens, 
                           bool include_hetatm,
                           bool include_water,
                           const String& selection,
                           const String& asa_abs, const String& asa_rel,
                           const String& asa_atom) {

  return ost::mol::alg::Accessibility(ent, probe_radius, include_hydrogens, 
                                      include_hetatm, include_water, 
                                      selection, asa_abs, asa_rel, asa_atom);
}

} // ns




void export_accessibility() {

    def("Accessibility", WrapAccessibilityHandle, (arg("ent"), 
                                                   arg("probe_radius")=1.4,
                                                   arg("include_hydrogens")=false,
                                                   arg("include_hetatm")=false,
                                                   arg("include_water")=false,
                                                   arg("selection")="",
                                                   arg("asa_abs")="asaAbs",
                                                   arg("asa_rel")="asaRel",
                                                   arg("asa_atom")="asaAtom"));

    def("Accessibility", WrapAccessibilityView, (arg("ent"), 
                                                 arg("probe_radius")=1.4,
                                                 arg("include_hydrogens")=false,
                                                 arg("include_hetatm")=false,
                                                 arg("include_water")=false,
                                                 arg("selection")="",
                                                 arg("asa_abs")="asaAbs",
                                                 arg("asa_rel")="asaRel",
                                                 arg("asa_atom")="asaAtom"));
}

