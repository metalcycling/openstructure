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
#include <stdexcept>
#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>

using namespace boost::python;

#include <ost/mol/alg/molck.hh>

using namespace ost::mol::alg;

namespace {

object MolckSettingsInitWrapper(tuple args, dict kwargs){

  object self = args[0];
  args = tuple(args.slice(1,_));

  bool rm_unk_atoms = false;
  if(kwargs.contains("rm_unk_atoms")){
    rm_unk_atoms = extract<bool>(kwargs["rm_unk_atoms"]);
    kwargs["rm_unk_atoms"].del();
  }

  bool rm_non_std = false;
  if(kwargs.contains("rm_non_std")){
    rm_non_std = extract<bool>(kwargs["rm_non_std"]);
    kwargs["rm_non_std"].del();
  }

  bool rm_hyd_atoms = true;
  if(kwargs.contains("rm_hyd_atoms")){
    rm_hyd_atoms = extract<bool>(kwargs["rm_hyd_atoms"]);
    kwargs["rm_hyd_atoms"].del();
  }

  bool rm_oxt_atoms = false;
  if(kwargs.contains("rm_oxt_atoms")){
    rm_oxt_atoms = extract<bool>(kwargs["rm_oxt_atoms"]);
    kwargs["rm_oxt_atoms"].del();
  }

  bool rm_zero_occ_atoms = false;
  if(kwargs.contains("rm_zero_occ_atoms")){
    rm_zero_occ_atoms = extract<bool>(kwargs["rm_zero_occ_atoms"]);
    kwargs["rm_zero_occ_atoms"].del();
  }

  bool colored = false;
  if(kwargs.contains("colored")){
    colored = extract<bool>(kwargs["colored"]);
    kwargs["colored"].del();
  }

  bool map_nonstd_res = true;
  if(kwargs.contains("map_nonstd_res")){
    map_nonstd_res = extract<bool>(kwargs["map_nonstd_res"]);
    kwargs["map_nonstd_res"].del();
  }

  bool assign_elem = true;
  if(kwargs.contains("assign_elem")){
    assign_elem = extract<bool>(kwargs["assign_elem"]);
    kwargs["assign_elem"].del();
  }

  if(len(kwargs) > 0){
    std::stringstream ss;
    ss << "Invalid keywords observed when setting up MolckSettings! ";
    ss << "Or did you pass the same keyword twice? ";
    ss << "Valid keywords are: rm_unk_atoms, rm_non_std, rm_hyd_atoms, ";
    ss << "rm_oxt_atoms, rm_zero_occ_atoms, colored, map_nonstd_res, ";
    ss << "assign_elem!";
    throw std::invalid_argument(ss.str());
  }


  return self.attr("__init__")(rm_unk_atoms,
                               rm_non_std,
                               rm_hyd_atoms,
                               rm_oxt_atoms,
                               rm_zero_occ_atoms,
                               colored,
                               map_nonstd_res,
                               assign_elem);
}}

void export_Molck()
{
  class_<MolckSettings>("MolckSettings", no_init)
    .def("__init__", raw_function(MolckSettingsInitWrapper))
    .def(init<bool , bool, bool, bool, bool, bool, bool, bool>())
    .def("ToString", &MolckSettings::ToString)
    .def("__repr__", &MolckSettings::ToString)
    .def("__str__", &MolckSettings::ToString)
    .def_readwrite("rm_unk_atoms", &MolckSettings::rm_unk_atoms)
    .def_readwrite("rm_non_std", &MolckSettings::rm_non_std)
    .def_readwrite("rm_hyd_atoms", &MolckSettings::rm_hyd_atoms)
    .def_readwrite("rm_oxt_atoms", &MolckSettings::rm_oxt_atoms)
    .def_readwrite("rm_zero_occ_atoms", &MolckSettings::rm_zero_occ_atoms)
    .def_readwrite("colored", &MolckSettings::colored)
    .def_readwrite("map_nonstd_res", &MolckSettings::map_nonstd_res)
    .def_readwrite("assign_elem", &MolckSettings::assign_elem);

  def("MapNonStandardResidues", &MapNonStandardResidues, (arg("ent"),
                                                          arg("lib"),
                                                          arg("log_diags")=false));

  def("RemoveAtoms", &RemoveAtoms, (arg("ent"),
                                    arg("lib"),
                                    arg("rm_unk_atoms")=false,
                                    arg("rm_non_std")=false,
                                    arg("rm_hyd_atoms")=true,
                                    arg("rm_oxt_atoms")=false,
                                    arg("rm_zero_occ_atoms")=false,
                                    arg("colored")=false));

  def("CleanUpElementColumn", &CleanUpElementColumn, (arg("ent"), arg("lib")));

  def("Molck", &Molck, (arg("ent"), arg("lib"), arg("settings"), arg("prune")=true));
}
