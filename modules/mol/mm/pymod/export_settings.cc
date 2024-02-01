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
#include <ost/mol/mm/settings.hh>
#include <OpenMM.h> //for definition of Integrator
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

using namespace boost::python;


void export_Settings()
{

  enum_<ost::mol::mm::Platform>("Platform")
    .value("Reference", ost::mol::mm::Reference)
    .value("OpenCL", ost::mol::mm::OpenCL)
    .value("CUDA", ost::mol::mm::CUDA)
    .value("CPU", ost::mol::mm::CPU)
  ;

  enum_<ost::mol::mm::NonbondedMethod>("NonbondedMethod")
    .value("NoCutoff", ost::mol::mm::NoCutoff)
    .value("CutoffNonPeriodic", ost::mol::mm::CutoffNonPeriodic)
    .value("CutoffPeriodic", ost::mol::mm::CutoffPeriodic)
    .value("Ewald", ost::mol::mm::Ewald)
    .value("PME", ost::mol::mm::PME)
  ;

  //export the termini exceptions object
  class_<ost::mol::mm::TerminiExceptions>("TerminiExceptions",init<>())
    .def("SetException", &ost::mol::mm::TerminiExceptions::SetException)
    .def("HasException", &ost::mol::mm::TerminiExceptions::HasException)
    .def("GetException", &ost::mol::mm::TerminiExceptions::GetException)
  ;

  //export the settings

  class_<ost::mol::mm::Settings>("Settings",init<>())
    .def_readwrite("add_bonds",&ost::mol::mm::Settings::add_bonds)
    .def_readwrite("add_angles",&ost::mol::mm::Settings::add_angles)
    .def_readwrite("add_dihedrals",&ost::mol::mm::Settings::add_dihedrals)
    .def_readwrite("add_impropers",&ost::mol::mm::Settings::add_impropers)
    .def_readwrite("add_cmaps",&ost::mol::mm::Settings::add_cmaps)
    .def_readwrite("add_exclusions",&ost::mol::mm::Settings::add_exclusions)
    .def_readwrite("add_nonbonded",&ost::mol::mm::Settings::add_nonbonded)
    .def_readwrite("add_gbsa",&ost::mol::mm::Settings::add_gbsa)
    .def_readwrite("constrain_hbonds",&ost::mol::mm::Settings::constrain_hbonds)
    .def_readwrite("constrain_bonds",&ost::mol::mm::Settings::constrain_bonds)
    .def_readwrite("rigid_water",&ost::mol::mm::Settings::rigid_water)
    .def_readwrite("strict_interactions",&ost::mol::mm::Settings::strict_interactions)
    .def_readwrite("ideal_bond_length_constraints",&ost::mol::mm::Settings::ideal_bond_length_constraints)
    .def_readwrite("fix_heavy_atoms",&ost::mol::mm::Settings::fix_heavy_atoms)
    .def_readwrite("kill_electrostatics",&ost::mol::mm::Settings::kill_electrostatics)
    .def_readwrite("generate_disulfid_bonds",&ost::mol::mm::Settings::generate_disulfid_bonds)
    .def_readwrite("nonbonded_method",&ost::mol::mm::Settings::nonbonded_method)
    .def_readwrite("nonbonded_cutoff",&ost::mol::mm::Settings::nonbonded_cutoff)
    .def_readwrite("remove_cmm_motion",&ost::mol::mm::Settings::remove_cmm_motion)
    .def_readwrite("periodic_box_extents",&ost::mol::mm::Settings::periodic_box_extents)
    .def_readwrite("init_temperature",&ost::mol::mm::Settings::init_temperature)
    .def_readwrite("forcefield",&ost::mol::mm::Settings::forcefield)
    .def_readwrite("termini_exceptions",&ost::mol::mm::Settings::termini_exceptions)
    .def_readwrite("platform",&ost::mol::mm::Settings::platform)
    .def_readwrite("reference_properties",&ost::mol::mm::Settings::reference_properties)
    .def_readwrite("cpu_properties",&ost::mol::mm::Settings::cpu_properties)
    .def_readwrite("cuda_properties",&ost::mol::mm::Settings::cuda_properties)
    .def_readwrite("opencl_properties",&ost::mol::mm::Settings::opencl_properties)
    .def_readwrite("add_thermostat",&ost::mol::mm::Settings::add_thermostat)
    .def_readwrite("thermostat_temperature",&ost::mol::mm::Settings::thermostat_temperature)
    .def_readwrite("thermostat_collision_frequency",&ost::mol::mm::Settings::thermostat_collision_frequency)
    .def_readwrite("add_barostat",&ost::mol::mm::Settings::add_barostat)
    .def_readwrite("barostat_temperature",&ost::mol::mm::Settings::barostat_temperature)
    .def_readwrite("barostat_pressure",&ost::mol::mm::Settings::barostat_pressure)
    .def_readwrite("barostat_frequency",&ost::mol::mm::Settings::barostat_frequency)
    .def_readwrite("integrator",&ost::mol::mm::Settings::integrator)
    .def_readwrite("solvent_dielectric",&ost::mol::mm::Settings::solvent_dielectric)
    .def_readwrite("solute_dielectric",&ost::mol::mm::Settings::solute_dielectric)
    .def_readwrite("reaction_field_dielectric",&ost::mol::mm::Settings::reaction_field_dielectric)
    .def_readwrite("use_dispersion_correction",&ost::mol::mm::Settings::use_dispersion_correction)
    .def_readwrite("keep_ff_specific_naming",&ost::mol::mm::Settings::keep_ff_specific_naming)
    .def_readwrite("openmm_plugin_directory",&ost::mol::mm::Settings::openmm_plugin_directory)
    .def_readwrite("custom_plugin_directory",&ost::mol::mm::Settings::custom_plugin_directory)

  ;

  class_<ost::mol::mm::PropertyMap>("PropertyMap", no_init)
    .def(map_indexing_suite<ost::mol::mm::PropertyMap>()) 
  ;

  boost::python::register_ptr_to_python<ost::mol::mm::SettingsPtr>();
  boost::python::register_ptr_to_python<ost::mol::mm::TerminiExceptionsPtr>();

}

