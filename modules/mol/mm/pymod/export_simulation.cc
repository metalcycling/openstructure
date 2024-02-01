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
#include <ost/mol/mm/simulation.hh>

using namespace boost::python;



namespace{


template<typename T>
std::vector<T> ListToVec(const boost::python::list& l){
  std::vector<T> vec;

  for (int i = 0; i < boost::python::len(l); ++i){
    vec.push_back(boost::python::extract<T>(l[i]));
  }
  return vec;
}

void WrapAddPositionConstraints(ost::mol::mm::SimulationPtr p, const boost::python::list& l){
  std::vector<uint> vec  = ListToVec<uint>(l);
  p->AddPositionConstraints(vec);
}

}


void export_Simulation()
{

  class_<ost::mol::mm::Simulation>("Simulation",no_init)
    .def(init<const ost::mol::EntityHandle, const ost::mol::mm::SettingsPtr>())
    .def(init<const ost::mol::mm::TopologyPtr,const ost::mol::EntityHandle&,const ost::mol::mm::SettingsPtr>())
    .def("Save",&ost::mol::mm::Simulation::Save,(arg("filename")))
    .def("Load",&ost::mol::mm::Simulation::Load,(arg("filename")))
    .staticmethod("Load")
    .def("IsPlatformAvailable",&ost::mol::mm::Simulation::IsPlatformAvailable,
         (arg("settings")))
    .staticmethod("IsPlatformAvailable")
    .def("Steps",&ost::mol::mm::Simulation::Steps,(arg("steps")))
    .def("GetPositions",&ost::mol::mm::Simulation::GetPositions,(arg("enforce_periodic_box")=false, arg("in_angstrom")=true))
    .def("GetVelocities",&ost::mol::mm::Simulation::GetVelocities)
    .def("GetForces",&ost::mol::mm::Simulation::GetForces)
    .def("SetPositions",&ost::mol::mm::Simulation::SetPositions,(arg("positions"),arg("in_angstrom")=true))
    .def("SetVelocities",&ost::mol::mm::Simulation::SetVelocities)
    .def("GetEntity",&ost::mol::mm::Simulation::GetEntity)
    .def("ApplySD",&ost::mol::mm::Simulation::ApplySD, (arg("tolerance")=100.0,arg("max_iterations")=1000))
    .def("ApplyLBFGS",&ost::mol::mm::Simulation::ApplyLBFGS,(arg("tolerance")=1.0,arg("max_iterations")=1000))
    .def("GetEnergy",&ost::mol::mm::Simulation::GetEnergy)                                                                             
    .def("GetPotentialEnergy",&ost::mol::mm::Simulation::GetPotentialEnergy)                                                                             
    .def("GetKineticEnergy",&ost::mol::mm::Simulation::GetKineticEnergy) 
    .def("Register",&ost::mol::mm::Simulation::Register) 
    .def("UpdatePositions",&ost::mol::mm::Simulation::UpdatePositions, (arg("enforce_periodic_box")=false))
    .def("GetTopology",&ost::mol::mm::Simulation::GetTopology)
    .def("ResetHarmonicBond",&ost::mol::mm::Simulation::ResetHarmonicBond,(arg("bond_idx"),arg("bond_length"),arg("force_constant"))) 
    .def("ResetHarmonicAngle",&ost::mol::mm::Simulation::ResetHarmonicAngle,(arg("angle_idx"),arg("angle"),arg("force_constant")))
    .def("ResetUreyBradleyAngle",&ost::mol::mm::Simulation::ResetUreyBradleyAngle,(arg("angle_idx"),arg("angle"),arg("angle_force_constant"),arg("bond_length"),arg("bond_force_constant")))
    .def("ResetPeriodicDihedral",&ost::mol::mm::Simulation::ResetPeriodicDihedral,(arg("dihedral_idx"),arg("multiplicity"),arg("phase"),arg("force_constant")))
    .def("ResetPeriodicImproper",&ost::mol::mm::Simulation::ResetPeriodicImproper,(arg("improper_idx"),arg("multiplicity"),arg("phase"),arg("force_constant")))
    .def("ResetHarmonicImproper",&ost::mol::mm::Simulation::ResetHarmonicImproper,(arg("improper_idx"),arg("phase"),arg("force_constant")))
    .def("ResetLJPair",&ost::mol::mm::Simulation::ResetLJPair,(arg("ljpair_idx"),arg("sigma"),arg("epsilon")))
    .def("ResetDistanceConstraint",&ost::mol::mm::Simulation::ResetDistanceConstraint,(arg("constraint_idx"),arg("length")))  
    .def("ResetHarmonicPositionRestraint",&ost::mol::mm::Simulation::ResetHarmonicPositionRestraint,(arg("restraint_idx"),arg("ref_position"),arg("k"),arg("x_scale")=1.0,arg("y_scale")=1.0,arg("z_scale")=1.0))                                                                     
    .def("ResetHarmonicDistanceRestraint",&ost::mol::mm::Simulation::ResetHarmonicDistanceRestraint,(arg("restraint_idx"),arg("length"),arg("force_constant")))
    .def("ResetLJ",&ost::mol::mm::Simulation::ResetLJ,(arg("particle_idx"),arg("sigma"),arg("epsilon")))
    .def("ResetGBSA",&ost::mol::mm::Simulation::ResetGBSA,(arg("particle_idx"),arg("radius"),arg("scaling")))
    .def("ResetCharge",&ost::mol::mm::Simulation::ResetCharge,(arg("particle_idx"),arg("charge")))
    .def("ResetMass",&ost::mol::mm::Simulation::ResetMass,(arg("particle_idx"),arg("mass")))
    .def("AddPositionConstraint",&ost::mol::mm::Simulation::AddPositionConstraint,(arg("particle_idx")))
    .def("AddPositionConstraints",&WrapAddPositionConstraints,(arg("particle_idx")))
    .def("ResetPositionConstraints",&ost::mol::mm::Simulation::ResetPositionConstraints)
    .def("GetPeriodicBoxExtents",&ost::mol::mm::Simulation::GetPeriodicBoxExtents)
    .def("SetPeriodicBoxExtents",&ost::mol::mm::Simulation::SetPeriodicBoxExtents,(arg("extents")))
  ;

  boost::python::register_ptr_to_python<ost::mol::mm::SimulationPtr>();
}
