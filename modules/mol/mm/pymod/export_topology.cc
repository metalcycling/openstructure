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
#include <ost/mol/mm/topology.hh>
#include <ost/mol/mm/topology_creator.hh>
#include <ost/mol/entity_handle.hh>
#include <vector>

using namespace boost::python;
namespace{

  template<typename T>
  int Get(T index, int i){
    return index[i];
  }

  template<typename T>
  void Set(T& index, int i, int v){
    index[i] = v;
  }

  template<typename T>
  std::vector<T> ListToVec(const boost::python::list& l){

    std::vector<T> vec;
    for (int i = 0; i < boost::python::len(l); ++i){
      vec.push_back(boost::python::extract<T>(l[i]));
    }
    return vec;
  }

  template<typename T>
  boost::python::list VecToList(std::vector<T>& vec){
    boost::python::list l;
    for(typename std::vector<T>::iterator it=vec.begin();it!=vec.end();++it){
      l.append(*it);
    }
    return l;
  }

  ost::mol::mm::TopologyPtr WrapTopologyConstructor(const boost::python::list& masses_list){
    std::vector<Real> masses_vector = ListToVec<Real>(masses_list);
    ost::mol::mm::TopologyPtr p(new ost::mol::mm::Topology(masses_vector));
    return p;
  }

  void WrapAddCMap(ost::mol::mm::TopologyPtr p, int index_one, int index_two, int index_three,
                   int index_four, int index_five, uint dimension, const boost::python::list& l){
    std::vector<Real> v = ListToVec<Real>(l);
    p->AddCMap(index_one,index_two,index_three,index_four,index_five,dimension,v);
  }

  void WrapSetSigmas(ost::mol::mm::TopologyPtr p, const boost::python::list& l){
    std::vector<Real> v = ListToVec<Real>(l);
    p->SetSigmas(v);
  }

  void WrapSetEpsilons(ost::mol::mm::TopologyPtr p, const boost::python::list& l){
    std::vector<Real> v = ListToVec<Real>(l);
    p->SetEpsilons(v);
  }

  void WrapSetGBSARadii(ost::mol::mm::TopologyPtr p, const boost::python::list& l){
    std::vector<Real> v = ListToVec<Real>(l);
    p->SetGBSARadii(v);
  }

  void WrapSetOBCScalings(ost::mol::mm::TopologyPtr p, const boost::python::list& l){
    std::vector<Real> v = ListToVec<Real>(l);
    p->SetOBCScalings(v);
  }

  void WrapSetCharges(ost::mol::mm::TopologyPtr p, const boost::python::list& l){
    std::vector<Real> v = ListToVec<Real>(l);
    p->SetCharges(v);
  }

  void WrapSetMasses(ost::mol::mm::TopologyPtr p, const boost::python::list& l){
    std::vector<Real> v = ListToVec<Real>(l);
    p->SetMasses(v);
  }


  boost::python::list WrapGetSigmas(ost::mol::mm::TopologyPtr p){
    std::vector<Real> v = p->GetSigmas();
    boost::python::list l = VecToList<Real>(v);
    return l;
  }

  boost::python::list WrapGetEpsilons(ost::mol::mm::TopologyPtr p){
    std::vector<Real> v = p->GetEpsilons();
    boost::python::list l = VecToList<Real>(v);
    return l;
  }

  boost::python::list WrapGetGBSARadii(ost::mol::mm::TopologyPtr p){
    std::vector<Real> v = p->GetGBSARadii();
    boost::python::list l = VecToList<Real>(v);
    return l;
  }

  boost::python::list WrapGetOBCScalings(ost::mol::mm::TopologyPtr p){
    std::vector<Real> v = p->GetOBCScalings();
    boost::python::list l = VecToList<Real>(v);
    return l;
  }

  boost::python::list WrapGetCharges(ost::mol::mm::TopologyPtr p){
    std::vector<Real> v = p->GetCharges();
    boost::python::list l = VecToList<Real>(v);
    return l;
  }

  boost::python::list WrapGetMasses(ost::mol::mm::TopologyPtr p){
    std::vector<Real> v = p->GetMasses();
    boost::python::list l = VecToList<Real>(v);
    return l;
  }


  boost::python::tuple WrapGetHarmonicBondParam(ost::mol::mm::TopologyPtr p,uint index){
    uint i1,i2;
    Real b,k;
    p->GetHarmonicBondParameters(index,i1,i2,b,k);
    return boost::python::make_tuple(i1,i2,b,k);
  }

  boost::python::tuple WrapGetHarmonicAngleParam(ost::mol::mm::TopologyPtr p,uint index){
    uint i1,i2,i3;
    Real a,k;
    p->GetHarmonicAngleParameters(index,i1,i2,i3,a,k);
    return boost::python::make_tuple(i1,i2,i3,a,k);
  }

  boost::python::tuple WrapGetUreyBradleyAngleParam(ost::mol::mm::TopologyPtr p,uint index){
    uint i1,i2,i3;
    Real a,a_k,b,b_k;
    p->GetUreyBradleyAngleParameters(index,i1,i2,i3,a,a_k,b,b_k);
    return boost::python::make_tuple(i1,i2,i3,a,a_k,b,b_k);
  }

  boost::python::tuple WrapGetPeriodicDihedralParam(ost::mol::mm::TopologyPtr p,uint index){
    uint i1,i2,i3,i4;
    int m;
    Real ph,k;
    p->GetPeriodicDihedralParameters(index,i1,i2,i3,i4,m,ph,k);
    return boost::python::make_tuple(i1,i2,i3,i4,m,ph,k);
  }

  boost::python::tuple WrapGetPeriodicImproperParam(ost::mol::mm::TopologyPtr p,uint index){
    uint i1,i2,i3,i4;
    int m;
    Real ph,k;
    p->GetPeriodicImproperParameters(index,i1,i2,i3,i4,m,ph,k);
    return boost::python::make_tuple(i1,i2,i3,i4,m,ph,k);
  }

  boost::python::tuple WrapGetHarmonicImproperParam(ost::mol::mm::TopologyPtr p,uint index){
    uint i1,i2,i3,i4;
    Real a,k;
    p->GetHarmonicImproperParameters(index,i1,i2,i3,i4,a,k);
    return boost::python::make_tuple(i1,i2,i3,i4,a,k);
  }

  boost::python::tuple WrapGetCMapParam(ost::mol::mm::TopologyPtr p, uint index){
    std::vector<Real> map;
    int dimension;
    uint i1,i2,i3,i4,i5;
    p->GetCMapParameters(index,i1,i2,i3,i4,i5,dimension,map);
    boost::python::list l = VecToList<Real>(map);
    return boost::python::make_tuple(i1,i2,i3,i4,i5,dimension,l);
  }

  boost::python::tuple WrapGetLJPairParam(ost::mol::mm::TopologyPtr p, uint index){
    uint i1,i2;
    Real s,e;
    p->GetLJPairParameters(index,i1,i2,s,e);
    return boost::python::make_tuple(i1,i2,s,e);
  }

  boost::python::tuple WrapGetDistanceConstraintParam(ost::mol::mm::TopologyPtr p,uint index){
    uint i1,i2;
    Real d;
    p->GetDistanceConstraintParameters(index,i1,i2,d);
    return boost::python::make_tuple(i1,i2,d);
  }

  boost::python::tuple WrapGetHarmonicPositionRestraintParam(ost::mol::mm::TopologyPtr p,uint index){
    uint atom_index;
    geom::Vec3 ref_pos;
    Real k;
    Real x_scale, y_scale, z_scale;
    p->GetHarmonicPositionRestraintParameters(index,atom_index,ref_pos,k,x_scale,y_scale,z_scale);
    return boost::python::make_tuple(atom_index,ref_pos,k,x_scale,y_scale,z_scale);
  }

  boost::python::tuple WrapGetHarmonicDistanceRestraintParam(ost::mol::mm::TopologyPtr p,uint index){
    uint i1,i2;
    Real l,k;
    p->GetHarmonicDistanceRestraintParameters(index,i1,i2,l,k);
    return boost::python::make_tuple(i1,i2,l,k);
  }

  boost::python::tuple WrapGetFGMDHBondDonorParam(ost::mol::mm::TopologyPtr p, uint index){
    uint i1, i2;
    Real d,k_d,a,k_a,b,k_b;
    p->GetFGMDHBondDonorParameters(index,i1,i2,d,k_d,a,k_a,b,k_b);
    return boost::python::make_tuple(i1,i2,d,k_d,a,k_a,b,k_b);
  }

  boost::python::tuple WrapGetFGMDHBondAcceptorParam(ost::mol::mm::TopologyPtr p, uint index){
    uint i1, i2;
    p->GetFGMDHBondAcceptorParameters(index,i1,i2);
    return boost::python::make_tuple(i1,i2);
  }

  void WrapSetCMapParameters(ost::mol::mm::TopologyPtr p, uint index, int dimension, boost::python::list l){
    std::vector<Real> v = ListToVec<Real>(l);
    p->SetCMapParameters(index,dimension,v);
  }

  boost::python::list GetHarmonicBondIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2){
    std::vector<uint> return_vec = p->GetHarmonicBondIndices(i1,i2);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetHarmonicBondIndicesSingleIndex(ost::mol::mm::TopologyPtr p, int i){
    std::vector<uint> return_vec = p->GetHarmonicBondIndices(i); 
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetHarmonicAngleIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2, uint i3){
    std::vector<uint> return_vec = p->GetHarmonicAngleIndices(i1,i2,i3);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetHarmonicAngleIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetHarmonicAngleIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetUreyBradleyAngleIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2, uint i3){
    std::vector<uint> return_vec = p->GetUreyBradleyAngleIndices(i1,i2,i3);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetUreyBradleyAngleIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetUreyBradleyAngleIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetPeriodicDihedralIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2, uint i3, int i4){
    std::vector<uint> return_vec = p->GetPeriodicDihedralIndices(i1,i2,i3,i4);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetPeriodicDihedralIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetPeriodicDihedralIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetPeriodicImproperIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2, uint i3, uint i4){
    std::vector<uint> return_vec = p->GetPeriodicImproperIndices(i1,i2,i3,i4);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetPeriodicImproperIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetPeriodicImproperIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetHarmonicImproperIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2, uint i3, uint i4){
    std::vector<uint> return_vec = p->GetHarmonicImproperIndices(i1,i2,i3,i4);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetHarmonicImproperIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetHarmonicImproperIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetCMapIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2, uint i3, uint i4, uint i5){
    std::vector<uint> return_vec = p->GetCMapIndices(i1,i2,i3,i4,i5);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetCMapIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetCMapIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetLJPairIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetLJPairIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetDistanceConstraintIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetDistanceConstraintIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetHarmonicDistanceRestraintIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2){
    std::vector<uint> return_vec = p->GetHarmonicDistanceRestraintIndices(i1,i2);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetHarmonicDistanceRestraintIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetHarmonicDistanceRestraintIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetFGMDHBondDonorIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2){
    std::vector<uint> return_vec = p->GetFGMDHBondDonorIndices(i1,i2);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetFGMDHBondDonorIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetFGMDHBondDonorIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetHarmonicPositionRestraintIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
  std::vector<uint> return_vec = p->GetHarmonicPositionRestraintIndices(i);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetFGMDHBondAcceptorIndices(ost::mol::mm::TopologyPtr p, uint i1, uint i2){
    std::vector<uint> return_vec = p->GetFGMDHBondAcceptorIndices(i1,i2);
    return VecToList<uint>(return_vec);
  }

  boost::python::list GetFGMDHBondAcceptorIndicesSingleIndex(ost::mol::mm::TopologyPtr p, uint i){
    std::vector<uint> return_vec = p->GetFGMDHBondAcceptorIndices(i);
    return VecToList<uint>(return_vec);
  }


  void MergeTop(ost::mol::mm::TopologyPtr top, ost::mol::mm::TopologyPtr other){
    top->Merge(other);
  }

  void MergeTopEnt(ost::mol::mm::TopologyPtr top, ost::mol::EntityHandle& ent, 
                   ost::mol::mm::TopologyPtr other, ost::mol::EntityHandle& other_ent){
    top->Merge(ent,other,other_ent);
  }
}

void export_Topology()
{

  class_<ost::mol::mm::TopologyCreator>("TopologyCreator",no_init)
    .def("Create",&ost::mol::mm::TopologyCreator::Create).staticmethod("Create")
  ;

  class_<ost::mol::mm::Topology>("Topology",no_init)

    .def("__init__",make_constructor(&WrapTopologyConstructor))
    .def("Save",&ost::mol::mm::Topology::Save,(arg("filename")))
    .def("Load",&ost::mol::mm::Topology::Load,(arg("filename"))).staticmethod("Load")
 
    //interaction adding functions
    .def("AddHarmonicBond",&ost::mol::mm::Topology::AddHarmonicBond,(arg("idx_one"),arg("idx_two"),arg("bond_length"),arg("force_constant")))
    .def("AddHarmonicAngle",&ost::mol::mm::Topology::AddHarmonicAngle,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("angle"),arg("force_constant")))
    .def("AddUreyBradleyAngle",&ost::mol::mm::Topology::AddUreyBradleyAngle,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("angle"),arg("angle_force_constant"),arg("bond_length"),arg("bond_force_constant")))
    .def("AddPeriodicDihedral",&ost::mol::mm::Topology::AddPeriodicDihedral,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("idx_four"),arg("multiplicity"),arg("phase"),arg("force_constant")))
    .def("AddPeriodicImproper",&ost::mol::mm::Topology::AddPeriodicImproper,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("idx_four"),arg("multiplicity"),arg("phase"),arg("force_constant")))
    .def("AddHarmonicImproper",&ost::mol::mm::Topology::AddHarmonicImproper,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("idx_four"),arg("angle"),arg("force_constant")))
    .def("AddCMap",&WrapAddCMap,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("idx_four"),arg("idx_five"),arg("dimension"),arg("map")))
    .def("AddLJPair",&ost::mol::mm::Topology::AddLJPair,(arg("idx_one"),arg("idx_two"),arg("sigma"),arg("epsilon")))
    .def("AddDistanceConstraint",&ost::mol::mm::Topology::AddDistanceConstraint,(arg("idx_one"),arg("idx_two"),arg("distance")))
    .def("AddExclusion",&ost::mol::mm::Topology::AddExclusion,(arg("idx_one"),arg("idx_two")))
    .def("AddPositionConstraint",&ost::mol::mm::Topology::AddPositionConstraint,(arg("idx")))
    .def("ResetPositionConstraints",&ost::mol::mm::Topology::ResetPositionConstraints)
    .def("ResetExclusions",&ost::mol::mm::Topology::ResetExclusions)
    .def("AddHarmonicPositionRestraint",&ost::mol::mm::Topology::AddHarmonicPositionRestraint,(arg("idx"),arg("ref_position"),arg("k"),arg("x_scale")=1.0,arg("y_scale")=1.0,arg("z_scale")=1.0))
    .def("AddHarmonicDistanceRestraint",&ost::mol::mm::Topology::AddHarmonicDistanceRestraint,(arg("idx_one"),arg("idx_two"),arg("length"),arg("force_constant")))
    .def("AddFGMDHBondDonor",&ost::mol::mm::Topology::AddFGMDHBondDonor,(arg("idx_one"),arg("idx_two"),arg("length"),arg("k_length"),arg("alpha"),arg("k_alpha"),arg("beta"),arg("k_beta")))
    .def("AddFGMDHBondAcceptor",&ost::mol::mm::Topology::AddFGMDHBondAcceptor,(arg("idx_one"),arg("idx_two")))
 
    //single atom parameter getter and setter functions
    .def("SetSigmas",&WrapSetSigmas,(arg("sigmas")))
    .def("SetSigma",&ost::mol::mm::Topology::SetSigma,(arg("particle_idx"),arg("sigma")))
    .def("SetEpsilons",&WrapSetEpsilons,(arg("epsilons")))
    .def("SetEpsilon",&ost::mol::mm::Topology::SetEpsilon,(arg("particle_idx"),arg("epsilon")))
    .def("SetGBSARadii",&WrapSetGBSARadii,(arg("gbsa_radii")))
    .def("SetGBSARadius",&ost::mol::mm::Topology::SetGBSARadius,(arg("particle_idx"),arg("radius")))
    .def("SetOBCScalings",&WrapSetOBCScalings,(arg("obc_scalings")))
    .def("SetOBCScaling",&ost::mol::mm::Topology::SetOBCScaling,(arg("particle_idx"),arg("obc_scaling")))
    .def("SetCharges",&WrapSetCharges,(arg("charges")))
    .def("SetCharge",&ost::mol::mm::Topology::SetCharge,(arg("particle_idx"),arg("charge")))
    .def("SetMasses",&WrapSetMasses,(arg("masses")))
    .def("SetMass",&ost::mol::mm::Topology::SetMass,(arg("particle_idx"),arg("mass")))
    .def("GetSigmas",&WrapGetSigmas)
    .def("GetEpsilons",&WrapGetEpsilons)
    .def("GetGBSARadii",&WrapGetGBSARadii)
    .def("GetOBCScalings",&WrapGetOBCScalings)
    .def("GetCharges",&WrapGetCharges)
    .def("GetMasses",&WrapGetMasses)
    .def("GetCharge",&ost::mol::mm::Topology::GetCharge,(arg("particle_idx")))
    .def("GetMass",&ost::mol::mm::Topology::GetMass,(arg("particle_idx")))
    .def("GetOBCScaling",&ost::mol::mm::Topology::GetOBCScaling,(arg("particle_idx")))
    .def("GetGBSARadius",&ost::mol::mm::Topology::GetGBSARadius,(arg("particle_idx")))

    //getter and setter functions for nonbonded fudge parameters
    .def("SetFudgeLJ",&ost::mol::mm::Topology::SetFudgeLJ,(arg("fudge_param")))
    .def("SetFudgeQQ",&ost::mol::mm::Topology::SetFudgeQQ,(arg("fudge_param")))
    .def("GetFudgeLJ",&ost::mol::mm::Topology::GetFudgeLJ)
    .def("GetFudgeQQ",&ost::mol::mm::Topology::GetFudgeQQ)
 
    //getter functions for interaction parameters
    .def("GetHarmonicBondParameters",&WrapGetHarmonicBondParam,(arg("interaction_idx")))
    .def("GetHarmonicAngleParameters",&WrapGetHarmonicAngleParam,(arg("interaction_idx")))
    .def("GetUreyBradleyAngleParameters",&WrapGetUreyBradleyAngleParam,(arg("interaction_idx")))
    .def("GetPeriodicDihedralParameters",&WrapGetPeriodicDihedralParam,(arg("interaction_idx")))
    .def("GetPeriodicImproperParameters",&WrapGetPeriodicImproperParam,(arg("interaction_idx")))
    .def("GetHarmonicImproperParameters",&WrapGetHarmonicImproperParam,(arg("interaction_idx")))
    .def("GetCMapParameters",&WrapGetCMapParam,(arg("interaction_idx")))
    .def("GetLJPairParameters",&WrapGetLJPairParam,(arg("interaction_idx")))
    .def("GetDistanceConstraintParameters",&WrapGetDistanceConstraintParam,(arg("interaction_idx")))
    .def("GetHarmonicPositionRestraintParameters",&WrapGetHarmonicPositionRestraintParam,(arg("interaction_idx")))
    .def("GetHarmonicDistanceRestraintParameters",&WrapGetHarmonicDistanceRestraintParam,(arg("interaction_idx")))
    .def("GetFGMDHBondDonorParameters",&WrapGetFGMDHBondDonorParam,(arg("interaction_idx")))
    .def("GetFGMDHBondAcceptorParameters",&WrapGetFGMDHBondAcceptorParam,(arg("interaction_idx")))

    //setter functions for interaction parameters
    .def("SetHarmonicBondParameters",&ost::mol::mm::Topology::SetHarmonicBondParameters,(arg("interaction_idx"),arg("bond_length"),arg("force_constant")))
    .def("SetHarmonicAngleParameters",&ost::mol::mm::Topology::SetHarmonicAngleParameters,(arg("interaction_idx"),arg("angle"),arg("force_constant")))
    .def("SetUreyBradleyAngleParameters",&ost::mol::mm::Topology::SetUreyBradleyAngleParameters,(arg("interaction_idx"),arg("angle"),arg("angle_force_constant"),arg("bond_length"),arg("bond_force_constant")))
    .def("SetPeriodicDihedralParameters",&ost::mol::mm::Topology::SetPeriodicDihedralParameters,(arg("interaction_idx"),arg("multiplicity"),arg("phase"),arg("force_constant")))
    .def("SetPeriodicImproperParameters",&ost::mol::mm::Topology::SetPeriodicImproperParameters,(arg("interaction_idx"),arg("multiplicity"),arg("phase"),arg("force_constant")))
    .def("SetHarmonicImproperParameters",&ost::mol::mm::Topology::SetHarmonicImproperParameters,(arg("interaction_idx"),arg("angle"),arg("force_constant")))
    .def("SetCMapParameters",&WrapSetCMapParameters,(arg("interaction_idx"),arg("map")))
    .def("SetLJPairParameters",&ost::mol::mm::Topology::SetLJPairParameters,(arg("interaction_idx"),arg("sigma"),arg("epsilon")))
    .def("SetDistanceConstraintParameters",&ost::mol::mm::Topology::SetDistanceConstraintParameters,(arg("interaction_idx"),arg("distance")))
    .def("SetHarmonicPositionRestraintParameters",&ost::mol::mm::Topology::SetHarmonicPositionRestraintParameters,(arg("interaction_idx"),arg("ref_position"),arg("k"),arg("x_scale")=1.0,arg("y_scale")=1.0,arg("z_scale")=1.0))
    .def("SetHarmonicDistanceRestraintParameters",&ost::mol::mm::Topology::SetHarmonicDistanceRestraintParameters,(arg("interaction_idx"),arg("length"),arg("force_constant")))
    .def("SetFGMDHBondDonorParameters",&ost::mol::mm::Topology::SetFGMDHBondDonorParameters,(arg("interaction_idx"),arg("length"),arg("k_length"),arg("alpha"),arg("k_alpha"),arg("beta"),arg("k_beta")))

    //functions to find interactions certain atoms are involved in
    .def("GetHarmonicBondIndices",&GetHarmonicBondIndices,(arg("idx_one"),arg("idx_two")))
    .def("GetHarmonicBondIndices",&GetHarmonicBondIndicesSingleIndex,(arg("idx")))
    .def("GetHarmonicAngleIndices",&GetHarmonicAngleIndices,(arg("idx_one"),arg("idx_two"),arg("idx_three")))
    .def("GetHarmonicAngleIndices",&GetHarmonicAngleIndicesSingleIndex,(arg("idx")))
    .def("GetUreyBradleyAngleIndices",&GetUreyBradleyAngleIndices,(arg("idx_one"),arg("idx_two"),arg("idx_three")))
    .def("GetUreyBradleyAngleIndices",&GetUreyBradleyAngleIndicesSingleIndex,(arg("idx")))
    .def("GetPeriodicDihedralIndices",&GetPeriodicDihedralIndices,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("idx_four")))
    .def("GetPeriodicDihedralIndices",&GetPeriodicDihedralIndicesSingleIndex,(arg("idx")))
    .def("GetPeriodicImproperIndices",&GetPeriodicImproperIndices,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("idx_four")))
    .def("GetPeriodicImproperIndices",&GetPeriodicImproperIndicesSingleIndex,(arg("idx")))
    .def("GetHarmonicImproperIndices",&GetHarmonicImproperIndices,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("idx_four")))
    .def("GetHarmonicImproperIndices",&GetHarmonicImproperIndicesSingleIndex,(arg("idx")))
    .def("GetDistanceConstraintIndex",&ost::mol::mm::Topology::GetDistanceConstraintIndex,(arg("idx_one"),arg("idx_two")))
    .def("GetDistanceConstraintIndices",&GetDistanceConstraintIndicesSingleIndex,(arg("idx")))    
    .def("GetCMapIndices",&GetCMapIndices,(arg("idx_one"),arg("idx_two"),arg("idx_three"),arg("idx_four"),arg("idx_five")))
    .def("GetCMapIndices",&GetCMapIndicesSingleIndex,(arg("idx")))
    .def("GetLJPairIndex",&ost::mol::mm::Topology::GetLJPairIndex,(arg("idx_one"),arg("idx_two")))
    .def("GetLJPairIndices",&GetLJPairIndicesSingleIndex,(arg("idx")))
    .def("GetHarmonicDistanceRestraintIndices",&GetHarmonicDistanceRestraintIndices,(arg("idx_one"),arg("idx_two")))
    .def("GetHarmonicDistanceRestraintIndices",&GetHarmonicDistanceRestraintIndicesSingleIndex,(arg("idx")))
    .def("GetHarmonicPositionRestraintIndices",&GetHarmonicPositionRestraintIndicesSingleIndex,(arg("idx")))
    .def("GetFGMDHBondDonorIndices",&GetFGMDHBondDonorIndices,(arg("idx_one"),arg("idx_two")))
    .def("GetFGMDHBondDonorIndices",&GetFGMDHBondDonorIndicesSingleIndex,(arg("idx")))
    .def("GetFGMDHBondAcceptorIndices",&GetFGMDHBondAcceptorIndices,(arg("idx_one"),arg("idx_two")))
    .def("GetFGMDHBondAcceptorIndices",&GetFGMDHBondAcceptorIndicesSingleIndex,(arg("idx")))


    //functions to get amount of data in topology
    .def("GetNumParticles",&ost::mol::mm::Topology::GetNumParticles)
    .def("GetNumHarmonicBonds",&ost::mol::mm::Topology::GetNumHarmonicBonds)
    .def("GetNumHarmonicAngles",&ost::mol::mm::Topology::GetNumHarmonicAngles)
    .def("GetNumUreyBradleyAngles",&ost::mol::mm::Topology::GetNumUreyBradleyAngles)
    .def("GetNumPeriodicDihedrals",&ost::mol::mm::Topology::GetNumPeriodicDihedrals)
    .def("GetNumPeriodicImpropers",&ost::mol::mm::Topology::GetNumPeriodicImpropers)
    .def("GetNumHarmonicImpropers",&ost::mol::mm::Topology::GetNumHarmonicImpropers)
    .def("GetNumCMaps",&ost::mol::mm::Topology::GetNumCMaps)
    .def("GetNumLJPairs",&ost::mol::mm::Topology::GetNumLJPairs)
    .def("GetNumDistanceConstraints",&ost::mol::mm::Topology::GetNumDistanceConstraints)
    .def("GetNumExclusions",&ost::mol::mm::Topology::GetNumExclusions)
    .def("GetNumHarmonicPositionRestraints",&ost::mol::mm::Topology::GetNumHarmonicPositionRestraints)
    .def("GetNumHarmonicDistanceRestraints",&ost::mol::mm::Topology::GetNumHarmonicDistanceRestraints)  
    .def("GetNumFGMDHBondDonors",&ost::mol::mm::Topology::GetNumFGMDHBondDonors)  
    .def("GetNumFGMDHBondAcceptors",&ost::mol::mm::Topology::GetNumFGMDHBondAcceptors)  

    .def("Merge",&MergeTop,(arg("other_topology")))
    .def("Merge",&MergeTopEnt,(arg("ent"),arg("other_topology"),arg("other_ent")))                                                           
  ;


  boost::python::register_ptr_to_python<ost::mol::mm::TopologyPtr>();
}
