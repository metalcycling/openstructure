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

#ifndef OST_MM_TOPOLOGY_HH
#define OST_MM_TOPOLOGY_HH

#include <vector>
#include <map>
#include <set>

#include <boost/shared_ptr.hpp>

#include <ost/message.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/mm/index.hh>
#include <ost/mol/mm/buildingblock.hh>
#include <ost/mol/mm/block_modifiers.hh>
#include <ost/mol/mm/forcefield.hh>
#include <ost/mol/mm/index.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/mol/residue_prop.hh>

#include <time.h>



namespace ost { namespace mol{ namespace mm{

class Topology;
typedef boost::shared_ptr<Topology> TopologyPtr;


class Topology{

public:

  Topology(const std::vector<Real>& masses);

  Topology() { } //should not be accessible from Python to avoid messing around
                 //with empty topology

  static TopologyPtr Load(const String& filename);

  void Save(const String& filename);

  uint AddHarmonicBond(uint index_one,
                       uint index_two, 
                       Real bond_length,
                       Real force_constant);

  uint AddHarmonicAngle(uint index_one,
                        uint index_two,
                        uint index_three,
                        Real angle,
                        Real force_constant);

  uint AddUreyBradleyAngle(uint index_one,
                           uint index_two,
                           uint index_three,
                           Real angle,
                           Real angle_force_constant,
                           Real bond_length,
                           Real bond_force_constant);

  uint AddPeriodicDihedral(uint index_one,
                           uint index_two,
                           uint index_three,
                           uint index_four,
                           int multiplicity,
                           Real phase,
                           Real force_constant);

  uint AddPeriodicImproper(uint index_one,
                           uint index_two,
                           uint index_three,
                           uint index_four,
                           int multiplicity,
                           Real phase,
                           Real force_constant);

  uint AddHarmonicImproper(uint index_one,
                           uint index_two,
                           uint index_three,
                           uint index_four,
                           Real angle,
                           Real force_constant);

  uint AddCMap(uint index_one,
               uint index_two,
               uint index_three,
               uint index_four,
               uint index_five,
               int dimension,
               std::vector<Real> values);

  uint AddLJPair(uint index_one, 
                 uint index_two,
                 Real sigma,
                 Real epsilon);


  uint AddDistanceConstraint(uint index_one, 
                             uint index_two,
                             Real distance);

  uint AddExclusion(uint index_one,
                    uint index_two);

  void AddPositionConstraint(uint index);

  void ResetPositionConstraints() { position_constraints_.clear(); }

  void ResetExclusions() { exclusions_.clear(); }

  uint AddHarmonicPositionRestraint(uint index, const geom::Vec3& ref_position, Real k, 
                                    Real x_scale = 1.0, Real y_scale = 1.0, Real z_scale = 1.0);

  uint AddHarmonicDistanceRestraint(uint index_one, uint index_two, 
                                    Real length, Real force_constant);

  uint AddFGMDHBondDonor(uint index_one, uint index_two,
                         Real length, Real k_length, Real alpha, 
                         Real k_alpha, Real beta, Real k_beta);

  uint AddFGMDHBondAcceptor(uint index_one, uint index_two);

  //Single atom parameters are expected to be set at once...
  void SetSigmas(const std::vector<Real>& sigmas);

  void SetSigma(uint index, Real sigma);

  void SetEpsilons(const std::vector<Real>& epsilons);

  void SetEpsilon(uint index, Real epsilon);

  void SetGBSARadii(const std::vector<Real>& gbsa_radii);

  void SetGBSARadius(uint index, Real radius);

  void SetOBCScalings(const std::vector<Real>& obc_scaling);

  void SetOBCScaling(uint index, Real scaling);

  void SetCharges(const std::vector<Real>& charges);

  void SetCharge(uint index, Real charge);

  void SetMasses(const std::vector<Real>& masses);

  void SetMass(uint index, Real mass);

  void SetFudgeQQ(Real fudge) { fudge_qq_ = fudge; }

  void SetFudgeLJ(Real fudge) { fudge_lj_ = fudge;} 


  void GetHarmonicBondParameters(uint index, uint& index_one, uint& index_two,
                                 Real& bond_length, Real& force_constant) const;

  void GetHarmonicAngleParameters(uint index, uint& index_one, uint& index_two, uint& index_three,
                                  Real& angle, Real& force_constant) const;

  void GetUreyBradleyAngleParameters(uint index, uint& index_one, uint& index_two, uint& index_three,
                                     Real& angle, Real& angle_force_constant, Real& bond_length, Real& bond_force_constant) const;

  void GetPeriodicDihedralParameters(uint index, uint& index_one, uint& index_two, uint& index_three, uint& index_four,
                                     int& multiplicity, Real& phase, Real& force_constant) const;

  void GetPeriodicImproperParameters(uint index, uint& index_one, uint& index_two, uint& index_three, uint& index_four,
                                     int& multiplicity, Real& phase, Real& force_constant) const;

  void GetHarmonicImproperParameters(uint index, uint& index_one, uint& index_two, uint& index_three, uint& index_four,
                                     Real& angle, Real& force_constant) const;

  void GetCMapParameters(uint index, uint& index_one, uint& index_two, uint& index_three, uint& index_four, uint& index_five,
                         int& dimension, std::vector<Real>& map) const;

  void GetLJPairParameters(uint index, uint& index_one, uint& index_two,
                           Real& sigma, Real& epsilon) const;

  void GetDistanceConstraintParameters(uint index, uint& index_one, uint& index_two,
                                       Real& distance) const;

  void GetHarmonicPositionRestraintParameters(uint index, uint& atom_index, geom::Vec3& ref_position,
                                              Real& k, Real& x_scale, Real& y_scale, Real& z_scale) const;

  void GetHarmonicDistanceRestraintParameters(uint index, uint& atom_one, uint& atom_two, Real& length,
                                              Real& force_constant) const;

  void GetFGMDHBondDonorParameters(uint index, uint& index_one, uint& index_two,
                                   Real& length, Real& k_length, Real& alpha, 
                                   Real& k_alpha, Real& beta, Real& k_beta) const;

  void GetFGMDHBondAcceptorParameters(uint index, uint& index_one, uint& index_two) const;

  void SetHarmonicBondParameters(uint index, const Real bond_length, const Real force_constant);

  void SetHarmonicAngleParameters(uint index, const Real angle, const Real force_constant);

  void SetUreyBradleyAngleParameters(uint index, const Real angle, const Real angle_force_constant, 
                                     const Real bond_length, const Real bond_force_constant);

  void SetPeriodicDihedralParameters(uint index, const int multiplicity, 
                                     const Real phase, const Real force_constant);

  void SetPeriodicImproperParameters(uint index, const int multiplicity, const Real phase, 
                                     const Real force_constant);

  void SetHarmonicImproperParameters(uint index, const Real angle, const Real force_constant);

  void SetCMapParameters(uint index, const int dimension, const std::vector<Real>& map);

  void SetLJPairParameters(uint index, const Real sigma, const Real epsilon);

  void SetDistanceConstraintParameters(uint index, const Real distance);

  void SetHarmonicPositionRestraintParameters(uint index, const geom::Vec3& ref_position, Real k, 
                                               Real x_scale = 1.0, Real y_scale = 1.0, Real z_scale = 1.0);

  void SetHarmonicDistanceRestraintParameters(uint index, Real length, Real force_constant);

  void SetFGMDHBondDonorParameters(uint index, Real length, Real k_length, 
                                   Real alpha, Real k_alpha, Real beta, 
                                   Real k_beta);

  const std::vector<std::pair<Index<2>, std::vector<Real> > >& GetHarmonicBonds() const { return harmonic_bonds_; }
 
  const std::vector<std::pair<Index<3>, std::vector<Real> > >& GetHarmonicAngles() const { return harmonic_angles_; }
 
  const std::vector<std::pair<Index<3>, std::vector<Real> > >& GetUreyBradleyAngles() const { return urey_bradley_angles_; }
 
  const std::vector<std::pair<Index<4>, std::vector<Real> > >& GetPeriodicDihedrals() const { return periodic_dihedrals_; }
 
  const std::vector<std::pair<Index<4>, std::vector<Real> > >& GetPeriodicImpropers() const { return periodic_impropers_; }
 
  const std::vector<std::pair<Index<4>, std::vector<Real> > >& GetHarmonicImpropers() const { return harmonic_impropers_; }
 
  const std::vector<std::pair<Index<5>, std::vector<Real> > >& GetCMaps() const { return cmaps_; }
 
  const std::vector<std::pair<Index<2>, std::vector<Real> > >& GetLJPairs() const { return lj_pairs_; }
 
  const std::vector<std::pair<Index<2>, std::vector<Real> > >& GetDistanceConstraints() const { return distance_constraints_; }

  const std::vector<std::pair<Index<1>, std::vector<Real> > >& GetHarmonicPositionRestraints() const { return harmonic_position_restraints_; }
 
  const std::vector<std::pair<Index<2>, std::vector<Real> > >& GetHarmonicDistanceRestraints() const{ return harmonic_distance_restraints_; }

  const std::vector<Index<2> >& GetExclusions() const { return exclusions_; }
 
  const std::vector<uint>& GetPositionConstraints() const { return position_constraints_; }

  const std::vector<std::pair<Index<2>, std::vector<Real> > >& GetFGMDHBondDonors() const { return fgmd_hbond_donors_; }

  const std::vector<Index<2> >& GetFGMDHBondAcceptors() const { return fgmd_hbond_acceptors_; }


  std::vector<Real> GetSigmas() const { return sigmas_; }

  std::vector<Real> GetEpsilons() const { return epsilons_; }
 
  std::vector<Real> GetGBSARadii() const { return gbsa_radii_; }

  std::vector<Real> GetOBCScalings() const { return obc_scaling_; }
 
  std::vector<Real> GetCharges() const { return charges_; }

  Real GetCharge(uint index) const;

  Real GetMass(uint index) const;

  Real GetSigma(uint index) const;

  Real GetEpsilon(uint index) const;

  Real GetGBSARadius(uint index) const;

  Real GetOBCScaling(uint index) const;

  std::vector<Real> GetMasses() const { return atom_masses_; }

  Real GetFudgeQQ() const { return fudge_qq_; }

  Real GetFudgeLJ() const { return fudge_lj_; }

  std::vector<uint> GetHarmonicBondIndices(uint index_one,
                                            uint index_two) const;

  std::vector<uint> GetHarmonicAngleIndices(uint index_one,
                                             uint index_two,
                                             uint index_three) const;

  std::vector<uint> GetUreyBradleyAngleIndices(uint index_one,
                                                uint index_two,
                                                uint index_three) const;

  std::vector<uint> GetPeriodicDihedralIndices(uint index_one,
                                                uint index_two,
                                                uint index_three,
                                                uint index_four) const;

  std::vector<uint> GetPeriodicImproperIndices(uint index_one,
                                                uint index_two,
                                                uint index_three,
                                                uint index_four) const;

  std::vector<uint> GetHarmonicImproperIndices(uint index_one,
                                                uint index_two,
                                                uint index_three,
                                                uint index_four) const;

  std::vector<uint> GetCMapIndices(uint index_one,
                                    uint index_two,
                                    uint index_three,
                                    uint index_four,
                                    uint index_five) const;

  int GetLJPairIndex(uint index_one,
                     uint index_two) const;

  int GetDistanceConstraintIndex(uint index_one,
                                 uint index_two) const;



  std::vector<uint> GetHarmonicDistanceRestraintIndices(uint index_one,
                                                        uint index_two) const;

  std::vector<uint> GetFGMDHBondDonorIndices(uint index_one,
                                             uint index_two) const;

  std::vector<uint> GetFGMDHBondAcceptorIndices(uint index_one,
                                                uint index_two) const;

  std::vector<uint> GetHarmonicBondIndices(uint atom_index) const;

  std::vector<uint> GetHarmonicAngleIndices(uint atom_index) const;

  std::vector<uint> GetUreyBradleyAngleIndices(uint atom_index) const;

  std::vector<uint> GetPeriodicDihedralIndices(uint atom_index) const;

  std::vector<uint> GetPeriodicImproperIndices(uint atom_index) const;

  std::vector<uint> GetHarmonicImproperIndices(uint atom_index) const;

  std::vector<uint> GetCMapIndices(uint atom_index) const;

  std::vector<uint> GetLJPairIndices(uint atom_index) const;

  std::vector<uint> GetDistanceConstraintIndices(uint atom_index) const;

  std::vector<uint> GetHarmonicPositionRestraintIndices(uint atom_index) const;

  std::vector<uint> GetHarmonicDistanceRestraintIndices(uint atom_index) const;

  std::vector<uint> GetFGMDHBondDonorIndices(uint atom_index) const;

  std::vector<uint> GetFGMDHBondAcceptorIndices(uint atom_index) const;

  uint GetNumParticles() { return num_particles_; }

  uint GetNumHarmonicBonds() { return harmonic_bonds_.size(); }

  uint GetNumHarmonicAngles() { return harmonic_angles_.size(); }

  uint GetNumUreyBradleyAngles() { return urey_bradley_angles_.size(); }

  uint GetNumPeriodicDihedrals() { return periodic_dihedrals_.size(); }

  uint GetNumPeriodicImpropers() { return periodic_impropers_.size(); }

  uint GetNumHarmonicImpropers() { return harmonic_impropers_.size(); }

  uint GetNumCMaps() { return cmaps_.size(); }

  uint GetNumLJPairs() { return lj_pairs_.size(); }

  uint GetNumDistanceConstraints() { return distance_constraints_.size(); }

  uint GetNumPositionConstraints() { return position_constraints_.size(); }

  uint GetNumHarmonicPositionRestraints() { return harmonic_position_restraints_.size(); }

  uint GetNumHarmonicDistanceRestraints() { return harmonic_distance_restraints_.size();}

  uint GetNumExclusions() { return exclusions_.size(); }

  uint GetNumFGMDHBondDonors() { return fgmd_hbond_donors_.size(); }

  uint GetNumFGMDHBondAcceptors() { return fgmd_hbond_acceptors_.size(); }

  void Merge(ost::mol::EntityHandle& ent, TopologyPtr other, const ost::mol::EntityHandle& other_ent);

  void Merge(TopologyPtr other);

  template <typename DS>
  void Serialize(DS& ds){


    uint num_items = 0;
    Index<2> actual_index;

    ds & num_particles_;
    ds & fudge_qq_;
    ds & fudge_lj_;

    if(ds.IsSource()){
      ds & num_items;
      atom_masses_ = std::vector<Real>(num_items);
      ds & num_items;
      sigmas_ = std::vector<Real>(num_items);
      ds & num_items;
      epsilons_ = std::vector<Real>(num_items);
      ds & num_items;
      gbsa_radii_ = std::vector<Real>(num_items);
      ds & num_items;
      obc_scaling_ = std::vector<Real>(num_items);
      ds & num_items;
      charges_ = std::vector<Real>(num_items);
      ds & num_items;
      position_constraints_ = std::vector<uint>(num_items);
    }
    else{
      num_items = atom_masses_.size();
      ds & num_items;
      num_items = sigmas_.size();
      ds & num_items;
      num_items = epsilons_.size();
      ds & num_items;
      num_items = gbsa_radii_.size();
      ds & num_items;
      num_items = obc_scaling_.size();
      ds & num_items;
      num_items = charges_.size();
      ds & num_items;
      num_items = position_constraints_.size();
      ds & num_items;
    }

    for(std::vector<Real>::iterator i = atom_masses_.begin();
        i != atom_masses_.end(); ++i){
      ds & *i;
    }

    for(std::vector<Real>::iterator i = sigmas_.begin();
        i != sigmas_.end(); ++i){
      ds & *i;
    }

    for(std::vector<Real>::iterator i = epsilons_.begin();
        i != epsilons_.end(); ++i){
      ds & *i;
    }

    for(std::vector<Real>::iterator i = gbsa_radii_.begin();
        i != gbsa_radii_.end(); ++i){
      ds & *i;
    }

    for(std::vector<Real>::iterator i = obc_scaling_.begin();
        i != obc_scaling_.end(); ++i){
      ds & *i;
    }

    for(std::vector<Real>::iterator i = charges_.begin();
        i != charges_.end(); ++i){
      ds & *i;
    }

    for(std::vector<uint>::iterator i = position_constraints_.begin();
        i != position_constraints_.end(); ++i){
      ds & *i;
    }

    if(ds.IsSource()){
      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        harmonic_bonds_.push_back(std::make_pair(Index<2>(),std::vector<Real>(2)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        harmonic_angles_.push_back(std::make_pair(Index<3>(),std::vector<Real>(2)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        urey_bradley_angles_.push_back(std::make_pair(Index<3>(),std::vector<Real>(4)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        periodic_dihedrals_.push_back(std::make_pair(Index<4>(),std::vector<Real>(3)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        periodic_impropers_.push_back(std::make_pair(Index<4>(),std::vector<Real>(3)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        harmonic_impropers_.push_back(std::make_pair(Index<4>(),std::vector<Real>(2)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        uint num_values = 0;
        ds & num_values;
        cmaps_.push_back(std::make_pair(Index<5>(),std::vector<Real>(num_values)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        lj_pairs_.push_back(std::make_pair(Index<2>(),std::vector<Real>(2)));
      } 

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        distance_constraints_.push_back(std::make_pair(Index<2>(),std::vector<Real>(1)));
      }

      ds & num_items;
      exclusions_ = std::vector<Index<2> >(num_items);

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        harmonic_position_restraints_.push_back(std::make_pair(Index<1>(),std::vector<Real>(7)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        harmonic_distance_restraints_.push_back(std::make_pair(Index<2>(),std::vector<Real>(2)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        fgmd_hbond_donors_.push_back(std::make_pair(Index<2>(),std::vector<Real>(6)));
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        fgmd_hbond_acceptors_.push_back(Index<2>());
      }
    }
    else{
      num_items = harmonic_bonds_.size();
      ds & num_items;
      num_items = harmonic_angles_.size();
      ds & num_items;
      num_items = urey_bradley_angles_.size();
      ds & num_items;
      num_items = periodic_dihedrals_.size();
      ds & num_items;
      num_items = periodic_impropers_.size();
      ds & num_items;
      num_items = harmonic_impropers_.size();
      ds & num_items;
      num_items = cmaps_.size();
      ds & num_items;
      for(uint i = 0; i < cmaps_.size(); ++i){
        num_items = cmaps_[i].second.size();
        ds & num_items;
      }
      num_items = lj_pairs_.size();
      ds & num_items;
      num_items = distance_constraints_.size();
      ds & num_items;
      num_items = exclusions_.size();
      ds & num_items;
      num_items = harmonic_position_restraints_.size();
      ds & num_items;
      num_items = harmonic_distance_restraints_.size();
      ds & num_items;
      num_items = fgmd_hbond_donors_.size();
      ds & num_items;
      num_items = fgmd_hbond_acceptors_.size();
      ds & num_items;
    }

    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::iterator i = harmonic_bonds_.begin();
        i != harmonic_bonds_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
    }

    for(std::vector<std::pair<Index<3>,std::vector<Real> > >::iterator i = harmonic_angles_.begin();
        i != harmonic_angles_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
    }

    for(std::vector<std::pair<Index<3>,std::vector<Real> > >::iterator i = urey_bradley_angles_.begin();
        i != urey_bradley_angles_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
      ds & i->second[2];
      ds & i->second[3];
    }

    for(std::vector<std::pair<Index<4>,std::vector<Real> > >::iterator i = periodic_dihedrals_.begin();
        i != periodic_dihedrals_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
      ds & i->second[2];
    }

    for(std::vector<std::pair<Index<4>,std::vector<Real> > >::iterator i = periodic_impropers_.begin();
        i != periodic_impropers_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
      ds & i->second[2];
    }

    for(std::vector<std::pair<Index<4>,std::vector<Real> > >::iterator i = harmonic_impropers_.begin();
        i != harmonic_impropers_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
    }

    for(std::vector<std::pair<Index<5>,std::vector<Real> > >::iterator i = cmaps_.begin();
        i != cmaps_.end(); ++i){
      ds & i->first;
      for(std::vector<Real>::iterator j = i->second.begin();
          j != i->second.end(); ++j){
        ds & (*j);
      }
    }

    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::iterator i = lj_pairs_.begin();
        i != lj_pairs_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
    }

    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::iterator i = distance_constraints_.begin();
        i != distance_constraints_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
    }

    for(std::vector<Index<2> >::iterator i = exclusions_.begin();
        i != exclusions_.end(); ++i){
      ds & (*i);
    }

    for(std::vector<std::pair<Index<1>,std::vector<Real> > >::iterator i = harmonic_position_restraints_.begin();
        i != harmonic_position_restraints_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
      ds & i->second[2];
      ds & i->second[3];
      ds & i->second[4];
      ds & i->second[5];
      ds & i->second[6];
    }

    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::iterator i = harmonic_distance_restraints_.begin();
        i != harmonic_distance_restraints_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
    }

    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::iterator i = fgmd_hbond_donors_.begin();
        i != fgmd_hbond_donors_.end(); ++i){
      ds & i->first;
      ds & i->second[0];
      ds & i->second[1];
      ds & i->second[2];
      ds & i->second[3];
      ds & i->second[4];
      ds & i->second[5];
    }

    for(std::vector<Index<2> >::iterator i = fgmd_hbond_acceptors_.begin();
        i != fgmd_hbond_acceptors_.end(); ++i){
      ds & (*i);
    }

    if(ds.IsSource()){
      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        ds & actual_index;
        added_lj_pairs_.insert(actual_index);
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        ds & actual_index;
        added_distance_constraints_.insert(actual_index);
      }

      ds & num_items;
      for(uint i = 0; i < num_items; ++i){
        ds & actual_index;
        added_exclusions_.insert(actual_index);
      }
    }
    else{
      num_items = added_lj_pairs_.size();
      ds & num_items;
      for(std::set<Index<2> >::iterator i = added_lj_pairs_.begin(); 
          i != added_lj_pairs_.end(); ++i){
        actual_index = *i;
        ds & actual_index;
      }
      num_items = added_distance_constraints_.size();
      ds & num_items;
      for(std::set<Index<2> >::iterator i = added_distance_constraints_.begin(); 
          i != added_distance_constraints_.end(); ++i){
        actual_index = *i;
        ds & actual_index;
      }
      num_items = added_exclusions_.size();
      ds & num_items;
      for(std::set<Index<2> >::iterator i = added_exclusions_.begin(); 
          i != added_exclusions_.end(); ++i){
        actual_index = *i;
        ds & actual_index;
      }
    }
  }

private:

  void MergeTop(TopologyPtr other);

  void MergeEnt(ost::mol::EntityHandle& ent, const ost::mol::EntityHandle& other_ent);

  void CheckTopToAdd(TopologyPtr other);

  void CheckEntToAdd(ost::mol::EntityHandle& ent, TopologyPtr other, 
                     const ost::mol::EntityHandle& other_ent);

  uint num_particles_;

  //fudge parameters for lj 1,4 pairs
  Real fudge_qq_;
  Real fudge_lj_;

  //single atom parameters
  std::vector<Real> atom_masses_;
  std::vector<Real> sigmas_;
  std::vector<Real> epsilons_;
  std::vector<Real> gbsa_radii_;
  std::vector<Real> obc_scaling_;
  std::vector<Real> charges_;
  std::vector<uint> position_constraints_;

  //interactions with multiple atoms involved
  std::vector<std::pair<Index<2>,std::vector<Real> > > harmonic_bonds_;
  std::vector<std::pair<Index<3>,std::vector<Real> > > harmonic_angles_;
  std::vector<std::pair<Index<3>,std::vector<Real> > > urey_bradley_angles_;
  std::vector<std::pair<Index<4>,std::vector<Real> > > periodic_dihedrals_;
  std::vector<std::pair<Index<4>,std::vector<Real> > > periodic_impropers_;
  std::vector<std::pair<Index<4>,std::vector<Real> > > harmonic_impropers_;
  std::vector<std::pair<Index<5>,std::vector<Real> > > cmaps_;
  std::vector<std::pair<Index<2>,std::vector<Real> > > lj_pairs_;
  std::vector<std::pair<Index<2>,std::vector<Real> > > distance_constraints_;
  std::vector<Index<2> > exclusions_;
  std::vector<std::pair<Index<1>,std::vector<Real> > > harmonic_position_restraints_;
  std::vector<std::pair<Index<2>,std::vector<Real> > > harmonic_distance_restraints_;
  std::vector<std::pair<Index<2>,std::vector<Real> > > fgmd_hbond_donors_;
  std::vector<Index<2> > fgmd_hbond_acceptors_;

  //the atoms of the interactions, that should be unique get tracked in here
  //note, that this is waste of memory, needs better implementation
  std::set<Index<2> > added_lj_pairs_;
  std::set<Index<2> > added_distance_constraints_;
  std::set<Index<2> > added_exclusions_;
};


}}} //ns

#endif
