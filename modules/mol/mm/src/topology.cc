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

#include <ost/mol/mm/topology.hh>


namespace ost{ namespace mol{ namespace mm{



Topology::Topology(const std::vector<Real>& masses){
  num_particles_ = masses.size();
  atom_masses_ = masses;
  fudge_lj_ = 1.0;
  fudge_qq_ = 1.0;
}

TopologyPtr Topology::Load(const String& filename){

  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open topology. File '"
       << filename << "' does not exist";
    throw ost::io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);
  TopologyPtr top_p(new Topology);
  ds >> *top_p;
  return top_p;
}

void Topology::Save(const String& filename){
  std::ofstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSink ds(stream);
  ds << *this;  
}

uint Topology::AddHarmonicBond(uint index_one,
                               uint index_two, 
                               Real bond_length,
                               Real force_constant){

  if(index_one >= num_particles_ || index_two >= num_particles_){
    throw ost::Error("Index of bond particle exceeds number of particles present in topology!");
  }
  Index<2> index(index_one,index_two);
  std::vector<Real> parameters;
  parameters.push_back(bond_length);
  parameters.push_back(force_constant);
  harmonic_bonds_.push_back(std::make_pair(index,parameters));
  return harmonic_bonds_.size()-1;
}

uint Topology::AddHarmonicAngle(uint index_one,
                                uint index_two,
                                uint index_three,
                                Real angle,
                                Real force_constant){

  if(uint(index_one) >= num_particles_ || uint(index_two) >= num_particles_ || 
     uint(index_three) >= num_particles_){
    throw ost::Error("Index of angle particle exceeds number of particles present in topology!");
  }
  Index<3> index(index_one,index_two,index_three); 
  std::vector<Real> parameters;
  parameters.push_back(angle);
  parameters.push_back(force_constant);
  harmonic_angles_.push_back(std::make_pair(index,parameters));
  return harmonic_angles_.size()-1;
}

uint Topology::AddUreyBradleyAngle(uint index_one,
                                   uint index_two,
                                   uint index_three,
                                   Real angle,
                                   Real angle_force_constant,
                                   Real bond_length,
                                   Real bond_force_constant){

  if(index_one >= num_particles_ || index_two >= num_particles_ || 
     index_three >= num_particles_){
    throw ost::Error("Index of angle particle exceeds number of particles present in topology!");
  }
  Index<3> index(index_one,index_two,index_three); 
  std::vector<Real> parameters;
  parameters.push_back(angle);
  parameters.push_back(angle_force_constant);
  parameters.push_back(bond_length);
  parameters.push_back(bond_force_constant);
  urey_bradley_angles_.push_back(std::make_pair(index,parameters));
  return urey_bradley_angles_.size()-1;
}

uint Topology::AddPeriodicDihedral(uint index_one,
                                   uint index_two,
                                   uint index_three,
                                   uint index_four,
                                   int multiplicity,
                                   Real phase,
                                   Real force_constant){

  if(index_one >= num_particles_ || index_two >= num_particles_||
     index_three >= num_particles_ || index_four >= num_particles_){
    throw ost::Error("Index of dihedral particle exceeds number of particles present in topology!");
  }
  Index<4> index(index_one,index_two,index_three,index_four);
  std::vector<Real> parameters;
  parameters.push_back(multiplicity);
  parameters.push_back(phase);
  parameters.push_back(force_constant);
  periodic_dihedrals_.push_back(std::make_pair(index,parameters));
  return periodic_dihedrals_.size()-1;
}

uint Topology::AddPeriodicImproper(uint index_one,
                                   uint index_two,
                                   uint index_three,
                                   uint index_four,
                                   int multiplicity,
                                   Real phase,
                                   Real force_constant){

  if(index_one >= num_particles_ || index_two >= num_particles_||
     index_three >= num_particles_ || index_four >= num_particles_){
    throw ost::Error("Index of improper particle exceeds number of particles present in topology!");
  }
  Index<4> index(index_one,index_two,index_three,index_four);
  std::vector<Real> parameters;
  parameters.push_back(multiplicity);
  parameters.push_back(phase);
  parameters.push_back(force_constant);
  periodic_impropers_.push_back(std::make_pair(index,parameters));
  return periodic_impropers_.size()-1;
}

uint Topology::AddHarmonicImproper(uint index_one,
                                   uint index_two,
                                   uint index_three,
                                   uint index_four,
                                   Real angle,
                                   Real force_constant){

  if(index_one >= num_particles_ || index_two >= num_particles_||
     index_three >= num_particles_ || index_four >= num_particles_){
    throw ost::Error("Index of improper particle exceeds number of particles present in topology!");
  }
  Index<4> index(index_one,index_two,index_three,index_four);
  std::vector<Real> parameters;
  parameters.push_back(angle);
  parameters.push_back(force_constant);
  harmonic_impropers_.push_back(std::make_pair(index,parameters));
  return harmonic_impropers_.size()-1;
}

uint Topology::AddCMap(uint index_one,
                       uint index_two,
                       uint index_three,
                       uint index_four,
                       uint index_five,
                       int dimension,
                       std::vector<Real> values){

  if(index_one >= num_particles_ || index_two >= num_particles_ ||
     index_three >= num_particles_ || index_four >= num_particles_ ||
     index_five >= num_particles_){
    throw ost::Error("Index of cmap particle exceeds number of particles present in topology!");
  }
  if(uint(dimension * dimension) != values.size()){
    throw ost::Error("Dimension of CMap is inconsistent with its values!");
  }
  Index<5> index(index_one,index_two,index_three,index_four,index_five);
  cmaps_.push_back(std::make_pair(index,values));
  return cmaps_.size()-1;
}

uint Topology::AddLJPair(uint index_one, 
                         uint index_two,
                         Real sigma,
                         Real epsilon){

  if(index_one >= num_particles_ || index_two >= num_particles_){
    throw ost::Error("Index of lj pair particle exceeds number of particles present in topology!");
  }
  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  if(added_lj_pairs_.find(index) != added_lj_pairs_.end() || added_lj_pairs_.find(reverse_index) != added_lj_pairs_.end()){
    throw ost::Error("This particular lj 1,4-interaction is already parametrized!");
  }
  added_lj_pairs_.insert(index);
  std::vector<Real> parameters;
  parameters.push_back(sigma);
  parameters.push_back(epsilon);
  lj_pairs_.push_back(std::make_pair(index,parameters));
  return lj_pairs_.size()-1;
}

uint Topology::AddExclusion(uint index_one,
                            uint index_two){
  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  if(index_one >= num_particles_ || index_two >= num_particles_){
    throw ost::Error("Index of exclusion particle exceeds number of particles present in topology!");
  }
  if(added_exclusions_.find(index) != added_exclusions_.end() || added_exclusions_.find(reverse_index) != added_exclusions_.end()){
    throw ost::Error("This particular exclusion has already been set!");
  }
  added_exclusions_.insert(index);
  exclusions_.push_back(index);
  return exclusions_.size()-1;
}

uint Topology::AddDistanceConstraint(uint index_one, 
                                     uint index_two,
                                     Real distance){

  if(index_one >= num_particles_ || index_two >= num_particles_){
    throw ost::Error("Index of distance constraint atom exceeds number of particles present in topology!");
  }
  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  if(added_distance_constraints_.find(index) != added_distance_constraints_.end() || added_distance_constraints_.find(reverse_index) != added_distance_constraints_.end()){
    throw ost::Error("This particular distance constraint is already parametrized!");
  }
  added_distance_constraints_.insert(index);
  std::vector<Real> parameters;
  parameters.push_back(distance);
  distance_constraints_.push_back(std::make_pair(index,parameters));
  return distance_constraints_.size()-1;
}

void Topology::AddPositionConstraint(uint index){
  if(index >= num_particles_){
    throw ost::Error("Index of position constraint exceeds number of particles present in topology!");
  }
  std::vector<uint>::iterator it = std::find(position_constraints_.begin(),position_constraints_.end(),index);
  if(it != position_constraints_.end()) return; //already present
  position_constraints_.push_back(index);
}

uint Topology::AddHarmonicPositionRestraint(uint ind, const geom::Vec3& ref_position, Real k, 
                                            Real x_scale, Real y_scale, Real z_scale){
  if(ind >= num_particles_){
    throw ost::Error("Index of harmonic distance constraint exceeds number of particles present in topology!");
  }
  Index<1> index(ind);
  std::vector<Real> parameters;
  parameters.push_back(ref_position[0]);
  parameters.push_back(ref_position[1]);
  parameters.push_back(ref_position[2]);
  parameters.push_back(k);
  parameters.push_back(x_scale);
  parameters.push_back(y_scale);
  parameters.push_back(z_scale);

  harmonic_position_restraints_.push_back(std::make_pair(index,parameters));
  return harmonic_position_restraints_.size()-1;
}

uint Topology::AddHarmonicDistanceRestraint(uint index_one, uint index_two,
                                            Real length, Real force_constant){
  if(index_one >= num_particles_ || index_two >= num_particles_){
    throw ost::Error("Index of distance restraint atom exceeds number of particles present in topology!");
  }
  Index<2> index(index_one,index_two);
  std::vector<Real> parameters;
  parameters.push_back(length);
  parameters.push_back(force_constant);
  harmonic_distance_restraints_.push_back(std::make_pair(index,parameters));
  return harmonic_distance_restraints_.size()-1;
}

uint Topology::AddFGMDHBondDonor(uint index_one, uint index_two,
                                 Real length, Real k_length, Real alpha, 
                                 Real k_alpha, Real beta, Real k_beta){
  if(index_one >= num_particles_ || index_two >= num_particles_){
    throw ost::Error("Index of fgmd hbond donor atom exceeds number of particles present in topology!");
  }
  Index<2> index(index_one,index_two);
  std::vector<Real> parameters;
  parameters.push_back(length);
  parameters.push_back(k_length);
  parameters.push_back(alpha);
  parameters.push_back(k_alpha);
  parameters.push_back(beta);
  parameters.push_back(k_beta);
  fgmd_hbond_donors_.push_back(std::make_pair(index,parameters));
  return fgmd_hbond_donors_.size()-1;
}

uint Topology::AddFGMDHBondAcceptor(uint index_one, uint index_two){
  if(index_one >= num_particles_ || index_two >= num_particles_){
    throw ost::Error("Index of fgmd hbond acceptor atom exceeds number of particles present in topology!");
  }
  Index<2> index(index_one,index_two);
  fgmd_hbond_acceptors_.push_back(index);
  return fgmd_hbond_acceptors_.size()-1;
}

void Topology::SetSigmas(const std::vector<Real>& sigmas){
  if(sigmas.size() != num_particles_){
    throw ost::Error("Expect the same number of sigma parameters than particles!");
  }
  sigmas_ = sigmas;
}

void Topology::SetSigma(uint index, Real sigma){
  if(sigmas_.empty()){
    throw ost::Error("Modification of sigma parameter requires initial assignment of sigma parameters in globo!");
  }
  if(index >= num_particles_){
    throw ost::Error("Given particle index exceeds number of available particles!");
  }
  sigmas_[index] = sigma;
}

void Topology::SetEpsilons(const std::vector<Real>& epsilons){
  if(epsilons.size() != num_particles_){
    throw ost::Error("Expect the same number of epsilon parameters than particles!");
  }
  epsilons_ = epsilons;
}

void Topology::SetEpsilon(uint index, Real epsilon){
  if(epsilons_.empty()){
    throw ost::Error("Modification of epsilon parameter requires initial assignment of epsilon parameters in globo!");
  }
  if(index >= num_particles_){
    throw ost::Error("Given particle index exceeds number of available particles!");
  }
  epsilons_[index] = epsilon;
}

void Topology::SetGBSARadii(const std::vector<Real>& gbsa_radii){
  if(gbsa_radii.size() != num_particles_){
    throw ost::Error("Expect the same number of gbsa parameters than particles!");
  }
  gbsa_radii_ = gbsa_radii;
}

void Topology::SetGBSARadius(uint index, Real radius){
  if(gbsa_radii_.empty()){
    throw ost::Error("Modification of gbsa radius requires initial assignment of gbsa radii in globo!");
  }
  if(index >= num_particles_){
    throw ost::Error("Given particle index exceeds number of available particles!");
  }
  gbsa_radii_[index] = radius;
}

void Topology::SetOBCScalings(const std::vector<Real>& obc_scaling){
  if(obc_scaling.size() != num_particles_){
    throw ost::Error("Expect the same number of gbsa parameters than particles!");
  }
  obc_scaling_ = obc_scaling;
}

void Topology::SetOBCScaling(uint index, Real scaling){
  if(obc_scaling_.empty()){
    throw ost::Error("Modification of obc scaling parameter requires initial assignment of obc scaling parameters in globo!");
  }
  if(index >= num_particles_){
    throw ost::Error("Given particle index exceeds number of available particles!");
  }
  obc_scaling_[index] = scaling;
}

void Topology::SetCharges(const std::vector<Real>& charges){
  if(charges.size() != num_particles_){
    throw ost::Error("Expect the same number of charges than particles!");
  }
  charges_ = charges;
}

void Topology::SetCharge(uint index, Real charge){
  if(charges_.empty()){
    throw ost::Error("Modification of charge parameter requires initial assignment of charges in globo!");
  }
  if(index >= num_particles_){
    throw ost::Error("Given particle index exceeds number of available particles!");
  }
  charges_[index] = charge;  
}

void Topology::SetMasses(const std::vector<Real>& masses){
  if(masses.size() != num_particles_){
    throw ost::Error("Expect masses vector to have the same number of elements than particles in the topologies entity!");
  }
  atom_masses_ = masses;
}

void Topology::SetMass(uint index, Real mass){
  if(atom_masses_.empty()){
    throw ost::Error("Modification of mass parameter requires initial assignment of masses in globo!");
  }
  if(index >= num_particles_){
    throw ost::Error("Given particle index exceeds number of available particles!");
  }
  atom_masses_[index] = mass;

}


void Topology::GetHarmonicBondParameters(uint index, uint& index_one, uint& index_two,
                                         Real& bond_length, Real& force_constant) const{
  if(index >= harmonic_bonds_.size()){
    throw ost::Error("Provided index exceeds number of harmonic bonds present in topology!");
  }
  bond_length = harmonic_bonds_[index].second[0];
  force_constant = harmonic_bonds_[index].second[1];
  index_one = harmonic_bonds_[index].first[0];
  index_two = harmonic_bonds_[index].first[1];
}

void Topology::GetHarmonicAngleParameters(uint index, uint& index_one, uint& index_two, uint& index_three,
                                          Real& angle, Real& force_constant) const{
  if(index >= harmonic_angles_.size()){
    throw ost::Error("Provided index exceeds number of harmonic angles present in topology!");
  }
  angle = harmonic_angles_[index].second[0];
  force_constant = harmonic_angles_[index].second[1];
  index_one = harmonic_angles_[index].first[0];
  index_two = harmonic_angles_[index].first[1];
  index_three = harmonic_angles_[index].first[2];
}

void Topology::GetUreyBradleyAngleParameters(uint index, uint& index_one, uint& index_two, uint& index_three,
                                             Real& angle, Real& angle_force_constant, Real& bond_length, Real& bond_force_constant) const{
  if(index >= urey_bradley_angles_.size()){
    throw ost::Error("Provided index exceeds number of urey_bradley angles present in topology!");
  }
  angle = urey_bradley_angles_[index].second[0];
  angle_force_constant = urey_bradley_angles_[index].second[1];
  bond_length = urey_bradley_angles_[index].second[2];
  bond_force_constant = urey_bradley_angles_[index].second[3];
  index_one = urey_bradley_angles_[index].first[0];
  index_two = urey_bradley_angles_[index].first[1];
  index_three = urey_bradley_angles_[index].first[2];
}

void Topology::GetPeriodicDihedralParameters(uint index, uint& index_one, uint& index_two, uint& index_three, uint& index_four,
                                             int& multiplicity, Real& phase, Real& force_constant) const{
  if(index >= periodic_dihedrals_.size()){
    throw ost::Error("Provided index exceeds number of periodic dihedrals present in topology!");
  }
  multiplicity = int(periodic_dihedrals_[index].second[0]);
  phase = periodic_dihedrals_[index].second[1];
  force_constant = periodic_dihedrals_[index].second[2];
  index_one = periodic_dihedrals_[index].first[0];
  index_two = periodic_dihedrals_[index].first[1];
  index_three = periodic_dihedrals_[index].first[2];
  index_four = periodic_dihedrals_[index].first[3];
}

void Topology::GetPeriodicImproperParameters(uint index, uint& index_one, uint& index_two, uint& index_three, uint& index_four,
                                             int& multiplicity, Real& phase, Real& force_constant) const{
  if(index >= periodic_impropers_.size()){
    throw ost::Error("Provided index exceeds number of periodic impropers present in topology!");
  }
  multiplicity = int(periodic_impropers_[index].second[0]);
  phase = periodic_impropers_[index].second[1];
  force_constant = periodic_impropers_[index].second[2];
  index_one = periodic_impropers_[index].first[0];
  index_two = periodic_impropers_[index].first[1];
  index_three = periodic_impropers_[index].first[2];
  index_four = periodic_impropers_[index].first[3];
}

void Topology::GetHarmonicImproperParameters(uint index, uint& index_one, uint& index_two, uint& index_three, uint& index_four,
                                             Real& angle, Real& force_constant) const{
  if(index >= harmonic_impropers_.size()){
    throw ost::Error("Provided index exceeds number of harmonic impropers present in topology!");
  }
  angle = harmonic_impropers_[index].second[0];
  force_constant = harmonic_impropers_[index].second[1];
  index_one = harmonic_impropers_[index].first[0];
  index_two = harmonic_impropers_[index].first[1];
  index_three = harmonic_impropers_[index].first[2];
  index_four = harmonic_impropers_[index].first[3];
}

void Topology::GetCMapParameters(uint index, uint& index_one, uint& index_two, uint& index_three, uint& index_four, uint& index_five,
                                 int& dimension, std::vector<Real>& map) const{
  if(index >= cmaps_.size()){
    throw ost::Error("Provided index exceeds number of cmaps present in topology!");
  }
  dimension = int(sqrt(cmaps_[index].second.size()));
  map = cmaps_[index].second;
  index_one = cmaps_[index].first[0];
  index_two = cmaps_[index].first[1];
  index_three = cmaps_[index].first[2];
  index_four = cmaps_[index].first[3];
  index_five = cmaps_[index].first[4];
  
}

void Topology::GetLJPairParameters(uint index, uint& index_one, uint& index_two,
                                   Real& sigma, Real& epsilon) const{
  if(index >= lj_pairs_.size()){
    throw ost::Error("Provided index exceeds number of lj pairs present in topology!");
  }
  sigma = lj_pairs_[index].second[0];
  epsilon = lj_pairs_[index].second[1];
  index_one = lj_pairs_[index].first[0];
  index_two = lj_pairs_[index].first[1];
}

void Topology::GetDistanceConstraintParameters(uint index, uint& index_one, uint& index_two,
                                               Real& distance) const{
  if(index >= distance_constraints_.size()){
    throw ost::Error("Provided index exceeds number of distance_constraints present in topology!");
  }
  distance = distance_constraints_[index].second[0];
  index_one = distance_constraints_[index].first[0];
  index_two = distance_constraints_[index].first[1];
}

void Topology::GetHarmonicPositionRestraintParameters(uint index, uint& atom_index, geom::Vec3& ref_position,
                                               Real& k, Real& x_scale, Real& y_scale, Real& z_scale) const{
  if(index >= harmonic_position_restraints_.size()){
    throw ost::Error("Provided index exceeds numnber of harmonic position constraints present in topology!");
  }
  atom_index = harmonic_position_restraints_[index].first[0];
  ref_position[0] = harmonic_position_restraints_[index].second[0];
  ref_position[1] = harmonic_position_restraints_[index].second[1];
  ref_position[2] = harmonic_position_restraints_[index].second[2];
  k = harmonic_position_restraints_[index].second[3];
  x_scale = harmonic_position_restraints_[index].second[4];
  y_scale = harmonic_position_restraints_[index].second[5];
  z_scale = harmonic_position_restraints_[index].second[6];
}

void Topology::GetHarmonicDistanceRestraintParameters(uint index, uint& atom_one, uint& atom_two, Real& length,
                                            Real& force_constant) const{
  if(index >= harmonic_distance_restraints_.size()){
    throw ost::Error("Provided index exceeds number of harmonic distance restraints present in topology!");
  }
  atom_one = harmonic_distance_restraints_[index].first[0];
  atom_two = harmonic_distance_restraints_[index].first[1];
  length = harmonic_distance_restraints_[index].second[0];
  force_constant = harmonic_distance_restraints_[index].second[1];
}

void Topology::GetFGMDHBondDonorParameters(uint index, uint& atom_one, uint& atom_two, Real& length,
                                           Real& k_length, Real& alpha, Real& k_alpha,
                                           Real& beta, Real& k_beta) const{
  if(index >= fgmd_hbond_donors_.size()){
    throw ost::Error("Provided index exceeds number of fgmd hbond donors present in topology!");
  }
  atom_one = harmonic_distance_restraints_[index].first[0];
  atom_two = harmonic_distance_restraints_[index].first[1];
  length = fgmd_hbond_donors_[index].second[0];
  k_length = fgmd_hbond_donors_[index].second[1];  
  alpha = fgmd_hbond_donors_[index].second[2];
  k_alpha = fgmd_hbond_donors_[index].second[3];
  beta = fgmd_hbond_donors_[index].second[4];
  k_beta = fgmd_hbond_donors_[index].second[5];
}

void Topology::GetFGMDHBondAcceptorParameters(uint index, uint& atom_one, uint& atom_two) const{
  if(index >= fgmd_hbond_acceptors_.size()){
    throw ost::Error("Provided index exceeds number of fgmd hbond acceptors present in topology!");
  }
  atom_one = fgmd_hbond_acceptors_[index][0];
  atom_two = fgmd_hbond_acceptors_[index][1];
}

void Topology::SetHarmonicBondParameters(uint index, const Real bond_length, const Real force_constant){
  if(index >= harmonic_bonds_.size()){
    throw ost::Error("Provided index exceeds number of harmonic bonds present in topology!");
  }
  harmonic_bonds_[index].second[0] = bond_length;
  harmonic_bonds_[index].second[1] = force_constant;
}

void Topology::SetHarmonicAngleParameters(uint index, const Real angle, const Real force_constant){
  if(index >= harmonic_angles_.size()){
    throw ost::Error("Provided index exceeds number of harmonic angles present in topology!");
  }
  harmonic_angles_[index].second[0] = angle;
  harmonic_angles_[index].second[1] = force_constant;
}

void Topology::SetUreyBradleyAngleParameters(uint index, const Real angle, const Real angle_force_constant, const Real bond_length, const Real bond_force_constant){
  if(index >= urey_bradley_angles_.size()){
    throw ost::Error("Provided index exceeds number of urey_bradley angles present in topology!");
  }
  urey_bradley_angles_[index].second[0] = angle;
  urey_bradley_angles_[index].second[1] = angle_force_constant;
  urey_bradley_angles_[index].second[2] = bond_length;
  urey_bradley_angles_[index].second[3] = bond_force_constant;
}

void Topology::SetPeriodicDihedralParameters(uint index, const int multiplicity, const Real phase, const Real force_constant){
  if(index >= periodic_dihedrals_.size()){
    throw ost::Error("Provided index exceeds number of periodic dihedrals present in topology!");
  }
  periodic_dihedrals_[index].second[0] = multiplicity;
  periodic_dihedrals_[index].second[1] = phase;
  periodic_dihedrals_[index].second[2] = force_constant;
}

void Topology::SetPeriodicImproperParameters(uint index, const int multiplicity, const Real phase, const Real force_constant){
  if(index >= periodic_impropers_.size()){
    throw ost::Error("Provided index exceeds number of periodic impropers present in topology!");
  }
  periodic_impropers_[index].second[0] = multiplicity;
  periodic_impropers_[index].second[1] = phase;
  periodic_impropers_[index].second[2] = force_constant;
}

void Topology::SetHarmonicImproperParameters(uint index, const Real angle, const Real force_constant){
  if(index >= harmonic_impropers_.size()){
    throw ost::Error("Provided index exceeds number of harmonic impropers present in topology!");
  }
  harmonic_impropers_[index].second[0] = angle;
  harmonic_impropers_[index].second[1] = force_constant;
}

void Topology::SetCMapParameters(uint index, const int dimension, const std::vector<Real>& map){
  if(index >= cmaps_.size()){
    throw ost::Error("Provided index exceeds number of cmaps present in topology!");
  }
  if(uint(dimension * dimension) != map.size()){
    throw ost::Error("Provided dimension and map are inconsistent!");
  }
  cmaps_[index].second = map;
}

void Topology::SetLJPairParameters(uint index, const Real sigma, const Real epsilon){
  if(index >= lj_pairs_.size()){
    throw ost::Error("Provided index exceeds number of lj pairs present in topology!");
  }
  lj_pairs_[index].second[0] = sigma;
  lj_pairs_[index].second[1] = epsilon;
}

void Topology::SetDistanceConstraintParameters(uint index, const Real distance){
  if(index >= distance_constraints_.size()){
    throw ost::Error("Provided index exceeds number of distance_constraints present in topology!");
  }
  distance_constraints_[index].second[0] = distance;
}

void Topology::SetHarmonicPositionRestraintParameters(uint index, const geom::Vec3& ref_position, Real k, 
                                                      Real x_scale, Real y_scale, Real z_scale){
  if(index >= harmonic_position_restraints_.size()){
    throw ost::Error("Provided index exceeds number of harmonic position restraints present in topology!");
  }
  harmonic_position_restraints_[index].second[0] = ref_position[0];
  harmonic_position_restraints_[index].second[1] = ref_position[1];
  harmonic_position_restraints_[index].second[2] = ref_position[2];
  harmonic_position_restraints_[index].second[3] = k;
  harmonic_position_restraints_[index].second[4] = x_scale;
  harmonic_position_restraints_[index].second[5] = y_scale;
  harmonic_position_restraints_[index].second[6] = z_scale;
}

void Topology::SetHarmonicDistanceRestraintParameters(uint index, Real length, 
                                                      Real force_constant){
  if(index >= harmonic_distance_restraints_.size()){
    throw ost::Error("Provided index exceeds number of harmonic distance restraints present in topology!");
  }
  harmonic_distance_restraints_[index].second[0] = length;
  harmonic_distance_restraints_[index].second[1] = force_constant; 
}

void Topology::SetFGMDHBondDonorParameters(uint index, Real length, Real k_length, 
                                           Real alpha, Real k_alpha, Real beta, Real k_beta){
  if(index >= fgmd_hbond_donors_.size()){
    throw ost::Error("Provided index exceeds number of fgmd donors present in topology!");
  }
  fgmd_hbond_donors_[index].second[0] = length;
  fgmd_hbond_donors_[index].second[1] = k_length;
  fgmd_hbond_donors_[index].second[2] = alpha; 
  fgmd_hbond_donors_[index].second[3] = k_alpha; 
  fgmd_hbond_donors_[index].second[4] = beta; 
  fgmd_hbond_donors_[index].second[5] = k_beta; 
}

Real Topology::GetMass(uint index) const{
  if(index >= num_particles_) throw ost::Error("Provided index exceeds number of particles present in topology!");
  return atom_masses_[index];
}

Real Topology::GetSigma(uint index) const{
  if(index >= num_particles_) throw ost::Error("Provided index exceeds number of particles present in topology!");
  if(sigmas_.empty()) throw ost::Error("No sigmas set!");
  return sigmas_[index];
}

Real Topology::GetEpsilon(uint index) const{
  if(index >= num_particles_) throw ost::Error("Provided index exceeds number of particles present in topology!");
  if(epsilons_.empty()) throw ost::Error("No epsilons set!");
  return epsilons_[index];
}

Real Topology::GetGBSARadius(uint index) const{
  if(index >= num_particles_) throw ost::Error("Provided index exceeds number of particles present in topology!");
  if(gbsa_radii_.empty()) throw ost::Error("No gbsa radii set!");
  return gbsa_radii_[index];
}

Real Topology::GetOBCScaling(uint index) const{
  if(index >= num_particles_) throw ost::Error("Provided index exceeds number of particles present in topology!");
  if(obc_scaling_.empty()) throw ost::Error("No obc scalings set!");
  return obc_scaling_[index];
}


Real Topology::GetCharge(uint index) const{
  if(index >= num_particles_) throw ost::Error("Provided index exceeds number of particles present in topology!");
  if(charges_.empty()) throw ost::Error("No charges set!");
  return charges_[index];
}

std::vector<uint> Topology::GetHarmonicBondIndices(uint index_one,
                                                    uint index_two) const{
  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_bonds_.size(); ++i){
    if(index == harmonic_bonds_[i].first || reverse_index == harmonic_bonds_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetHarmonicAngleIndices(uint index_one,
                                                     uint index_two,
                                                     uint index_three) const{
  Index<3> index(index_one,index_two,index_three);
  Index<3> reverse_index(index_three,index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_angles_.size(); ++i){
    if(index == harmonic_angles_[i].first || reverse_index == harmonic_angles_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetUreyBradleyAngleIndices(uint index_one,
                                                        uint index_two,
                                                        uint index_three) const{
  Index<3> index(index_one,index_two,index_three);
  Index<3> reverse_index(index_three,index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < urey_bradley_angles_.size(); ++i){
    if(index == urey_bradley_angles_[i].first || reverse_index == urey_bradley_angles_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetPeriodicDihedralIndices(uint index_one,
                                                        uint index_two,
                                                        uint index_three,
                                                        uint index_four) const{
  Index<4> index(index_one,index_two,index_three,index_four);
  Index<4> reverse_index(index_four,index_three,index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < periodic_dihedrals_.size(); ++i){
    if(index == periodic_dihedrals_[i].first || reverse_index == periodic_dihedrals_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetPeriodicImproperIndices(uint index_one,
                                                        uint index_two,
                                                        uint index_three,
                                                        uint index_four) const{
  Index<4> index(index_one,index_two,index_three,index_four);
  Index<4> reverse_index(index_four,index_three,index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < periodic_impropers_.size(); ++i){
    if(index == periodic_impropers_[i].first || reverse_index == periodic_impropers_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetHarmonicImproperIndices(uint index_one,
                                                        uint index_two,
                                                        uint index_three,
                                                        uint index_four) const{
  Index<4> index(index_one,index_two,index_three,index_four);
  Index<4> reverse_index(index_four,index_three,index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_impropers_.size(); ++i){
    if(index == harmonic_impropers_[i].first || reverse_index == harmonic_impropers_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetCMapIndices(uint index_one,
                                            uint index_two,
                                            uint index_three,
                                            uint index_four,
                                            uint index_five) const{
  //note, that in case of cmaps the order of the atoms matters
  Index<5> index(index_one,index_two,index_three,index_four,index_five);
  std::vector<uint> return_indices;
  for(uint i = 0; i < cmaps_.size(); ++i){
    if(index == cmaps_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

int Topology::GetLJPairIndex(uint index_one,
                                             uint index_two) const{
  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  for(uint i = 0; i < lj_pairs_.size(); ++i){
    if(index == lj_pairs_[i].first || reverse_index == lj_pairs_[i].first){
      return i;
    }
  }
  return -1;
}

int Topology::GetDistanceConstraintIndex(uint index_one,
                                         uint index_two) const{
  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  for(uint i = 0; i < distance_constraints_.size(); ++i){
    if(index == distance_constraints_[i].first || reverse_index == distance_constraints_[i].first){
      return i;
    }
  }
  return -1;
}

std::vector<uint> Topology::GetHarmonicPositionRestraintIndices(uint atom_index) const{
  Index<1> index(atom_index);
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_position_restraints_.size(); ++i){
    if(index == harmonic_position_restraints_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetHarmonicDistanceRestraintIndices(uint index_one, uint index_two) const{

  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_distance_restraints_.size(); ++i){
    if(index == harmonic_distance_restraints_[i].first || reverse_index == harmonic_distance_restraints_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetFGMDHBondDonorIndices(uint index_one, uint index_two) const{

  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < fgmd_hbond_donors_.size(); ++i){
    if(index == fgmd_hbond_donors_[i].first || reverse_index == fgmd_hbond_donors_[i].first){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetFGMDHBondAcceptorIndices(uint index_one, uint index_two) const{

  Index<2> index(index_one,index_two);
  Index<2> reverse_index(index_two,index_one);
  std::vector<uint> return_indices;
  for(uint i = 0; i < fgmd_hbond_acceptors_.size(); ++i){
    if(index == fgmd_hbond_acceptors_[i] || reverse_index == fgmd_hbond_acceptors_[i]){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetHarmonicBondIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_bonds_.size(); ++i){
    if(harmonic_bonds_[i].first[0] == atom_index || harmonic_bonds_[i].first[1] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetHarmonicAngleIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_angles_.size(); ++i){
    if(harmonic_angles_[i].first[0] == atom_index || harmonic_angles_[i].first[1] == atom_index ||
       harmonic_angles_[i].first[2] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetUreyBradleyAngleIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < urey_bradley_angles_.size(); ++i){
    if(urey_bradley_angles_[i].first[0] == atom_index || urey_bradley_angles_[i].first[1] == atom_index ||
       urey_bradley_angles_[i].first[2] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetPeriodicDihedralIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < periodic_dihedrals_.size(); ++i){
    if(periodic_dihedrals_[i].first[0] == atom_index || periodic_dihedrals_[i].first[1] == atom_index ||
       periodic_dihedrals_[i].first[2] == atom_index || periodic_dihedrals_[i].first[3] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetPeriodicImproperIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < periodic_impropers_.size(); ++i){
    if(periodic_impropers_[i].first[0] == atom_index || periodic_impropers_[i].first[1] == atom_index ||
       periodic_impropers_[i].first[2] == atom_index || periodic_impropers_[i].first[3] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetHarmonicImproperIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_impropers_.size(); ++i){
    if(harmonic_impropers_[i].first[0] == atom_index || harmonic_impropers_[i].first[1] == atom_index ||
       harmonic_impropers_[i].first[2] == atom_index || harmonic_impropers_[i].first[3] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetCMapIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < cmaps_.size(); ++i){
    if(cmaps_[i].first[0] == atom_index || cmaps_[i].first[1] == atom_index ||
       cmaps_[i].first[2] == atom_index || cmaps_[i].first[3] == atom_index ||
       cmaps_[i].first[4] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetLJPairIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < lj_pairs_.size(); ++i){
    if(lj_pairs_[i].first[0] == atom_index || lj_pairs_[i].first[1] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetDistanceConstraintIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < distance_constraints_.size(); ++i){
    if(distance_constraints_[i].first[0] == atom_index || distance_constraints_[i].first[1] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetHarmonicDistanceRestraintIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < harmonic_distance_restraints_.size(); ++i){
    if(harmonic_distance_restraints_[i].first[0] == atom_index || harmonic_distance_restraints_[i].first[1] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetFGMDHBondDonorIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < fgmd_hbond_donors_.size(); ++i){
    if(fgmd_hbond_donors_[i].first[0] == atom_index || fgmd_hbond_donors_[i].first[1] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

std::vector<uint> Topology::GetFGMDHBondAcceptorIndices(uint atom_index) const{
  std::vector<uint> return_indices;
  for(uint i = 0; i < fgmd_hbond_acceptors_.size(); ++i){
    if(fgmd_hbond_acceptors_[i][0] == atom_index || fgmd_hbond_acceptors_[i][1] == atom_index){
      return_indices.push_back(i);
    }
  }
  return return_indices;
}

void Topology::Merge(ost::mol::EntityHandle& ent, TopologyPtr other, 
                     const ost::mol::EntityHandle& other_ent){

  this->CheckEntToAdd(ent,other,other_ent);
  this->CheckTopToAdd(other);
  this->MergeEnt(ent, other_ent);
  this->MergeTop(other);
}

void Topology::Merge(TopologyPtr other){

  this->CheckTopToAdd(other);
  this->MergeTop(other);
}

void Topology::MergeTop(TopologyPtr other){

  uint old_num_particles = num_particles_;
  num_particles_ = num_particles_ + other->GetNumParticles();

  //let's map over masses
  atom_masses_.resize(old_num_particles + other->GetNumParticles());
  memcpy(&atom_masses_[old_num_particles],&other->atom_masses_[0],other->GetNumParticles() * sizeof(Real));

  //let's map over the other per atom parameters if they're defined in the destination topology
  if(!charges_.empty()){
    charges_.resize(old_num_particles + other->GetNumParticles());
    memcpy(&charges_[old_num_particles],&other->charges_[0],other->GetNumParticles() * sizeof(Real));
  }

  if(!sigmas_.empty()){
    sigmas_.resize(old_num_particles + other->GetNumParticles());
    memcpy(&sigmas_[old_num_particles],&other->sigmas_[0],other->GetNumParticles() * sizeof(Real));
  }

  if(!epsilons_.empty()){
    epsilons_.resize(old_num_particles + other->GetNumParticles());
    memcpy(&epsilons_[old_num_particles],&other->epsilons_[0],other->GetNumParticles() * sizeof(Real));
  }

  if(!gbsa_radii_.empty()){
    gbsa_radii_.resize(old_num_particles + other->GetNumParticles());
    memcpy(&gbsa_radii_[old_num_particles],&other->gbsa_radii_[0],other->GetNumParticles() * sizeof(Real));
  }

  if(!obc_scaling_.empty()){
    obc_scaling_.resize(old_num_particles + other->GetNumParticles());
    memcpy(&obc_scaling_[old_num_particles],&other->obc_scaling_[0],other->GetNumParticles() * sizeof(Real));
  }

  //finally, all the interactions get mapped over

  const std::vector<std::pair<Index<2>, std::vector<Real> > >& harmonic_bonds = other->GetHarmonicBonds();
  if(!harmonic_bonds.empty()){
    for(std::vector<std::pair<Index<2>, std::vector<Real> > >::const_iterator i = harmonic_bonds.begin();
        i != harmonic_bonds.end(); ++i){
      this->AddHarmonicBond(old_num_particles + i->first[0],
                            old_num_particles + i->first[1],
                            i->second[0],i->second[1]);
    }
  }

  const std::vector<std::pair<Index<3>, std::vector<Real> > >& harmonic_angles = other->GetHarmonicAngles();
  if(!harmonic_angles.empty()){
    for(std::vector<std::pair<Index<3>, std::vector<Real> > >::const_iterator i = harmonic_angles.begin();
        i != harmonic_angles.end(); ++i){
      this->AddHarmonicAngle(old_num_particles + i->first[0],
                             old_num_particles + i->first[1],
                             old_num_particles + i->first[2],
                             i->second[0],i->second[1]);
    }
  }

  const std::vector<std::pair<Index<3>, std::vector<Real> > >& urey_bradley_angles = other->GetUreyBradleyAngles();
  if(!urey_bradley_angles.empty()){
    for(std::vector<std::pair<Index<3>, std::vector<Real> > >::const_iterator i = urey_bradley_angles.begin();
        i != urey_bradley_angles.end(); ++i){
      this->AddUreyBradleyAngle(old_num_particles + i->first[0],
                                old_num_particles + i->first[1],
                                old_num_particles + i->first[2],
                                i->second[0],i->second[1],
                                i->second[2],i->second[3]);
    }
  }

  const std::vector<std::pair<Index<4>, std::vector<Real> > >& periodic_dihedrals = other->GetPeriodicDihedrals();
  if(!periodic_dihedrals.empty()){
    for(std::vector<std::pair<Index<4>, std::vector<Real> > >::const_iterator i = periodic_dihedrals.begin();
        i != periodic_dihedrals.end(); ++i){
      this->AddPeriodicDihedral(old_num_particles + i->first[0],
                                old_num_particles + i->first[1],
                                old_num_particles + i->first[2],
                                old_num_particles + i->first[3],
                                int(i->second[0]),i->second[1],i->second[2]);
    }
  }

  const std::vector<std::pair<Index<4>, std::vector<Real> > >& periodic_impropers = other->GetPeriodicImpropers();
  if(!periodic_impropers.empty()){
    for(std::vector<std::pair<Index<4>, std::vector<Real> > >::const_iterator i = periodic_impropers.begin();
        i != periodic_impropers.end(); ++i){
      this->AddPeriodicImproper(old_num_particles + i->first[0],
                                old_num_particles + i->first[1],
                                old_num_particles + i->first[2],
                                old_num_particles + i->first[3],
                                int(i->second[0]),i->second[1],i->second[2]);
    }
  }

  const std::vector<std::pair<Index<4>, std::vector<Real> > >& harmonic_impropers = other->GetHarmonicImpropers();
  if(!harmonic_impropers.empty()){
    for(std::vector<std::pair<Index<4>, std::vector<Real> > >::const_iterator i = harmonic_impropers.begin();
        i != harmonic_impropers.end(); ++i){
      this->AddHarmonicImproper(old_num_particles + i->first[0],
                                old_num_particles + i->first[1],
                                old_num_particles + i->first[2],
                                old_num_particles + i->first[3],
                                i->second[0],i->second[1]);
    }
  }

  const std::vector<std::pair<Index<5>, std::vector<Real> > >& cmaps = other->GetCMaps();
  if(!cmaps.empty()){
    for(std::vector<std::pair<Index<5>, std::vector<Real> > >::const_iterator i = cmaps.begin();
        i != cmaps.end(); ++i){
      uint dimension = sqrt(i->second.size());
      this->AddCMap(old_num_particles + i->first[0],
                    old_num_particles + i->first[1],
                    old_num_particles + i->first[2],
                    old_num_particles + i->first[3],
                    old_num_particles + i->first[4],
                    dimension,i->second);
    }
  }

  const std::vector<std::pair<Index<2>, std::vector<Real> > >& ljpairs = other->GetLJPairs();
  if(!ljpairs.empty()){
    for(std::vector<std::pair<Index<2>, std::vector<Real> > >::const_iterator i = ljpairs.begin();
        i != ljpairs.end(); ++i){
      this->AddLJPair(old_num_particles + i->first[0],
                      old_num_particles + i->first[1],
                      i->second[0],i->second[1]);
    }
  }


  const std::vector<std::pair<Index<2>, std::vector<Real> > >& distance_constraints = other->GetDistanceConstraints();
  if(!distance_constraints.empty()){
    for(std::vector<std::pair<Index<2>, std::vector<Real> > >::const_iterator i = distance_constraints.begin();
        i != distance_constraints.end(); ++i){
      this->AddDistanceConstraint(old_num_particles + i->first[0],
                                  old_num_particles + i->first[1],
                                  i->second[0]);
    }
  }

  const std::vector<Index<2> >& exclusions = other->GetExclusions();
  if(!exclusions.empty()){
    for(std::vector<Index<2> >::const_iterator i = exclusions.begin();
        i != exclusions.end(); ++i){
      this->AddExclusion(old_num_particles + (*i)[0],
                         old_num_particles + (*i)[1]);
    }
  }

  const std::vector<std::pair<Index<1>,std::vector<Real> > >& harmonic_position_restraints = other->GetHarmonicPositionRestraints();
  if(!harmonic_position_restraints.empty()){
    for(std::vector<std::pair<Index<1>,std::vector<Real> > >::const_iterator i = harmonic_position_restraints.begin();
        i != harmonic_position_restraints.end(); ++i){
      geom::Vec3 ref_pos(i->second[0], i->second[1], i->second[2]);
      this->AddHarmonicPositionRestraint(old_num_particles + i->first[0],
                                         ref_pos,
                                         i->second[3],
                                         i->second[4],
                                         i->second[5],
                                         i->second[6]);
    }
  }

  const std::vector<std::pair<Index<2>,std::vector<Real> > >& harmonic_distance_restraints = other->GetHarmonicDistanceRestraints();
  if(!harmonic_distance_restraints.empty()){
    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::const_iterator i = harmonic_distance_restraints.begin();
        i != harmonic_distance_restraints.end(); ++i){
      this->AddHarmonicDistanceRestraint(old_num_particles + i->first[0],
                                         old_num_particles + i->first[1],
                                         i->second[0],
                                         i->second[1]);
    }
  }

  const std::vector<std::pair<Index<2>,std::vector<Real> > >& donors = other->GetFGMDHBondDonors();
  if(!donors.empty()){
    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::const_iterator i = donors.begin();
        i != donors.end(); ++i){
      this->AddFGMDHBondDonor(old_num_particles + i->first[0],
                              old_num_particles + i->first[1],
                              i->second[0],
                              i->second[1],
                              i->second[2],
                              i->second[3],
                              i->second[4],
                              i->second[5]);
    }
  }

  const std::vector<Index<2> >& acceptors = other->GetFGMDHBondAcceptors();
  if(!acceptors.empty()){
    for(std::vector<Index<2> >::const_iterator i = acceptors.begin();
        i != acceptors.end(); ++i){
      this->AddFGMDHBondAcceptor(old_num_particles + (*i)[0],
                                 old_num_particles + (*i)[1]);
    }
  }

  const std::vector<uint>& position_constraints = other->GetPositionConstraints();
  if(!position_constraints.empty()){
    for(std::vector<uint>::const_iterator i = position_constraints.begin();
        i != position_constraints.end(); ++i){
      this->AddPositionConstraint(old_num_particles + (*i));
    }
  }

  for(std::set<Index<2> >::iterator i = other->added_lj_pairs_.begin();
      i != other->added_lj_pairs_.end(); ++i){
    added_lj_pairs_.insert(Index<2>(old_num_particles + (*i)[0],
                                    old_num_particles + (*i)[1]));
  }

  for(std::set<Index<2> >::iterator i = other->added_distance_constraints_.begin();
      i != other->added_distance_constraints_.end(); ++i){
    added_distance_constraints_.insert(Index<2>(old_num_particles + (*i)[0],
                                                old_num_particles + (*i)[1]));
  }

  for(std::set<Index<2> >::iterator i = other->added_exclusions_.begin();
      i != other->added_exclusions_.end(); ++i){
    added_exclusions_.insert(Index<2>(old_num_particles + (*i)[0],
                                      old_num_particles + (*i)[1]));
  }



}

void Topology::MergeEnt(ost::mol::EntityHandle& ent, const ost::mol::EntityHandle& other_ent){

  //mapper of hashcode from source atom to index in added_atom list
  std::map<long,int> index_mapper;
  ost::mol::AtomHandleList added_atoms;

  ost::mol::ChainHandleList other_chains = other_ent.GetChainList();

  //let's create an editor to copy over all chains, residues and atoms
  //from the other topologies entity
  ost::mol::XCSEditor ed = ent.EditXCS(ost::mol::BUFFERED_EDIT);
  for(ost::mol::ChainHandleList::iterator i = other_chains.begin(),
      e = other_chains.end(); i != e; ++i){
    ost::mol::ResidueHandleList res_list = i->GetResidueList();
    ost::mol::ChainHandle added_chain = ed.InsertChain(i->GetName());

    for(ost::mol::ResidueHandleList::iterator j = res_list.begin();
        j != res_list.end(); ++j){
      ost::mol::AtomHandleList atom_list = j->GetAtomList();
      ost::mol::ResidueHandle added_residue = ed.AppendResidue(added_chain,*j);

      for(ost::mol::AtomHandleList::iterator k = atom_list.begin();
          k != atom_list.end(); ++k){
        index_mapper[k->GetHashCode()] = added_atoms.size();
        added_atoms.push_back(ed.InsertAtom(added_residue,*k));
      }
    }
  }

  //let's rebuild the connectivity
  ost::mol::BondHandleList bond_list = other_ent.GetBondList();
  for(ost::mol::BondHandleList::iterator i = bond_list.begin();
      i != bond_list.end(); ++i){
    ed.Connect(added_atoms[index_mapper[i->GetFirst().GetHashCode()]],
               added_atoms[index_mapper[i->GetSecond().GetHashCode()]]);
  }
  ed.UpdateICS();
}


void Topology::CheckTopToAdd(TopologyPtr other){

  if(other->fudge_lj_ != fudge_lj_ || other->fudge_qq_ != fudge_qq_){
    throw ost::Error("Expect the fudge parameters of source topology to consistent with the fudge parameters from the destination topology!");
  }

  if(!charges_.empty()){
    if(other->charges_.empty()){
      throw ost::Error("Cannot merge topology without charges into a topology with defined charges!");
    }
  }

  if(!sigmas_.empty()){
    if(other->sigmas_.empty()){
      throw ost::Error("Cannot merge topology without lj sigmas into a topology with defined sigmas!");
    }
  }

  if(!epsilons_.empty()){
    if(other->epsilons_.empty()){
      throw ost::Error("Cannot merge topology without lj epsilons into a topology with defined epsilons!");
    }
  }

  if(!gbsa_radii_.empty()){
    if(other->gbsa_radii_.empty()){
      throw ost::Error("Cannot merge topology without gbsa radii into a topology with defined radii!");
    }
  }

  if(!obc_scaling_.empty()){
    if(other->obc_scaling_.empty()){
      throw ost::Error("Cannot merge topology without obc scaling into a topology with defined scaling!");
    }
  }
}

void Topology::CheckEntToAdd(ost::mol::EntityHandle& ent, TopologyPtr other,
                             const ost::mol::EntityHandle& other_ent){

  //check whether the particle numbers match
  if(num_particles_ != static_cast<uint>(ent.GetAtomCount())){
    throw ost::Error("Num Particles is not consistent with num atoms of provided ent!");
  }

  if(other->GetNumParticles() != static_cast<uint>(other_ent.GetAtomCount())){
    throw ost::Error("Num Particles is not consistent with num atoms of provided other_ent!");
  } 

  //check whether there is a chain from the new entity is already 
  //present in the actual entity 
  ost::mol::ChainHandleList other_chains = other_ent.GetChainList();
  for(ost::mol::ChainHandleList::iterator i = other_chains.begin(), 
      e = other_chains.end(); i != e; ++i){
    if(ent.FindChain(i->GetName()).IsValid()){
      std::stringstream ss;
      ss << "Chain with name \"" << i->GetName() << "\" from other topology";
      ss << "is already present in destination topology!";
      throw ost::Error(ss.str());
    }
  }
}

}}}//ns
