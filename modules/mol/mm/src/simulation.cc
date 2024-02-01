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

#include <ost/mol/mm/simulation.hh>
#include <OpenMM.h>

namespace ost{ namespace mol{ namespace mm{

Simulation::Simulation(const ost::mol::EntityHandle& handle, 
                       const SettingsPtr settings){

  //note, that ent_ will be "completed" inside this function!
  //(hydrogens and shit)
  
  ost::mol::EntityHandle ent = handle.Copy();
  TopologyPtr top = TopologyCreator::Create(ent,settings);
  this->Init(top, ent, settings);
}

Simulation::Simulation(const TopologyPtr top,
                       const ost::mol::EntityHandle& handle,
                       const SettingsPtr settings){

  if(static_cast<uint>(handle.GetAtomCount()) != top->GetNumParticles()){
    throw ost::Error("Number of atoms in entity must be consistent with number of particles in topology!");
  }
  ost::mol::EntityHandle ent = handle.Copy();
  this->Init(top, ent, settings);
}

void Simulation::Save(const String& filename){

  std::ofstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSink ds(stream);

  ds << *top_;

  uint num_chains;
  uint num_residues;
  uint num_atoms;
  uint num_bonded_atoms;
  Real bfac;
  Real occ;
  bool is_hetatm;
  String chain_name;
  String res_name;
  int resnum_num;
  char resnum_code;
  String atom_name;
  String atom_element;
  uint atom_index;

  num_chains = ent_.GetChainCount();
  ds & num_chains;
  ost::mol::ChainHandleList chain_list = ent_.GetChainList();
  for(ost::mol::ChainHandleList::iterator i = chain_list.begin();
      i != chain_list.end(); ++i){
    chain_name = i->GetName();
    num_residues = i->GetResidueCount();
    ds & chain_name;
    ds & num_residues;
    ost::mol::ResidueHandleList res_list = i->GetResidueList();
    for(ost::mol::ResidueHandleList::iterator j = res_list.begin();
        j != res_list.end(); ++j){
      res_name = j->GetKey();
      resnum_num = j->GetNumber().GetNum();
      resnum_code = j->GetNumber().GetInsCode(); 
      num_atoms = j->GetAtomCount();
      ds & res_name;
      ds & resnum_num;
      ds & resnum_code;
      ds & num_atoms;
      ost::mol::AtomHandleList atom_list = j->GetAtomList();
      for(ost::mol::AtomHandleList::iterator k = atom_list.begin();
          k != atom_list.end(); ++k){
        atom_name = k->GetName();
        atom_element = k->GetElement();
        bfac = k->GetBFactor();
        occ = k->GetOccupancy();
        is_hetatm = k->IsHetAtom();
        ds & atom_name;
        ds & atom_element;
        ds & bfac;
        ds & occ;
        ds & is_hetatm;
      }
    }
  }

  ost::mol::AtomHandleList bonded_atoms;

  std::map<long,int> atom_indices;
  int actual_index = 0;
  for(ost::mol::AtomHandleList::const_iterator i = atom_list_.begin(), 
      e = atom_list_.end(); i != e; ++i){
    atom_indices[i->GetHashCode()] = actual_index;
    ++actual_index;
  }

  for(ost::mol::AtomHandleList::iterator i = atom_list_.begin();
      i != atom_list_.end(); ++i){
    bonded_atoms = i->GetBondPartners();
    num_bonded_atoms = bonded_atoms.size();
    ds & num_bonded_atoms;
    for(ost::mol::AtomHandleList::iterator j = bonded_atoms.begin();
        j != bonded_atoms.end(); ++j){
      atom_index = atom_indices[j->GetHashCode()];
      ds & atom_index;
    }
  }

  context_->createCheckpoint(stream);  
}

SimulationPtr Simulation::Load(const String& filename, SettingsPtr settings){
  
  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open simulation File '"
       << filename << "' does not exist";
    throw ost::io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSource ds(stream);

  SimulationPtr sim_ptr(new Simulation);
  TopologyPtr top(new Topology);
  ds >> *top;

  uint num_chains;
  uint num_residues;
  uint num_atoms;
  uint num_bonded_atoms;
  Real bfac;
  Real occ;
  bool is_hetatm;
  String chain_name;
  String res_name;
  int resnum_num;
  char resnum_code;
  String atom_name;
  String atom_element;
  uint atom_index;

  ost::mol::EntityHandle ent = ost::mol::CreateEntity();
  ost::mol::XCSEditor ed = ent.EditXCS();
  ds & num_chains;
  for(uint i = 0; i < num_chains; ++i){
    ds & chain_name;
    ds & num_residues;
    ost::mol::ChainHandle chain = ed.InsertChain(chain_name);
    for(uint j = 0; j < num_residues; ++j){
      ds & res_name;
      ds & resnum_num;
      ds & resnum_code;
      ds & num_atoms;
      ost::mol::ResNum num(resnum_num,resnum_code);
      ost::mol::ResidueHandle res = ed.AppendResidue(chain,res_name,num);
      for(uint k = 0; k < num_atoms; ++k){
        ds & atom_name;
        ds & atom_element;
        ds & bfac;
        ds & occ;
        ds & is_hetatm;
        ed.InsertAtom(res, atom_name, geom::Vec3(0.0,0.0,0.0),
                      atom_element, occ, bfac, is_hetatm);
      }
    }
  }

  sim_ptr->Init(top, ent, settings);

  for(uint i = 0; i < sim_ptr->atom_list_.size(); ++i){
    ds & num_bonded_atoms;
    for(uint j = 0; j < num_bonded_atoms; ++j){
      ds & atom_index;
      ed.Connect(sim_ptr->atom_list_[i], sim_ptr->atom_list_[atom_index]);
    }
  }

  // also loads the positions that have been set in the context
  // they get mapped over to the attached entity
  sim_ptr->context_->loadCheckpoint(stream);
  sim_ptr->UpdatePositions();

  return sim_ptr;
}


bool Simulation::IsPlatformAvailable(const SettingsPtr settings) {
  // load settings-dependent plugin directories
  EnsurePluginsLoaded(settings->openmm_plugin_directory);
  EnsurePluginsLoaded(settings->custom_plugin_directory);
  // check if OpenMM platform exists by checking all (this is fast..)
  const Platform platform_id = settings->platform;
  for (int i = 0; i < OpenMM::Platform::getNumPlatforms(); ++i) {
    OpenMM::Platform& platform = OpenMM::Platform::getPlatform(i);
    if (   (platform_id == Reference && platform.getName() == "Reference")
        || (platform_id == OpenCL && platform.getName() == "OpenCL")
        || (platform_id == CUDA && platform.getName() == "CUDA")
        || (platform_id == CPU && platform.getName() == "CPU")) {
      return true;
    }
  }
  return false;
}

void Simulation::EnsurePluginsLoaded(const String& plugin_path) {
  // note: this is guaranteed to be constructed on first use only
  static std::set<String> already_loaded;
  if (already_loaded.find(plugin_path) == already_loaded.end()) {
    // not loaded yet: load directory, but only once!
    OpenMM::Platform::loadPluginsFromDirectory(plugin_path);
    already_loaded.insert(plugin_path);
  }
}


void Simulation::Init(const TopologyPtr top,
                      const ost::mol::EntityHandle& ent,
                      const SettingsPtr settings){

  top_ = top;
  ent_ = ent;
  atom_list_ = ent_.GetAtomList();

  if(!settings->integrator){
    //user did not specify an integrator, so let's just use a standard integrator
    settings->integrator = IntegratorPtr(new OpenMM::VerletIntegrator(0.001));
  }
  integrator_ = settings->integrator;

  system_ = SystemCreator::Create(top_,settings,system_force_mapper_); 

  //setting up the context, which combines the system with an integrator
  //to proceed in time, but first we have to load the proper platform

  EnsurePluginsLoaded(settings->openmm_plugin_directory);
  EnsurePluginsLoaded(settings->custom_plugin_directory);
  OpenMM::Platform* platform;
  std::map<String,String> context_properties;

  switch(settings->platform){
    case Reference:{
      platform = &OpenMM::Platform::getPlatformByName("Reference");
      for(PropertyMap::iterator i = settings->reference_properties.begin();
          i != settings->reference_properties.end(); ++i){
        context_properties[i->first] = i->second;
      }
      break;
    }
    case OpenCL:{
      platform = &OpenMM::Platform::getPlatformByName("OpenCL");
      for(PropertyMap::iterator i = settings->opencl_properties.begin();
          i != settings->opencl_properties.end(); ++i){
        context_properties[i->first] = i->second;
      }
      break;
    }
    case CUDA:{
      platform = &OpenMM::Platform::getPlatformByName("CUDA");
      for(PropertyMap::iterator i = settings->cuda_properties.begin();
          i != settings->cuda_properties.end(); ++i){
        context_properties[i->first] = i->second;
      }
      break;
    }
    case CPU:{
      platform = &OpenMM::Platform::getPlatformByName("CPU");
      for(PropertyMap::iterator i = settings->cpu_properties.begin();
          i != settings->cpu_properties.end(); ++i){
        context_properties[i->first] = i->second;
      }
      break;
    }
    default:{
      throw ost::Error("Invalid platform encountered when settings up simualation!");
    }
  }

  context_ = ContextPtr(new OpenMM::Context(*system_,*integrator_,*platform,context_properties));

  std::vector<OpenMM::Vec3> positions;
  geom::Vec3 ost_vec;
  OpenMM::Vec3 open_mm_vec;
  for(ost::mol::AtomHandleList::iterator i = atom_list_.begin();
      i!=atom_list_.end();++i){
    ost_vec = i->GetPos();
    open_mm_vec[0] = ost_vec[0]/10;
    open_mm_vec[1] = ost_vec[1]/10;
    open_mm_vec[2] = ost_vec[2]/10;
    positions.push_back(open_mm_vec);
  }

  context_->setPositions(positions);

  //make sure the context satisfies the distance constraints
  context_->applyConstraints(0.00001);

  if(settings->init_temperature > 0.0){
    context_->setVelocitiesToTemperature(settings->init_temperature);
  }
}

geom::Vec3List Simulation::GetPositions(bool enforce_periodic_box, bool in_angstrom){
  geom::Vec3List return_vec;
  StateExtractor::ExtractPositions(context_,return_vec, enforce_periodic_box, in_angstrom);
  return return_vec;
}

geom::Vec3List Simulation::GetVelocities(){
  geom::Vec3List return_vec;
  StateExtractor::ExtractVelocities(context_, return_vec);
  return return_vec;
}

geom::Vec3List Simulation::GetForces(){
  geom::Vec3List return_vec;
  StateExtractor::ExtractForces(context_, return_vec);
  return return_vec;
}

void Simulation::SetPositions(const geom::Vec3List& positions, bool in_angstrom){
  if(top_->GetNumParticles() != positions.size()){
    throw ost::Error("Number of positions does not correspond to number of particles in topology!");
  }
  std::vector<OpenMM::Vec3> openmm_positions(top_->GetNumParticles());
  OpenMM::Vec3* actual_pos = &openmm_positions[0];
  if(in_angstrom){
    for(geom::Vec3List::const_iterator i = positions.begin();
        i != positions.end(); ++i){
      (*actual_pos)[0] = (*i)[0]*0.1;
      (*actual_pos)[1] = (*i)[1]*0.1;
      (*actual_pos)[2] = (*i)[2]*0.1;
      ++actual_pos;
    }
  }
  else{
    for(geom::Vec3List::const_iterator i = positions.begin();
        i != positions.end(); ++i){
      (*actual_pos)[0] = (*i)[0];
      (*actual_pos)[1] = (*i)[1];
      (*actual_pos)[2] = (*i)[2];
      ++actual_pos;
    }
  }
  context_->setPositions(openmm_positions);
}

void Simulation::SetVelocities(geom::Vec3List& velocities){
  if(top_->GetNumParticles() != velocities.size()){
    throw ost::Error("Number of velocities does not correspond to number of particles in topology!");
  }
  std::vector<OpenMM::Vec3> openmm_velocities(top_->GetNumParticles());
  OpenMM::Vec3* actual_vel = &openmm_velocities[0];
  for(geom::Vec3List::iterator i = velocities.begin();
      i != velocities.end(); ++i){
    (*actual_vel)[0] = (*i)[0];
    (*actual_vel)[1] = (*i)[1];
    (*actual_vel)[2] = (*i)[2];
    ++actual_vel;
  }
  context_->setVelocities(openmm_velocities);
}

void Simulation::UpdatePositions(bool enforce_periodic_box){
  if(top_->GetNumParticles() != static_cast<uint>(ent_.GetAtomCount())){
    throw ost::Error("Num particles in topology and num atoms in entity are not consistent!");
  }
  geom::Vec3List positions;
  StateExtractor::ExtractPositions(context_, positions, enforce_periodic_box, 
                                   true);
  ost::mol::XCSEditor ed = ent_.EditXCS(ost::mol::BUFFERED_EDIT);
  for(uint i = 0; i < atom_list_.size(); ++i) {
    ed.SetAtomPos(atom_list_[i], positions[i]);
  }
}

bool Simulation::ApplySD(Real tolerance, int max_iterations){
  return Steep::Minimize(*context_,tolerance, max_iterations);
}

void Simulation::ApplyLBFGS(Real tolerance, int max_iterations){
  OpenMM::LocalEnergyMinimizer::minimize(*context_, tolerance, max_iterations);
}

void Simulation::Steps(int steps){

  int time_to_next_notification;

  while(steps > 0){
    time_to_next_notification = this->TimeToNextNotification();
    //we can just do the remaining steps and subtract them from
    //time_to_notify_
    if(steps < time_to_next_notification){
      integrator_->step(steps);
      for(std::vector<int>::iterator i = time_to_notify_.begin();
          i != time_to_notify_.end(); ++i){
        (*i) -= steps;
      } 
      steps = 0;
    }else{
      integrator_->step(time_to_next_notification);
      for(std::vector<int>::iterator i = time_to_notify_.begin();
          i != time_to_notify_.end(); ++i){
        (*i) -= time_to_next_notification;
      } 
      for(uint i = 0; i < time_to_notify_.size(); ++i){
        if(time_to_notify_[i] == 0){
          observers_[i]->Notify();
          time_to_notify_[i] = observers_[i]->Rhythm();
        }
      }
      steps -= time_to_next_notification;
    }
  }
}

Real Simulation::GetEnergy(){
  Real energy;
  StateExtractor::ExtractEnergy(context_,energy);
  return energy; 
}

Real Simulation::GetKineticEnergy(){
  Real energy;
  StateExtractor::ExtractKineticEnergy(context_,energy);
  return energy;
}

Real Simulation::GetPotentialEnergy(){
  Real energy;
  StateExtractor::ExtractPotentialEnergy(context_,energy);
  return energy;
}

void Simulation::Register(ObserverPtr o){
  observers_.push_back(o);
  time_to_notify_.push_back(o->Rhythm());
  o->Init(context_,top_,ent_);
}

int Simulation::TimeToNextNotification(){
  if(observers_.empty()){
    return std::numeric_limits<int>::max();
  }
  int minimal_time = std::numeric_limits<int>::max();
  for(std::vector<int>::iterator i = time_to_notify_.begin();
      i != time_to_notify_.end(); ++i){
    if(*i < minimal_time){
      minimal_time = *i;
    }
  }
  return minimal_time;
}

void Simulation::ResetHarmonicBond(uint index, Real bond_length, Real force_constant){
  if(system_force_mapper_.find(HarmonicBond) == system_force_mapper_.end()){
    throw ost::Error("There is no harmonic bond force in the simulation!");
  }
  uint force_index = system_force_mapper_[HarmonicBond];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::HarmonicBondForce* bond_force_ptr = reinterpret_cast<OpenMM::HarmonicBondForce*>(&force); 
  if(index >= uint(bond_force_ptr->getNumBonds())){
    throw ost::Error("Provided bond index exceeds number of bonds!");
  }
  double dummy_one, dummy_two;
  int particle1, particle2;
  //we know the parameters, but not the atom indices in the force object, we have to extract it from
  //the force itself
  bond_force_ptr->getBondParameters(index,particle1,particle2,dummy_one,dummy_two);
  bond_force_ptr->setBondParameters(index,particle1,particle2,bond_length,force_constant);
  bond_force_ptr->updateParametersInContext(*context_);
  top_->SetHarmonicBondParameters(index,bond_length,force_constant);
}

void Simulation::ResetHarmonicAngle(uint index, Real angle, Real force_constant){
  if(system_force_mapper_.find(HarmonicAngle) == system_force_mapper_.end()){
    throw ost::Error("There is no harmonic angle force in the simulation!");
  }
  uint force_index = system_force_mapper_[HarmonicAngle];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::HarmonicAngleForce* angle_force_ptr = reinterpret_cast<OpenMM::HarmonicAngleForce*>(&force); 
  if(index >= uint(angle_force_ptr->getNumAngles())){
    throw ost::Error("Provided angle index exceeds number of angles!");
  }
  double dummy_one, dummy_two;
  int particle1, particle2, particle3;
  //we know the parameters, but not the atom indices in the force object, we have to extract it from
  //the force itself
  angle_force_ptr->getAngleParameters(index,particle1,particle2,particle3,dummy_one,dummy_two);
  angle_force_ptr->setAngleParameters(index,particle1,particle2,particle3,angle,force_constant);
  angle_force_ptr->updateParametersInContext(*context_);
  top_->SetHarmonicAngleParameters(index, angle, force_constant);
}

void Simulation::ResetUreyBradleyAngle(uint index, Real angle, Real angle_force_constant, Real bond_length, Real bond_force_constant){
  if(system_force_mapper_.find(UreyBradleyAngle) == system_force_mapper_.end()){
    throw ost::Error("There is no urey bradley angle force in the simulation!");
  }
  uint force_index = system_force_mapper_[UreyBradleyAngle];
  OpenMM::Force& angle_force = system_->getForce(force_index);
  OpenMM::Force& bond_force = system_->getForce(force_index+1);

  OpenMM::HarmonicAngleForce* angle_force_ptr = reinterpret_cast<OpenMM::HarmonicAngleForce*>(&angle_force); 
  OpenMM::HarmonicBondForce* bond_force_ptr = reinterpret_cast<OpenMM::HarmonicBondForce*>(&bond_force); 
  if(index >= uint(angle_force_ptr->getNumAngles())){
    throw ost::Error("Provided angle index exceeds number of angles!");
  }
  double dummy_one, dummy_two;
  int particle1, particle2, particle3;
  //we know the parameters, but not the atom indices in the force object, we have to extract it from
  //the force itself
  angle_force_ptr->getAngleParameters(index,particle1,particle2,particle3,dummy_one,dummy_two);
  angle_force_ptr->setAngleParameters(index,particle1,particle2,particle3,angle,angle_force_constant);
  bond_force_ptr->setBondParameters(index,particle1,particle3,bond_length,bond_force_constant);
  angle_force_ptr->updateParametersInContext(*context_);
  bond_force_ptr->updateParametersInContext(*context_);
  top_->SetUreyBradleyAngleParameters(index,angle,angle_force_constant,bond_length,bond_force_constant);
}

void Simulation::ResetPeriodicDihedral(uint index, int multiplicity, Real phase, Real force_constant){
  if(system_force_mapper_.find(PeriodicDihedral) == system_force_mapper_.end()){
    throw ost::Error("There is no periodic dihedral force in the simulation!");
  }
  uint force_index = system_force_mapper_[PeriodicDihedral];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::PeriodicTorsionForce* dihedral_force_ptr = reinterpret_cast<OpenMM::PeriodicTorsionForce*>(&force); 
  if(index >= uint(dihedral_force_ptr->getNumTorsions())){
    throw ost::Error("Provided dihedral index exceeds number of dihedrals!");
  }
  int dummy_one;
  double dummy_two, dummy_three;
  int particle1, particle2, particle3, particle4;
  //we know the parameters, but not the atom indices in the force object, we have to extract it from
  //the force itself
  dihedral_force_ptr->getTorsionParameters(index,particle1,particle2,particle3,particle4,dummy_one,dummy_two,dummy_three);
  dihedral_force_ptr->setTorsionParameters(index,particle1,particle2,particle3,particle4,multiplicity,phase,force_constant);
  dihedral_force_ptr->updateParametersInContext(*context_);
  top_->SetPeriodicDihedralParameters(index, multiplicity, phase, force_constant);
}

void Simulation::ResetPeriodicImproper(uint index, int multiplicity, Real phase, Real force_constant){
  if(system_force_mapper_.find(PeriodicImproper) == system_force_mapper_.end()){
    throw ost::Error("There is no periodic improper force in the simulation!");
  }
  uint force_index = system_force_mapper_[PeriodicImproper];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::PeriodicTorsionForce* improper_ptr = reinterpret_cast<OpenMM::PeriodicTorsionForce*>(&force); 
  if(index >= uint(improper_ptr->getNumTorsions())){
    throw ost::Error("Provided improper index exceeds number of impropers!");
  }
  int dummy_one;
  double dummy_two, dummy_three;
  int particle1, particle2, particle3, particle4;
  //we know the parameters, but not the atom indices in the force object, we have to extract it from
  //the force itself
  improper_ptr->getTorsionParameters(index,particle1,particle2,particle3,particle4,dummy_one,dummy_two,dummy_three);
  improper_ptr->setTorsionParameters(index,particle1,particle2,particle3,particle4,multiplicity,phase,force_constant);
  improper_ptr->updateParametersInContext(*context_);
  top_->SetPeriodicImproperParameters(index, multiplicity, phase, force_constant);
}

void Simulation::ResetHarmonicImproper(uint index, Real phase, Real force_constant){
  if(system_force_mapper_.find(HarmonicImproper) == system_force_mapper_.end()){
    throw ost::Error("There is no harmonic improper force in the simulation!");
  }
  uint force_index = system_force_mapper_[HarmonicImproper];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::CustomTorsionForce* improper_ptr = reinterpret_cast<OpenMM::CustomTorsionForce*>(&force); 
  if(index >= uint(improper_ptr->getNumTorsions())){
    throw ost::Error("Provided improper index exceeds number of impropers!");
  }
  std::vector<double> parameters;
  parameters.push_back(0.0);
  parameters.push_back(0.0);
  int particle1, particle2, particle3, particle4;
  //we know the parameters, but not the atom indices in the force object, we have to extract it from
  //the force itself
  improper_ptr->getTorsionParameters(index,particle1,particle2,particle3,particle4,parameters);
  parameters[0] = phase;
  parameters[1] = force_constant;
  improper_ptr->setTorsionParameters(index,particle1,particle2,particle3,particle4,parameters);
  improper_ptr->updateParametersInContext(*context_);
  top_->SetHarmonicImproperParameters(index, phase, force_constant);
}

void Simulation::ResetLJPair(uint index, Real sigma, Real epsilon){
  if(system_force_mapper_.find(LJ) == system_force_mapper_.end()){
    throw ost::Error("There is no nonbonded force in the simulation!");
  }
  uint force_index = system_force_mapper_[LJ];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::NonbondedForce* nonbonded_ptr = reinterpret_cast<OpenMM::NonbondedForce*>(&force); 
  if(index >= uint(nonbonded_ptr->getNumExceptions())){
    throw ost::Error("Provided index exceeds number of defined pairs!");
  }
  double charge_prod, dummy_one, dummy_two;
  int particle1, particle2;
  //we know the parameters, but not the atom indices in the force object, we have to extract it from
  //the force itself
  nonbonded_ptr->getExceptionParameters(index,particle1,particle2,charge_prod,dummy_one,dummy_two);
  nonbonded_ptr->setExceptionParameters(index,particle1,particle2,charge_prod,sigma,epsilon);
  nonbonded_ptr->updateParametersInContext(*context_);
  top_->SetLJPairParameters(index, sigma, epsilon);
}

void Simulation::ResetDistanceConstraint(uint index, Real constraint_length){
  if(index >= uint(system_->getNumConstraints())){
    throw ost::Error("Provided index exceeds number of defined constraints!");
  }
  double dummy;
  int particle1, particle2;
  system_->getConstraintParameters(index,particle1,particle2,dummy);
  system_->setConstraintParameters(index,particle1,particle2,constraint_length);
  this->ReinitializeContext();
  top_->SetDistanceConstraintParameters(index, constraint_length);
}

void Simulation::AddPositionConstraint(uint index){
  if(index >= top_->GetNumParticles()){
    throw ost::Error("Provided index exceeds number of atoms!");
  }
  system_->setParticleMass(index,0.0);
  this->ReinitializeContext();
  top_->AddPositionConstraint(index);
}

void Simulation::AddPositionConstraints(const std::vector<uint>& index){

  for(std::vector<uint>::const_iterator i = index.begin(); 
      i != index.end(); ++i){
    if(*i >= top_->GetNumParticles()){
      throw ost::Error("Provided index exceeds number of atoms!");
    }
    system_->setParticleMass(*i,0.0);
    top_->AddPositionConstraint(*i);
  }
  this->ReinitializeContext();
}

void Simulation::ResetPositionConstraints(){
  std::vector<Real> original_masses = top_->GetMasses();
  for(uint i = 0; i < original_masses.size(); ++i){
    system_->setParticleMass(i,original_masses[i]);
  }
  top_->ResetPositionConstraints();
  this->ReinitializeContext();
}

void Simulation::ResetHarmonicPositionRestraint(uint index, const geom::Vec3& ref_position, Real k,
                                                Real x_scale, Real y_scale, Real z_scale){
  if(index >= uint(top_->GetNumHarmonicPositionRestraints())){
    throw ost::Error("Provided index exceeds number of defined harmonic position restraints!");
  }
  if(system_force_mapper_.find(HarmonicPositionRestraint) == system_force_mapper_.end()){
    throw ost::Error("There is no harmonic position restraint in the simulation!");
  }
  uint force_index = system_force_mapper_[HarmonicPositionRestraint];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::CustomExternalForce* restraint_ptr = reinterpret_cast<OpenMM::CustomExternalForce*>(&force); 
  if(index >= uint(restraint_ptr->getNumParticles())){
    throw ost::Error("Provided index exceeds number of defined harmonic position restraints!");
  }
  int particle;
  std::vector<double> parameters;
  restraint_ptr->getParticleParameters(index,particle,parameters);
  parameters[0] = ref_position[0] * 0.1;
  parameters[1] = ref_position[1] * 0.1;
  parameters[2] = ref_position[2] * 0.1;
  parameters[3] = k;
  parameters[4] = x_scale;
  parameters[5] = y_scale;
  parameters[6] = z_scale;
  restraint_ptr->setParticleParameters(index,particle,parameters);
  restraint_ptr->updateParametersInContext(*context_);
  top_->SetHarmonicPositionRestraintParameters(index,ref_position,k,x_scale,y_scale,z_scale);
}

void Simulation::ResetHarmonicDistanceRestraint(uint index, Real length, Real force_constant){
  if(index >= top_->GetNumHarmonicDistanceRestraints()){
    throw ost::Error("Provided index exceeds number of defined harmonic distance restraints!");
  }
  if(system_force_mapper_.find(HarmonicDistanceRestraint) == system_force_mapper_.end()){
    throw ost::Error("There is no harmonic distance restraint in the simulation!");
  }
  uint force_index = system_force_mapper_[HarmonicDistanceRestraint];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::HarmonicBondForce* restraint_ptr = reinterpret_cast<OpenMM::HarmonicBondForce*>(&force); 
  if(index >= uint(restraint_ptr->getNumBonds())){
    throw ost::Error("Provided index exceeds number of defined harmonic distance restraints!");
  }
  int particle_one,particle_two;
  double dummy_one, dummy_two;
  restraint_ptr->getBondParameters(int(index),particle_one,particle_two,dummy_one,dummy_two);
  restraint_ptr->setBondParameters(int(index),particle_one,particle_two,length,2*force_constant);
  restraint_ptr->updateParametersInContext(*context_);
  top_->SetHarmonicDistanceRestraintParameters(index,length,force_constant);
}

void Simulation::ResetLJ(uint index, Real sigma, Real epsilon){
  if(index >= top_->GetNumParticles()){
    throw ost::Error("Provided index exceeds number of atoms!");
  }
  if(system_force_mapper_.find(LJ) == system_force_mapper_.end()){
    throw ost::Error("There is no nonbonded force in the simulation!");
  }
  uint force_index = system_force_mapper_[LJ];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::NonbondedForce* nonbonded_ptr = reinterpret_cast<OpenMM::NonbondedForce*>(&force); 

  double charge, dummy_one, dummy_two;

  nonbonded_ptr->getParticleParameters(index,charge,dummy_one,dummy_two);
  nonbonded_ptr->setParticleParameters(index,charge,sigma,epsilon);
  nonbonded_ptr->updateParametersInContext(*context_);
  top_->SetSigma(index, sigma);
  top_->SetEpsilon(index, epsilon);
}

void Simulation::ResetGBSA(uint index, Real radius, Real scaling){
  if(index >= top_->GetNumParticles()){
    throw ost::Error("Provided index exceeds number of atoms!");
  }
  if(system_force_mapper_.find(GBSA) == system_force_mapper_.end()){
    throw ost::Error("There is no gbsa force in the simulation!");
  }
  uint force_index = system_force_mapper_[GBSA];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::GBSAOBCForce* gbsa_ptr = reinterpret_cast<OpenMM::GBSAOBCForce*>(&force); 

  double charge, dummy_one, dummy_two;

  gbsa_ptr->getParticleParameters(index,charge,dummy_one,dummy_two);
  gbsa_ptr->setParticleParameters(index,charge,radius,scaling);
  gbsa_ptr->updateParametersInContext(*context_);
  top_->SetGBSARadius(index, radius);
  top_->SetOBCScaling(index, scaling);
}

void Simulation::ResetCharge(uint index, Real charge){
  if(index >= top_->GetNumParticles()){
    throw ost::Error("Provided index exceeds number of atoms!");
  }
  if(system_force_mapper_.find(LJ) == system_force_mapper_.end()){
    throw ost::Error("There is no nonbonded force in the simulation!");
  }
  uint force_index = system_force_mapper_[LJ];
  OpenMM::Force& force = system_->getForce(force_index);
  OpenMM::NonbondedForce* nonbonded_ptr = reinterpret_cast<OpenMM::NonbondedForce*>(&force); 

  double dummy, sigma, epsilon;

  //first update the parameters of the single particle
  nonbonded_ptr->getParticleParameters(index,dummy,sigma,epsilon);
  nonbonded_ptr->setParticleParameters(index,charge,sigma,epsilon);
  
  //the charge product has to be updated in the exceptions (1,4 interactions)
   std::vector<uint> lj_pair_indices = top_->GetLJPairIndices(index);
   Real old_charge = top_->GetCharge(index);
   double charge_prod;
   Real new_charge_prod;
   int particle_one, particle_two;

  for(std::vector<uint>::iterator i = lj_pair_indices.begin();
      i != lj_pair_indices.end(); ++i){
    nonbonded_ptr->getExceptionParameters(*i,particle_one,particle_two,charge_prod,sigma,epsilon);
    new_charge_prod = charge_prod/old_charge*charge; // fudge_qq * q1 * q2
    nonbonded_ptr->setExceptionParameters(*i,particle_one,particle_two,new_charge_prod,sigma,epsilon);
  } 

  nonbonded_ptr->updateParametersInContext(*context_);
  top_->SetCharge(index, charge);

  //the charge also affects gbsa force...
  if(system_force_mapper_.find(GBSA) != system_force_mapper_.end()){
    force_index = system_force_mapper_[GBSA];
    OpenMM::Force& force = system_->getForce(force_index);
    OpenMM::GBSAOBCForce* gbsa_ptr = reinterpret_cast<OpenMM::GBSAOBCForce*>(&force); 

    double radius,scaling;

    gbsa_ptr->getParticleParameters(index,dummy,radius,scaling);
    gbsa_ptr->setParticleParameters(index,charge,radius,scaling);
    gbsa_ptr->updateParametersInContext(*context_); 
  }
}

void Simulation::ResetMass(uint index, Real mass){
  if(index >= top_->GetNumParticles()){
    throw ost::Error("Provided index exceeds number of atoms!");
  }
  system_->setParticleMass(index,mass);
  this->ReinitializeContext();
  top_->SetMass(index,mass);
}

geom::Vec3 Simulation::GetPeriodicBoxExtents(){
  geom::Vec3 vec;
  OpenMM::Vec3 ucell_a, ucell_b, ucell_c;
  OpenMM::State openmm_state = context_->getState(0); //create minimal state
                                                      //to only extract periodic
                                                      //box information
  openmm_state.getPeriodicBoxVectors(ucell_a,ucell_b,ucell_c);
  vec[0] = 10 * sqrt(ucell_a.dot(ucell_a));
  vec[1] = 10 * sqrt(ucell_a.dot(ucell_b));
  vec[2] = 10 * sqrt(ucell_a.dot(ucell_c));
  return vec;
}

void Simulation::SetPeriodicBoxExtents(geom::Vec3& vec){
  OpenMM::Vec3 ucell_a(vec[0]/10.0,0.0,0.0);
  OpenMM::Vec3 ucell_b(0.0,vec[1]/10.0,0.0);
  OpenMM::Vec3 ucell_c(0.0,0.0,vec[2]/10.0);
  context_->setPeriodicBoxVectors(ucell_a,ucell_b,ucell_c);
}

void Simulation::ReinitializeContext() {
  // reinitializing requires to reset all those things!
  // Be aware, state of random number generators etc might not be
  // preserved!

  // openmm uses bitmasks to selectively extract data from the context.
  // The context data might change with different OpenMM versions.
  // Passing -1 automagically sets all bits to true => ALL available
  // data gets extracted. 
  int true_bitmask = -1;
  OpenMM::State state = context_->getState(true_bitmask);
  context_->reinitialize();
  context_->setState(state); 
}

}}}

