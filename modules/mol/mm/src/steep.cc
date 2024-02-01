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

#include <ost/mol/mm/steep.hh>
#include <OpenMM.h>

namespace ost{ namespace mol{ namespace mm{

bool Steep::Minimize(OpenMM::Context& context, double tolerance, uint max_iterations) {
  const OpenMM::System& system = context.getSystem();
  if(system.getNumConstraints() > 0){
    throw ost::Error("Cannot handle constraints in steepest descent, use lbgfs instead!");
  }
  if( tolerance <= 0.0 ){
    throw ost::Error("tolerance value must be greater 0.0!"); //also important to avoid zero division error when calculating normalization factor for forces
  }

  bool converged = false;

  OpenMM::State state = context.getState(OpenMM::State::Energy | OpenMM::State::Forces | OpenMM::State::Positions);
  std::vector<OpenMM::Vec3> positions = state.getPositions();
  std::vector<OpenMM::Vec3> new_positions = positions;
  std::vector<OpenMM::Vec3> normalized_forces = state.getForces();
  double energy = state.getPotentialEnergy();
  double h = 0.01; 
  double new_energy,max_length,factor;
  bool normalize_forces = true;
  std::vector<uint> non_zero_mass_particles;
  for(int i = 0; i < system.getNumParticles(); ++i){
    if(system.getParticleMass(i) != 0){
      non_zero_mass_particles.push_back(i);
    }
  }
  for(uint i = 0; i < max_iterations; ++i){
    if(normalize_forces){
      const std::vector<OpenMM::Vec3>& forces = state.getForces();
      max_length = 0.0;
      for(std::vector<uint>::iterator j  = non_zero_mass_particles.begin(); 
          j != non_zero_mass_particles.end(); ++j){
        normalized_forces[*j] = forces[*j];
        max_length = std::max(max_length,forces[*j].dot(forces[*j]));
      }
      max_length = sqrt(max_length);
      if(max_length <= tolerance){
        converged = true;
        break;
      }
      factor = 1.0/max_length;
      for(std::vector<uint>::iterator j = non_zero_mass_particles.begin(); 
          j != non_zero_mass_particles.end(); ++j){
        normalized_forces[*j] *= factor;
      }
    }
    for(std::vector<uint>::iterator j = non_zero_mass_particles.begin(); 
        j != non_zero_mass_particles.end(); ++j){
      new_positions[*j] = positions[*j] + normalized_forces[*j] * h;
    }
    context.setPositions(new_positions);
    //context.computeVirtualSites();
    state = context.getState(OpenMM::State::Forces | OpenMM::State::Energy);
    new_energy = state.getPotentialEnergy();
    if(new_energy < energy){
      std::swap(positions,new_positions);
      energy = new_energy;
      h *= 1.2;
      normalize_forces = true;
    }
    else {
      h *= 0.5;
      normalize_forces = false;
    }  
  }
  context.setPositions(positions);
  return converged;
}
}}}//ns
