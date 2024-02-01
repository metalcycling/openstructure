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

#include <ost/mol/mm/system_creator.hh>
#include <OpenMM.h>

namespace ost{ namespace mol{ namespace mm{

SystemPtr SystemCreator::Create(const TopologyPtr top, 
                                const SettingsPtr settings,
                                std::map<FuncType,uint>& mapper){

  uint mapper_index = 0;

  SystemPtr sys(new OpenMM::System);

  //add particles to the system
  std::vector<Real> masses = top->GetMasses();

  //set masses of position constraints to zero
  const std::vector<uint>& position_constraints = top->GetPositionConstraints();
  if(!position_constraints.empty()){
    for(std::vector<uint>::const_iterator i = position_constraints.begin();
        i != position_constraints.end(); ++i){
      masses[*i] = 0.0;
    }
  }  

  for(std::vector<Real>::iterator i = masses.begin(); i != masses.end(); ++i){
    sys->addParticle(*i);
  }

  //parameters of all added interactions will temporarily be stored in here
  std::vector<Real> parameters;

  const std::vector<std::pair<Index<2>,std::vector<Real> > >& harmonic_bonds = top->GetHarmonicBonds();
  if(!harmonic_bonds.empty()){
    OpenMM::HarmonicBondForce& harmonic_bond_force = *new OpenMM::HarmonicBondForce();
    sys->addForce(&harmonic_bond_force);
    mapper[HarmonicBond] = mapper_index++;
    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::const_iterator i = harmonic_bonds.begin();
        i != harmonic_bonds.end(); ++i){
      harmonic_bond_force.addBond(i->first[0],i->first[1],i->second[0],i->second[1]);
    }
  }

  const std::vector<std::pair<Index<3>, std::vector<Real> > >& harmonic_angles = top->GetHarmonicAngles();
  if(!harmonic_angles.empty()){
    OpenMM::HarmonicAngleForce& harmonic_angle_force = *new OpenMM::HarmonicAngleForce();
    sys->addForce(&harmonic_angle_force);
    mapper[HarmonicAngle] = mapper_index++;
    for(std::vector<std::pair<Index<3>,std::vector<Real> > >::const_iterator i = harmonic_angles.begin();
        i != harmonic_angles.end(); ++i){
      harmonic_angle_force.addAngle(i->first[0],i->first[1],i->first[2],
                                    i->second[0],i->second[1]);
    }
  }
  
  const std::vector<std::pair<Index<3>,std::vector<Real> > >& urey_bradley_angles = top->GetUreyBradleyAngles();
  if(!urey_bradley_angles.empty()){
    //urey bradley is a mixture between a harmonic angle and harmonic bond term,
    //that gets split up
    OpenMM::HarmonicAngleForce& urey_angle_force = *new OpenMM::HarmonicAngleForce();
    OpenMM::HarmonicBondForce& urey_bond_force = *new OpenMM::HarmonicBondForce();
    sys->addForce(&urey_angle_force);
    sys->addForce(&urey_bond_force);
    mapper[UreyBradleyAngle] = mapper_index++;
    ++mapper_index;
    for(std::vector<std::pair<Index<3>,std::vector<Real> > >::const_iterator i = urey_bradley_angles.begin();
        i != urey_bradley_angles.end(); ++i){
      urey_angle_force.addAngle(i->first[0],i->first[1],i->first[2],
                                i->second[0],i->second[1]);
      urey_bond_force.addBond(i->first[0],i->first[2],i->second[2],i->second[3]);
    }
  }

  const std::vector<std::pair<Index<4>,std::vector<Real> > >& periodic_dihedrals = top->GetPeriodicDihedrals();
  if(!periodic_dihedrals.empty()){
    OpenMM::PeriodicTorsionForce& periodic_torsion_force = *new OpenMM::PeriodicTorsionForce();
    sys->addForce(&periodic_torsion_force);
    mapper[PeriodicDihedral] = mapper_index++;
    for(std::vector<std::pair<Index<4>,std::vector<Real> > >::const_iterator i = periodic_dihedrals.begin();
        i != periodic_dihedrals.end(); ++i){
      periodic_torsion_force.addTorsion(i->first[0],i->first[1],i->first[2],i->first[3],
                                        i->second[0],i->second[1],i->second[2]);
    }
  }

  const std::vector<std::pair<Index<4>,std::vector<Real> > >& periodic_impropers = top->GetPeriodicImpropers();
  if(!periodic_impropers.empty()){
    OpenMM::PeriodicTorsionForce& periodic_improper_force = *new OpenMM::PeriodicTorsionForce();       
    sys->addForce(&periodic_improper_force);
    mapper[PeriodicImproper] = mapper_index++;
    for(std::vector<std::pair<Index<4>,std::vector<Real> > >::const_iterator i = periodic_impropers.begin();
        i != periodic_impropers.end(); ++i){
      periodic_improper_force.addTorsion(i->first[0],i->first[1],i->first[2],i->first[3],
                                        i->second[0],i->second[1],i->second[2]);
    }
  }

  const std::vector<std::pair<Index<4>,std::vector<Real> > >& harmonic_impropers = top->GetHarmonicImpropers();
  if(!harmonic_impropers.empty()){
    std::vector<double> custom_parameters;
    custom_parameters.push_back(0.0);
    custom_parameters.push_back(0.0);
    OpenMM::CustomTorsionForce& harmonic_improper_force = *new OpenMM::CustomTorsionForce("0.5*k*(theta-b)^2");
    harmonic_improper_force.addPerTorsionParameter("b");
    harmonic_improper_force.addPerTorsionParameter("k"); 
    sys->addForce(&harmonic_improper_force);
    mapper[HarmonicImproper] = mapper_index++;
    for(std::vector<std::pair<Index<4>,std::vector<Real> > >::const_iterator i = harmonic_impropers.begin();
        i != harmonic_impropers.end(); ++i){
      custom_parameters[0] = i->second[0];
      custom_parameters[1] = i->second[1];
      harmonic_improper_force.addTorsion(i->first[0],i->first[1],i->first[2],i->first[3],
                                         custom_parameters);
    }
  }

  const std::vector<std::pair<Index<5>,std::vector<Real> > >& cmaps = top->GetCMaps();
  if(!cmaps.empty()){
    std::vector<std::vector<Real> > unique_cmap_parameters;
    std::vector<Real> current_parameters;
    std::vector<int> map_indices;
    bool cmap_present;

    for(std::vector<std::pair<Index<5>,std::vector<Real> > >::const_iterator i = cmaps.begin();
        i != cmaps.end(); ++i){
      cmap_present = false;
      current_parameters = i->second;
      for(uint j = 0; j < unique_cmap_parameters.size(); ++j){
        if(current_parameters == unique_cmap_parameters[j]){
          cmap_present = true;
          map_indices.push_back(j);
          break;
        }
      }
      if(!cmap_present){
        unique_cmap_parameters.push_back(current_parameters);
        map_indices.push_back(unique_cmap_parameters.size()-1);
      }
    } 

    OpenMM::CMAPTorsionForce& cmap_force = *new OpenMM::CMAPTorsionForce();       
    sys->addForce(&cmap_force);
    mapper[CMap] = mapper_index++;

    for(std::vector<std::vector<Real> >::iterator u_cmap = unique_cmap_parameters.begin();
        u_cmap != unique_cmap_parameters.end(); ++u_cmap){
      uint size = sqrt(u_cmap->size());
      if(size * size != u_cmap->size()){
        throw ost::Error("Expect cmap to be quadratic!");
      }

      //in charmm/gromacs the cmap has the following order regarding angle bins:
      //(phi1,psi1),(phi1,psi2)...
      //openmm requires following order:
      //(phi1,psi1),(phi2,psi1)...

      //in charmm/gromacs the range of the angles in the map is [-180,180],[-pi,pi]
      //respectively
      //openmm requires the range to be [0,2pi], the map has to be shifted...

      //=>we need basically two "transformation"

      std::vector<double> transformed_map;
      for(uint i = 0; i < size; ++i){
        for(uint j = 0; j < size; ++j){
          transformed_map.push_back((*u_cmap)[size*((j+size/2)%size) + ((i+size/2)%size)]);
        }
      }
      cmap_force.addMap(size, transformed_map);
    }

    //add the cmaps
    int atom1,atom2,atom3,atom4,atom5;
    for(uint i = 0; i < cmaps.size(); ++i){
      atom1 = cmaps[i].first[0];
      atom2 = cmaps[i].first[1];
      atom3 = cmaps[i].first[2];
      atom4 = cmaps[i].first[3];
      atom5 = cmaps[i].first[4];
      cmap_force.addTorsion(map_indices[i],atom1,atom2,atom3,atom4,atom2,atom3,atom4,atom5);
    }
  }

  const std::vector<std::pair<Index<1>,std::vector<Real> > >& harmonic_position_restraints = top->GetHarmonicPositionRestraints();
  if(!harmonic_position_restraints.empty()){
    OpenMM::CustomExternalForce& position_restraint_force = *new OpenMM::CustomExternalForce("k_force*(x_scale*(x-x0)^2+y_scale*(y-y0)^2+z_scale*(z-z0)^2)");
    position_restraint_force.addPerParticleParameter("x0");
    position_restraint_force.addPerParticleParameter("y0");
    position_restraint_force.addPerParticleParameter("z0");
    position_restraint_force.addPerParticleParameter("k_force");
    position_restraint_force.addPerParticleParameter("x_scale");
    position_restraint_force.addPerParticleParameter("y_scale");
    position_restraint_force.addPerParticleParameter("z_scale");
    sys->addForce(&position_restraint_force);
    mapper[HarmonicPositionRestraint] = mapper_index++;
    for(std::vector<std::pair<Index<1>,std::vector<Real> > >::const_iterator i = harmonic_position_restraints.begin();
        i != harmonic_position_restraints.end(); ++i){
      //stupid cast
      std::vector<double> parameters;
      for(std::vector<Real>::const_iterator j = i->second.begin();
          j != i->second.end(); ++j){
        parameters.push_back(*j);
      }
      parameters[0] *= 0.1;
      parameters[1] *= 0.1;
      parameters[2] *= 0.1;
      
      position_restraint_force.addParticle(i->first[0], parameters);
    }
  } 

  const std::vector<std::pair<Index<2>,std::vector<Real> > > harmonic_distance_restraints = top->GetHarmonicDistanceRestraints();
  if(!harmonic_distance_restraints.empty()){
    OpenMM::HarmonicBondForce& harmonic_distance_restraint_force = *new OpenMM::HarmonicBondForce();
    sys->addForce(&harmonic_distance_restraint_force);
    mapper[HarmonicDistanceRestraint] = mapper_index++;
    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::const_iterator i = harmonic_distance_restraints.begin();
        i != harmonic_distance_restraints.end(); ++i){
      //note, that the restraint formula is k(b-b0)^2 instead of 1/2*k(b-b0)^2 => multiply by two!
      harmonic_distance_restraint_force.addBond(i->first[0],i->first[1],i->second[0],2*i->second[1]);
    }
  }

  const std::vector<std::pair<Index<2>,std::vector<Real> > > fgmd_hbond_donors = top->GetFGMDHBondDonors();
  const std::vector<Index<2> > fgmd_hbond_acceptors = top->GetFGMDHBondAcceptors();
  if(!fgmd_hbond_donors.empty() && !fgmd_hbond_acceptors.empty()){
    String distance_term = "k_length*(distance(a1,d1)-k_length)^2";
    String alpha_term = "k_alpha*(angle(d2,d1,a1)-alpha)^2";
    String beta_term = "k_beta*(angle(d1,a1,a2)-beta)^2";
    String f = distance_term + "+" + alpha_term + "+" + beta_term; 
    OpenMM::CustomHbondForce& fgmd_hbond_force = *new OpenMM::CustomHbondForce(f);

    fgmd_hbond_force.addPerDonorParameter("length");
    fgmd_hbond_force.addPerDonorParameter("k_length");
    fgmd_hbond_force.addPerDonorParameter("alpha");
    fgmd_hbond_force.addPerDonorParameter("k_alpha");
    fgmd_hbond_force.addPerDonorParameter("beta");
    fgmd_hbond_force.addPerDonorParameter("k_beta");

    fgmd_hbond_force.setCutoffDistance(3.0);
    fgmd_hbond_force.setNonbondedMethod(OpenMM::CustomHbondForce::CutoffPeriodic);

    switch(settings->nonbonded_method){
      case NoCutoff:{
        fgmd_hbond_force.setNonbondedMethod(OpenMM::CustomHbondForce::NoCutoff);
        break;
      }
      case CutoffNonPeriodic:{
        fgmd_hbond_force.setNonbondedMethod(OpenMM::CustomHbondForce::CutoffNonPeriodic);
        break;        
      }
      case CutoffPeriodic:{
        fgmd_hbond_force.setNonbondedMethod(OpenMM::CustomHbondForce::CutoffPeriodic);
        break;        
      }
      case Ewald:{
        fgmd_hbond_force.setNonbondedMethod(OpenMM::CustomHbondForce::CutoffPeriodic);
        break;        
      }
      case PME:{
        fgmd_hbond_force.setNonbondedMethod(OpenMM::CustomHbondForce::CutoffPeriodic);
        break;        
      }
    }

    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::const_iterator i = fgmd_hbond_donors.begin();
        i != fgmd_hbond_donors.end(); ++i){
      //stupid cast
      std::vector<double> parameters;
      for(std::vector<Real>::const_iterator j = i->second.begin();
          j != i->second.end(); ++j){
        parameters.push_back(*j);
      }
      fgmd_hbond_force.addDonor(i->first[0],i->first[1],-1, parameters);
    }
    std::vector<double> dummy_vec;
    for(std::vector<Index<2> >::const_iterator i = fgmd_hbond_acceptors.begin();
        i != fgmd_hbond_acceptors.end(); ++i){
      fgmd_hbond_force.addAcceptor((*i)[0],(*i)[1],-1,dummy_vec);
    }
    sys->addForce(&fgmd_hbond_force);
  }

  std::vector<Real> sigmas = top->GetSigmas();
  std::vector<Real> epsilons = top->GetEpsilons();
  std::vector<Real> charges = top->GetCharges();
  if(!sigmas.empty() && !epsilons.empty() && !charges.empty()){
    OpenMM::NonbondedForce& nonbonded_force = *new OpenMM::NonbondedForce();
    switch(settings->nonbonded_method){
      case NoCutoff:{
        nonbonded_force.setNonbondedMethod(OpenMM::NonbondedForce::NoCutoff);
        break;
      }
      case CutoffNonPeriodic:{
        nonbonded_force.setNonbondedMethod(OpenMM::NonbondedForce::CutoffNonPeriodic);
        break;        
      }
      case CutoffPeriodic:{
        nonbonded_force.setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
        break;        
      }
      case Ewald:{
        nonbonded_force.setNonbondedMethod(OpenMM::NonbondedForce::Ewald);
        break;        
      }
      case PME:{
        nonbonded_force.setNonbondedMethod(OpenMM::NonbondedForce::PME);
        break;        
      }
    }
    //we can set the cutoff in any case, since it will have no effect when
    //nonbonded method is NoCutoff
    nonbonded_force.setCutoffDistance(settings->nonbonded_cutoff/10.0);
    nonbonded_force.setReactionFieldDielectric(settings->reaction_field_dielectric);
    nonbonded_force.setUseDispersionCorrection(settings->use_dispersion_correction);
    sys->addForce(&nonbonded_force);
    mapper[LJ] = mapper_index++;

    //we first add all single particles to the force
    for(uint i = 0; i < charges.size(); ++i){
      nonbonded_force.addParticle(charges[i], sigmas[i], epsilons[i]);
    }
    //take care of the 1-4 interactions
    const std::vector<std::pair<Index<2>,std::vector<Real> > >& lj_pairs = top->GetLJPairs();
    if(!lj_pairs.empty()){
      Real fudge_lj = top->GetFudgeLJ();
      Real fudge_qq = top->GetFudgeQQ();
      Real charge;
      for(std::vector<std::pair<Index<2>,std::vector<Real> > >::const_iterator i = lj_pairs.begin();
          i != lj_pairs.end(); ++i){
        charge = fudge_qq * charges[i->first[0]] * charges[i->first[1]];
        nonbonded_force.addException(i->first[0],i->first[1],charge,
                                     i->second[0],fudge_lj*i->second[1]);
      }
    }
    //take care of the exclusions
    const std::vector<Index<2> >& exclusions = top->GetExclusions();
    if(!exclusions.empty()){
      for(std::vector<Index<2> >::const_iterator i = exclusions.begin(); i != exclusions.end(); ++i){
        nonbonded_force.addException((*i)[0],(*i)[1],0.0,1.0,0.0);
      }
    }
  }
  std::vector<Real> gbsa_radii = top->GetGBSARadii();
  std::vector<Real> gbsa_scaling = top->GetOBCScalings();
  if(!gbsa_radii.empty() && !gbsa_scaling.empty()){
    std::vector<Real> charges = top->GetCharges();
    if(charges.empty())
      throw ost::Error("GBSA force requires atom charges to be set!");
    OpenMM::GBSAOBCForce& gbsa_force = *new OpenMM::GBSAOBCForce();
    sys->addForce(&gbsa_force);
    if(settings->nonbonded_method == NoCutoff){
      gbsa_force.setNonbondedMethod(OpenMM::GBSAOBCForce::NoCutoff);
    }
    else if(settings->nonbonded_method == CutoffNonPeriodic){
      gbsa_force.setNonbondedMethod(OpenMM::GBSAOBCForce::CutoffNonPeriodic);
    }
    else if(settings->nonbonded_method == CutoffPeriodic){
      gbsa_force.setNonbondedMethod(OpenMM::GBSAOBCForce::CutoffPeriodic);
    }
    else{
      throw ost::Error("Can only use Nonbonded methods NoCutoff, CutoffNonPeriodic,CutoffPeriodic when using GBSA");
    }
    mapper[GBSA] = mapper_index++;
    gbsa_force.setCutoffDistance(settings->nonbonded_cutoff/10);
    gbsa_force.setSolventDielectric(settings->solvent_dielectric);
    gbsa_force.setSoluteDielectric(settings->solute_dielectric);
    //it should be guaranteed from the topology object, that the vectors have the same size 
    for(uint i = 0; i < gbsa_radii.size(); ++i){
      gbsa_force.addParticle(charges[i],
                             gbsa_radii[i],
                             gbsa_scaling[i]);
    }
  }
  //set constraints
  std::vector<std::pair<Index<2>,std::vector<Real> > > distance_constraints = top->GetDistanceConstraints();
  if(!distance_constraints.empty()){
    for(std::vector<std::pair<Index<2>,std::vector<Real> > >::iterator i = distance_constraints.begin();
        i != distance_constraints.end(); ++i){
      sys->addConstraint(i->first[0],i->first[1],i->second[0]);
    }
  }
  //set periodic box extents
  if(geom::Length(settings->periodic_box_extents) > 0){
    OpenMM::Vec3 x_vec(settings->periodic_box_extents[0]/10,0.0,0.0);
    OpenMM::Vec3 y_vec(0.0,settings->periodic_box_extents[1]/10,0.0);
    OpenMM::Vec3 z_vec(0.0,0.0,settings->periodic_box_extents[2]/10);
    sys->setDefaultPeriodicBoxVectors(x_vec,y_vec,z_vec);
  }
  else if(settings->add_nonbonded && 
           !(settings->nonbonded_method == NoCutoff ||
            settings->nonbonded_method == CutoffNonPeriodic)){
    throw ost::Error("Chosen nonbonded method requires to define the periodic box extents!");
  }
  //set the CMMMotion removing force if required
  if(settings->remove_cmm_motion){
    if(settings->cmm_frequency != settings->cmm_frequency){
      throw ost::Error("Cannot set CMM Remover without valid cmm_frequency set!");
    }
    sys->addForce(new OpenMM::CMMotionRemover(settings->cmm_frequency));
  }
  if(settings->add_thermostat){
    if(settings->thermostat_temperature != settings->thermostat_temperature ||
       settings->thermostat_collision_frequency != settings->thermostat_collision_frequency){
      throw ost::Error("Cannot set thermostat without defined temperature and collision frequency!");
    }
    OpenMM::AndersenThermostat& thermostat = *new OpenMM::AndersenThermostat(settings->thermostat_temperature, 
                                                                             settings->thermostat_collision_frequency);
    sys->addForce(&thermostat);
  }
  if(settings->add_barostat){
    if(geom::Length(settings->periodic_box_extents) == 0){
      throw ost::Error("If you set a barostat, you also have to set nonzero periodic box extents!");
    }
    if(settings->barostat_temperature != settings->barostat_temperature ||
       settings->barostat_pressure != settings->barostat_pressure ||
       settings->barostat_frequency != settings->barostat_frequency){
      throw ost::Error("Cannot set barostat without defined temperature, pressure, and frequency parameters!");
    }
    OpenMM::MonteCarloBarostat& barostat = *new OpenMM::MonteCarloBarostat(settings->barostat_pressure,
                                                                           settings->barostat_temperature,
                                                                           settings->barostat_frequency);
    sys->addForce(&barostat);
  }

  return sys;
}
}}} //ns
