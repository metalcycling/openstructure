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

#include <boost/test/unit_test.hpp>
#include <ost/mol/mm/settings.hh>
#include <ost/mol/mm/forcefield.hh>
#include <ost/mol/mm/interaction.hh>
#include <ost/mol/mm/buildingblock.hh>
#include <ost/mol/mm/gromacs_block_modifiers.hh>
#include <ost/message.hh>
#include <ost/mol/builder.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/mol/mm/simulation.hh>
#include <ost/mol/mm/topology.hh>
#include <OpenMM.h>


using namespace ost::mol::mm;

BOOST_AUTO_TEST_SUITE( mol_mm );

BOOST_AUTO_TEST_CASE(test_simulation_basics){

  String pdb_name = "1CRN.pdb";

  ost::io::PDBReader reader(pdb_name, ost::io::IOProfile());
  ost::mol::EntityHandle test_ent = ost::mol::CreateEntity();
  reader.Import(test_ent);

  ost::conop::ProcessorPtr processor(new ost::conop::HeuristicProcessor);
  processor->Process(test_ent);

  SettingsPtr settings(new Settings);
  ForcefieldPtr forcefield = Forcefield::Load("CHARMM27.dat");
  settings->forcefield = forcefield;
  settings->add_gbsa = true;

  TopologyPtr top = TopologyCreator::Create(test_ent,settings);
  //lets add interactions, that are not set when using charmm
  top->AddHarmonicAngle(1,2,3,10.0,10.0);
  top->AddPeriodicImproper(1,2,3,4,2,10.0,10.0);
  top->AddDistanceConstraint(1,2,10.0);
  top->AddHarmonicPositionRestraint(1,geom::Vec3(0.0,0.0,0.0),10.0);
  top->AddHarmonicDistanceRestraint(1,2,10.0,10.0);

  settings->integrator = IntegratorPtr(new OpenMM::VerletIntegrator(0.002));

  Simulation sim(top, test_ent, settings);


  //we check, wether the reset functions have the desired effect on the topology
  //we cannot really check the effect on the openmm system...
  sim.ResetHarmonicBond(42,5.0,6.0);
  sim.ResetHarmonicAngle(0,8.0,9.0);
  sim.ResetUreyBradleyAngle(42,5.0,6.0,7.0,8.0);
  sim.ResetPeriodicDihedral(42,4,5.0,6.0);
  sim.ResetPeriodicImproper(0,5,8.0,9.0);
  sim.ResetHarmonicImproper(42,5.0,6.0);
  sim.ResetLJPair(42,5.0,6.0);
  sim.ResetDistanceConstraint(0,10.0);
  sim.ResetHarmonicPositionRestraint(0,geom::Vec3(1.0,1.0,1.0),5.0);
  sim.ResetHarmonicDistanceRestraint(0,5.0,10.0);
  sim.ResetLJ(42,2.0,3.0);
  sim.ResetGBSA(42,5.0,6.0);
  sim.ResetCharge(42,100.0);
  sim.ResetMass(42,200.0);

  uint ui1,ui2,ui3,ui4;
  int i1;
  Real r1,r2,r3,r4;
  geom::Vec3 pos;


  top->GetHarmonicBondParameters(42,ui1,ui2,r1,r2);
  BOOST_CHECK(r1 == 5.0 && r2 == 6.0);

  top->GetHarmonicAngleParameters(0,ui1,ui2,ui3,r1,r2);
  BOOST_CHECK(r1 == 8.0 && r2 == 9.0);

  top->GetUreyBradleyAngleParameters(42,ui1,ui2,ui3,r1,r2,r3,r4);
  BOOST_CHECK(r1 == 5.0 && r2 == 6.0 && r3 == 7.0 && r4 == 8.0);

  top->GetPeriodicDihedralParameters(42,ui1,ui2,ui3,ui4,i1,r1,r2);
  BOOST_CHECK(i1 == 4 && r1 == 5.0 && r2 == 6.0);

  top->GetPeriodicImproperParameters(0,ui1,ui2,ui3,ui4,i1,r1,r2);
  BOOST_CHECK(i1 == 5 && r1 == 8.0 && r2 == 9.0);

  top->GetLJPairParameters(42,ui1,ui2,r1,r2);
  BOOST_CHECK(r1 == 5.0 && r2 == 6.0);

  top->GetDistanceConstraintParameters(0,ui1,ui2,r1);
  BOOST_CHECK(r1 == 10.0);

  top->GetHarmonicPositionRestraintParameters(0, ui1, pos,r1, r2, r3, r3);
  BOOST_CHECK(pos == geom::Vec3(1.0,1.0,1.0) && r1 == 5.0);

  top->GetHarmonicDistanceRestraintParameters(0,ui1,ui2, r1, r2);
  BOOST_CHECK(r1 == 5.0 && r2 == 10.0);

  r1 = top->GetSigma(42);
  r2 = top->GetEpsilon(42);
  BOOST_CHECK(r1 == 2.0 && r2 == 3.0);

  r1 = top->GetGBSARadius(42);
  r2 = top->GetOBCScaling(42);
  BOOST_CHECK(r1 == 5.0 && r2 == 6.0);

  r1 = top->GetCharge(42);
  r2 = top->GetMass(42);
  BOOST_CHECK(r1 == 100.0 && r2 == 200.0);
}


BOOST_AUTO_TEST_CASE(test_simulation_energy_calculations){


  String pdb_name = "1CRN.pdb";

  ost::io::PDBReader reader(pdb_name, ost::io::IOProfile());
  ost::mol::EntityHandle test_ent = ost::mol::CreateEntity();
  reader.Import(test_ent);

  ost::conop::ProcessorPtr processor(new ost::conop::HeuristicProcessor);
  processor->Process(test_ent);

  SettingsPtr settings(new Settings);
  ForcefieldPtr forcefield = Forcefield::Load("CHARMM27.dat");
  settings->forcefield = forcefield;

  TopologyPtr top = TopologyCreator::Create(test_ent,settings);
  //lets add interactions, that are not set when using charmm

  settings->integrator = IntegratorPtr(new OpenMM::VerletIntegrator(0.002));

  Simulation sim(top, test_ent,settings);

  sim.ApplySD(1.0, 200);

  //the energy values have been compared with the results from the
  //CHARMM27 forcefield in gromacs

  //differenecs in energy(kJ/mol):
  //harmonic_bond:      0.00021344800001
  //urey_bradley_angle: 0.00021319300003
  //periodic_dihedral:  0.00018275799993
  //harmonic_improper:  0.00001410630000
  //cmap:               0.00192255100001
  //nonbonded:          0.00204422999968

  //You'll most likely get different results when repeating the comparison with gromacs
  //this is because of the precision.
  //Right now we have full real precision in the positions.
  //When you write down the pdb, precision will be lowered to 3 digits.
  //Another digit gets lost, when you create a gromacs topology...
  //Gromacs writes down the positions using nm instead of A, but
  //also has a precision of 3 digits.

  // Threshold is relatively permissive (1e-1->0.1% difference) depending on
  // OpenMM version / Hardware you might get slightly different results. 
  // Error propagates when performing 200 steps steepest descent (minimizes from
  // roughly 1981 down to -2381).
  BOOST_CHECK_CLOSE(sim.GetPotentialEnergy(), -2381.24, Real(1e-1));
}



BOOST_AUTO_TEST_SUITE_END();
