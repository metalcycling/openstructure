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
#include <ost/mol/mm/topology.hh>
#include <ost/mol/mm/topology_creator.hh>
#include <ost/mol/mm/settings.hh>
#include <ost/mol/mm/forcefield.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/conop/heuristic.hh>
#include <ost/message.hh>
#include <ost/mol/builder.hh>
#include <ost/mol/xcs_editor.hh>


using namespace ost::mol::mm;

BOOST_AUTO_TEST_SUITE( mol_mm );

BOOST_AUTO_TEST_CASE(test_topology_basics){

  String pdb_name = "1CRN.pdb";

  ost::io::PDBReader reader(pdb_name, ost::io::IOProfile());
  ost::mol::EntityHandle test_ent = ost::mol::CreateEntity();
  reader.Import(test_ent);

  ost::conop::ProcessorPtr processor(new ost::conop::HeuristicProcessor);
  processor->Process(test_ent);

  //check initialisation without settings
  std::vector<Real> zero_real(0);
  std::vector<Real> lot_of_reals(test_ent.GetAtomCount());

  TopologyPtr top(new Topology(lot_of_reals));

  uint ui1(0), ui2(0), ui3(0), ui4(0), ui5(0);
  int i1(0), i5(0);  
  Real r1(0.0), r2(0.0), r3(0.0), r4(0.0);
  geom::Vec3 zero_vec(0.0,0.0,0.0);
  geom::Vec3 test_vec(0.0,1.0,2.0);

  //check wether errors geth thrown, when atom indices are too big
  BOOST_CHECK_THROW(top->AddHarmonicBond(700,1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicBond(0,700,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicAngle(700,1,1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicAngle(1,700,1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicAngle(1,1,700,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddUreyBradleyAngle(700,1,1,r1,r2,r3,r4),ost::Error);
  BOOST_CHECK_THROW(top->AddUreyBradleyAngle(1,700,1,r1,r2,r3,r4),ost::Error);
  BOOST_CHECK_THROW(top->AddUreyBradleyAngle(1,1,700,r1,r2,r3,r4),ost::Error);
  BOOST_CHECK_THROW(top->AddPeriodicDihedral(700,1,1,1,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddPeriodicDihedral(1,700,1,1,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddPeriodicDihedral(1,1,700,1,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddPeriodicDihedral(1,1,1,700,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddPeriodicImproper(700,1,1,1,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddPeriodicImproper(1,700,1,1,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddPeriodicImproper(1,1,700,1,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddPeriodicImproper(1,1,1,700,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicImproper(700,1,1,1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicImproper(1,700,1,1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicImproper(1,1,700,1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicImproper(1,1,1,700,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddCMap(700,1,1,1,1,0,zero_real),ost::Error);
  BOOST_CHECK_THROW(top->AddCMap(1,700,1,1,1,0,zero_real),ost::Error);
  BOOST_CHECK_THROW(top->AddCMap(1,1,700,1,1,0,zero_real),ost::Error);
  BOOST_CHECK_THROW(top->AddCMap(1,1,1,700,1,0,zero_real),ost::Error);
  BOOST_CHECK_THROW(top->AddCMap(1,1,1,1,700,0,zero_real),ost::Error);
  BOOST_CHECK_THROW(top->AddCMap(1,1,1,1,1,0,lot_of_reals),ost::Error); //additionally checks error, when dimension and cmap are inconsistent
  BOOST_CHECK_THROW(top->AddLJPair(1,700,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddLJPair(700,1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddExclusion(1,700),ost::Error);
  BOOST_CHECK_THROW(top->AddExclusion(700,1),ost::Error);
  BOOST_CHECK_THROW(top->AddDistanceConstraint(1,700,r1),ost::Error);
  BOOST_CHECK_THROW(top->AddDistanceConstraint(700,1,r1),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicPositionRestraint(700,zero_vec,r1,r2,r3,r4),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicDistanceRestraint(700,1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->AddHarmonicDistanceRestraint(1,700,r1,r2),ost::Error);

  //check whether errors get thrown when single atom parameters are inconsistent with number of atoms in topology
  BOOST_CHECK_THROW(top->SetSigmas(zero_real),ost::Error);
  BOOST_CHECK_THROW(top->SetEpsilons(zero_real),ost::Error);
  BOOST_CHECK_THROW(top->SetGBSARadii(zero_real),ost::Error);
  BOOST_CHECK_THROW(top->SetOBCScalings(zero_real),ost::Error);
  BOOST_CHECK_THROW(top->SetMasses(zero_real),ost::Error);
  BOOST_CHECK_THROW(top->SetCharges(zero_real),ost::Error);
  
  //Single atom parameters are not allowed to be modified, without being set initially
  BOOST_CHECK_THROW(top->SetSigma(1,r1),ost::Error);
  BOOST_CHECK_THROW(top->SetEpsilon(1,r1),ost::Error);
  BOOST_CHECK_THROW(top->SetGBSARadius(1,r1),ost::Error);
  BOOST_CHECK_THROW(top->SetOBCScaling(1,r1),ost::Error);
  BOOST_CHECK_THROW(top->SetCharge(1,r1),ost::Error); 

  //Let's set the single atom parameters initially
  top->SetSigmas(lot_of_reals);
  top->SetEpsilons(lot_of_reals);
  top->SetGBSARadii(lot_of_reals);
  top->SetOBCScalings(lot_of_reals);
  top->SetCharges(lot_of_reals);  

  //do the index checks for them
  BOOST_CHECK_THROW(top->SetSigma(1000,1.0),ost::Error);
  BOOST_CHECK_THROW(top->SetEpsilon(1000,1.0),ost::Error);
  BOOST_CHECK_THROW(top->SetGBSARadius(1000,1.0),ost::Error);
  BOOST_CHECK_THROW(top->SetOBCScaling(1000,1.0),ost::Error);
  BOOST_CHECK_THROW(top->SetMass(1000,1.0),ost::Error);
  BOOST_CHECK_THROW(top->SetCharge(1000,1.0),ost::Error);   

  //similar stupid checks for SetParameter functions
  BOOST_CHECK_THROW(top->SetHarmonicBondParameters(700,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->SetHarmonicAngleParameters(3000,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->SetUreyBradleyAngleParameters(3000,r1,r2,r3,r4),ost::Error);
  BOOST_CHECK_THROW(top->SetPeriodicDihedralParameters(3000,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->SetPeriodicImproperParameters(3000,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->SetHarmonicImproperParameters(3000,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->SetCMapParameters(3000,0,zero_real),ost::Error);
  BOOST_CHECK_THROW(top->SetLJPairParameters(3000,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->SetDistanceConstraintParameters(3000,r1),ost::Error);
  BOOST_CHECK_THROW(top->SetHarmonicPositionRestraintParameters(3000,zero_vec,r1),ost::Error);
  BOOST_CHECK_THROW(top->SetHarmonicDistanceRestraintParameters(3000,r1,r2),ost::Error);

  BOOST_CHECK_THROW(top->GetHarmonicBondParameters(700,ui1,ui2,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->GetHarmonicAngleParameters(3000,ui1,ui2,ui3,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->GetUreyBradleyAngleParameters(3000,ui1,ui2,ui3,r1,r2,r3,r4),ost::Error);
  BOOST_CHECK_THROW(top->GetPeriodicDihedralParameters(3000,ui1,ui2,ui3,ui4,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->GetPeriodicImproperParameters(3000,ui1,ui2,ui3,ui4,i1,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->GetHarmonicImproperParameters(3000,ui1,ui2,ui3,ui4,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->GetCMapParameters(3000,ui1,ui2,ui3,ui4,ui5,i1,zero_real),ost::Error);
  BOOST_CHECK_THROW(top->GetLJPairParameters(3000,ui1,ui2,r1,r2),ost::Error);
  BOOST_CHECK_THROW(top->GetDistanceConstraintParameters(3000,ui1,ui2,r1),ost::Error);
  BOOST_CHECK_THROW(top->GetHarmonicPositionRestraintParameters(3000,ui1,zero_vec,r1,r2,r3,r4),ost::Error);
  BOOST_CHECK_THROW(top->GetHarmonicDistanceRestraintParameters(3000,ui1,ui2,r1,r2),ost::Error);

  //let's add some interactions to the topology
  for(int i = 0; i < 10; ++i){
    top->AddHarmonicBond(0,1,r1,r2);
    top->AddHarmonicAngle(0,1,2,r1,r2);
    top->AddUreyBradleyAngle(0,1,2,r1,r2,r3,r4);
    top->AddPeriodicDihedral(0,1,2,3,i1,r1,r2);
    top->AddPeriodicImproper(0,1,2,3,i1,r1,r2);
    top->AddHarmonicImproper(0,1,2,3,r1,r2);
    top->AddCMap(0,1,2,3,4,i1,zero_real);
    top->AddHarmonicPositionRestraint(0,zero_vec,r1,r2,r3,r4);
    top->AddHarmonicDistanceRestraint(0,1,r1,r2);
  }

  //following interactions are unique for a pair of atoms
  top->AddLJPair(0,1,r1,r2);
  top->AddDistanceConstraint(0,1,r1);
  top->AddExclusion(0,1);

  //errors should be thrown, when previously added interactions are added again
  BOOST_CHECK_THROW(top->AddLJPair(0,1,r1,r2), ost::Error);
  BOOST_CHECK_THROW(top->AddDistanceConstraint(0,1,r1), ost::Error);
  BOOST_CHECK_THROW(top->AddExclusion(0,1), ost::Error);
  BOOST_CHECK_THROW(top->AddLJPair(1,0,r1,r2), ost::Error);
  BOOST_CHECK_THROW(top->AddDistanceConstraint(1,0,r1), ost::Error);
  BOOST_CHECK_THROW(top->AddExclusion(1,0), ost::Error);

  //check set and get parameter stuff, but also the GetIndices functionality
  std::vector<uint> ind1,ind2;
  top->SetHarmonicBondParameters(0,0.5,15.0);
  top->GetHarmonicBondParameters(0,ui1,ui2,r1,r2);
  BOOST_CHECK(r1 == 0.5 && r2 == 15.0);
  ind1 = top->GetHarmonicBondIndices(0,1);
  ind2 = top->GetHarmonicBondIndices(1,0);
  BOOST_CHECK(ind1 == ind2 && ind1.size() == 10);
  BOOST_CHECK(top->GetHarmonicBondIndices(0).size() == 10);
  BOOST_CHECK(top->GetHarmonicBondIndices(1).size() == 10);
  BOOST_CHECK(top->GetHarmonicBondIndices(2).size() == 0);


  top->SetHarmonicAngleParameters(0,1.0,10.0);
  top->GetHarmonicAngleParameters(0,ui1,ui2,ui3,r1,r2);
  BOOST_CHECK(r1 == 1.0 && r2 == 10.0);
  ind1 = top->GetHarmonicAngleIndices(0,1,2);
  ind2 = top->GetHarmonicAngleIndices(2,1,0);
  BOOST_CHECK(ind1 == ind2 && ind1.size() == 10);
  BOOST_CHECK(top->GetHarmonicAngleIndices(0).size() == 10);
  BOOST_CHECK(top->GetHarmonicAngleIndices(1).size() == 10);
  BOOST_CHECK(top->GetHarmonicAngleIndices(2).size() == 10);
  BOOST_CHECK(top->GetHarmonicAngleIndices(3).size() == 0);

  top->SetUreyBradleyAngleParameters(0,0.5,0.5,0.5,0.5);
  top->GetUreyBradleyAngleParameters(0,ui1,ui2,ui3,r1,r2,r3,r4);
  BOOST_CHECK( r1 == 0.5 && r2 == 0.5 && r3 == 0.5 && r3 == 0.5 );
  ind1 = top->GetUreyBradleyAngleIndices(0,1,2);
  ind2 = top->GetUreyBradleyAngleIndices(2,1,0);
  BOOST_CHECK(ind1 == ind2 && ind1.size() == 10);
  BOOST_CHECK(top->GetUreyBradleyAngleIndices(0).size() == 10);
  BOOST_CHECK(top->GetUreyBradleyAngleIndices(1).size() == 10);
  BOOST_CHECK(top->GetUreyBradleyAngleIndices(2).size() == 10);
  BOOST_CHECK(top->GetUreyBradleyAngleIndices(3).size() == 0);


  top->SetPeriodicDihedralParameters(0,2,5.0,0.5);
  top->GetPeriodicDihedralParameters(0,ui1,ui2,ui3,ui4,i1,r1,r2);
  BOOST_CHECK( i5 == 0 && r1 == 5.0 && r2 == 0.5 );
  ind1 = top->GetPeriodicDihedralIndices(0,1,2,3);
  ind2 = top->GetPeriodicDihedralIndices(3,2,1,0);
  BOOST_CHECK(ind1 == ind2 && ind1.size() == 10);
  BOOST_CHECK(top->GetPeriodicDihedralIndices(0).size() == 10);
  BOOST_CHECK(top->GetPeriodicDihedralIndices(1).size() == 10);
  BOOST_CHECK(top->GetPeriodicDihedralIndices(2).size() == 10);
  BOOST_CHECK(top->GetPeriodicDihedralIndices(3).size() == 10);
  BOOST_CHECK(top->GetPeriodicDihedralIndices(4).size() == 0);

  top->SetPeriodicImproperParameters(6,2,2.0,1.0);
  top->GetPeriodicImproperParameters(6,ui1,ui2,ui3,ui4,i1,r1,r2);
  BOOST_CHECK( i5 == 0 && r1 == 2.0 && r2 == 1.0 );
  ind1 = top->GetPeriodicImproperIndices(0,1,2,3);
  ind2 = top->GetPeriodicImproperIndices(3,2,1,0);
  BOOST_CHECK(ind1 == ind2 && ind1.size() == 10);
  BOOST_CHECK(top->GetPeriodicImproperIndices(0).size() == 10);
  BOOST_CHECK(top->GetPeriodicImproperIndices(1).size() == 10);
  BOOST_CHECK(top->GetPeriodicImproperIndices(2).size() == 10);
  BOOST_CHECK(top->GetPeriodicImproperIndices(3).size() == 10);
  BOOST_CHECK(top->GetPeriodicImproperIndices(4).size() == 0);

  top->SetHarmonicImproperParameters(0,4.0,4.0);
  top->GetHarmonicImproperParameters(0,ui1,ui2,ui3,ui4,r1,r2);
  BOOST_CHECK(r1 == 4.0 && r2 == 4.0);
  ind1 = top->GetHarmonicImproperIndices(0,1,2,3);
  ind2 = top->GetHarmonicImproperIndices(3,2,1,0);
  BOOST_CHECK(ind1 == ind2 && ind1.size() == 10);
  BOOST_CHECK(top->GetHarmonicImproperIndices(0).size() == 10);
  BOOST_CHECK(top->GetHarmonicImproperIndices(1).size() == 10);
  BOOST_CHECK(top->GetHarmonicImproperIndices(2).size() == 10);
  BOOST_CHECK(top->GetHarmonicImproperIndices(3).size() == 10);
  BOOST_CHECK(top->GetHarmonicImproperIndices(4).size() == 0);

  std::vector<Real> one_real(1);
  one_real[0] = 42;
  top->SetCMapParameters(0,1,one_real);
  std::vector<Real> temp_real;
  top->GetCMapParameters(0,ui1,ui2,ui3,ui4,ui5,i1,temp_real);
  BOOST_CHECK( ui1 == 0 && temp_real.size() == 1 && temp_real[0] == 42);
  ind1 = top->GetCMapIndices(0,1,2,3,4);
  ind2 = top->GetCMapIndices(4,3,2,1,0);
  BOOST_CHECK(ind1.size() == 10);
  BOOST_CHECK(ind2.size() == 0); //cmaps are sensitive for the ordering!!
  BOOST_CHECK(top->GetCMapIndices(0).size() == 10);
  BOOST_CHECK(top->GetCMapIndices(1).size() == 10);
  BOOST_CHECK(top->GetCMapIndices(2).size() == 10);
  BOOST_CHECK(top->GetCMapIndices(3).size() == 10);
  BOOST_CHECK(top->GetCMapIndices(4).size() == 10);
  BOOST_CHECK(top->GetCMapIndices(5).size() == 0);

  top->SetLJPairParameters(0,3.5,4.5);
  top->GetLJPairParameters(0,ui1,ui2,r1,r2);
  BOOST_CHECK( r1 == 3.5 && r2 == 4.5 );
  i1 = top->GetLJPairIndex(0,1);
  i5 = top->GetLJPairIndex(1,0);
  BOOST_CHECK(i1 == i5 && i1 != -1);
  BOOST_CHECK(top->GetLJPairIndices(0).size() == 1);
  BOOST_CHECK(top->GetLJPairIndices(1).size() == 1);
  BOOST_CHECK(top->GetLJPairIndices(2).size() == 0);

  top->SetDistanceConstraintParameters(0,42.0);
  top->GetDistanceConstraintParameters(0,ui1,ui2,r1);
  BOOST_CHECK( r1 == 42.0 );
  i1 = top->GetDistanceConstraintIndex(0,1);
  i5 = top->GetDistanceConstraintIndex(1,0);
  BOOST_CHECK(i1 == i5 && i1 != -1);

  top->SetHarmonicPositionRestraintParameters(0,test_vec,3.0);
  top->GetHarmonicPositionRestraintParameters(0,ui1,zero_vec,r1,r2,r3,r4);
  BOOST_CHECK(zero_vec == test_vec);
  BOOST_CHECK(r1 == 3.0);
  BOOST_CHECK(r2 == 1.0); //should be the same for r3 and r4
  //no reverse checking necessary
  BOOST_CHECK(top->GetHarmonicPositionRestraintIndices(0).size() == 10);
  BOOST_CHECK(top->GetHarmonicPositionRestraintIndices(1).size() == 0);

  top->SetHarmonicDistanceRestraintParameters(0,2.0,2.0);
  top->GetHarmonicDistanceRestraintParameters(0,ui1,ui2,r1,r2);
  BOOST_CHECK(r1 == 2.0 && r2 == 2.0);
  ind1 = top->GetHarmonicDistanceRestraintIndices(0,1);
  ind2 = top->GetHarmonicDistanceRestraintIndices(1,0);
  BOOST_CHECK(ind1 == ind2 && ind1.size() == 10);
  BOOST_CHECK(top->GetHarmonicDistanceRestraintIndices(0).size() == 10);
  BOOST_CHECK(top->GetHarmonicDistanceRestraintIndices(1).size() == 10);
  BOOST_CHECK(top->GetHarmonicDistanceRestraintIndices(2).size() == 0);

  top->AddPositionConstraint(5);
  BOOST_CHECK(top->GetNumPositionConstraints() == 1);
  top->ResetPositionConstraints();
  BOOST_CHECK(top->GetNumPositionConstraints() == 0);
}

BOOST_AUTO_TEST_CASE(test_topology_merge){

  String pdb_name = "1CRN.pdb";
  ost::io::PDBReader reader(pdb_name, ost::io::IOProfile());
  ost::mol::EntityHandle test_ent = ost::mol::CreateEntity();
  ost::mol::EntityHandle temp;

  reader.Import(test_ent);

  ForcefieldPtr ff = Forcefield::Load("CHARMM27.dat");
  SettingsPtr settings = SettingsPtr(new Settings);
  settings->add_gbsa = true;
  settings->forcefield = ff;

  ost::mol::EntityHandle new_ent = ost::mol::Builder()
                                   .Chain("A")
                                     .Residue("ONE")
                                       .Atom("A",geom::Vec3(0,0,0))
                                       .Atom("B",geom::Vec3(1,1,1))
                                     .Residue("TWO")
                                       .Atom("A",geom::Vec3(2,2,2))
                                       .Atom("B",geom::Vec3(3,3,3))
                                     .Residue("THREE")
                                       .Atom("A",geom::Vec3(4,4,4))
                                       .Atom("B",geom::Vec3(5,5,5));

  //check whether error gets thrown when chain with same name is already present
  std::vector<Real> real_vec(6);
  TopologyPtr merge_top_one(new Topology(real_vec));
  merge_top_one->SetCharges(real_vec);
  merge_top_one->SetOBCScalings(real_vec);
  merge_top_one->SetGBSARadii(real_vec);
  merge_top_one->SetSigmas(real_vec);
  merge_top_one->SetEpsilons(real_vec);
  merge_top_one->SetFudgeQQ(1.0);
  merge_top_one->SetFudgeLJ(1.0);
  TopologyPtr top_one = TopologyCreator::Create(test_ent,settings);
  temp = test_ent.Copy();
  BOOST_CHECK_THROW(top_one->Merge(temp,merge_top_one,new_ent),ost::Error); 
  ost::mol::XCSEditor ed = new_ent.EditXCS();
  ed.RenameChain(new_ent.GetChainList()[0],"B");
  TopologyPtr merge_top_one_two(new Topology(real_vec)); // I know...
  merge_top_one_two->SetCharges(real_vec);
  merge_top_one_two->SetOBCScalings(real_vec);
  merge_top_one_two->SetGBSARadii(real_vec);
  merge_top_one_two->SetSigmas(real_vec);
  merge_top_one_two->SetEpsilons(real_vec);
  merge_top_one_two->SetFudgeQQ(1.0);
  merge_top_one_two->SetFudgeLJ(1.0);
  BOOST_CHECK_NO_THROW(top_one->Merge(temp,merge_top_one_two,new_ent));


  //check whether error gets thrown, when fudge parameters are inconsistent
  TopologyPtr merge_top_two(new Topology(real_vec));
  merge_top_two->SetCharges(real_vec);
  merge_top_two->SetOBCScalings(real_vec);
  merge_top_two->SetGBSARadii(real_vec);
  merge_top_two->SetSigmas(real_vec);
  merge_top_two->SetEpsilons(real_vec);
  merge_top_two->SetFudgeQQ(42.0);
  merge_top_two->SetFudgeLJ(42.0);
  TopologyPtr top_two = TopologyCreator::Create(test_ent,settings);
  temp = test_ent.Copy();
  BOOST_CHECK_THROW(top_two->Merge(temp,merge_top_two,new_ent),ost::Error); 
  merge_top_two->SetFudgeQQ(1.0);
  merge_top_two->SetFudgeLJ(1.0);
  BOOST_CHECK_NO_THROW(top_two->Merge(temp,merge_top_two,new_ent));


  //check whether error gets thrown when charges are not set
  TopologyPtr merge_top_three(new Topology(real_vec));
  merge_top_three->SetOBCScalings(real_vec);
  merge_top_three->SetGBSARadii(real_vec);
  merge_top_three->SetSigmas(real_vec);
  merge_top_three->SetEpsilons(real_vec);
  merge_top_three->SetFudgeQQ(1.0);
  merge_top_three->SetFudgeLJ(1.0);
  TopologyPtr top_three = TopologyCreator::Create(test_ent,settings);
  temp = test_ent.Copy();
  BOOST_CHECK_THROW(top_three->Merge(temp,merge_top_three,new_ent),ost::Error); 
  TopologyPtr top_three_two = TopologyCreator::Create(test_ent,settings);
  merge_top_three->SetCharges(real_vec);
  BOOST_CHECK_NO_THROW(top_three_two->Merge(temp,merge_top_three,new_ent));

  //check whether error gets thrown when obc scaling factors are not set
  TopologyPtr merge_top_four(new Topology(real_vec));
  merge_top_four->SetCharges(real_vec);
  merge_top_four->SetGBSARadii(real_vec);
  merge_top_four->SetSigmas(real_vec);
  merge_top_four->SetEpsilons(real_vec);
  merge_top_four->SetFudgeQQ(1.0);
  merge_top_four->SetFudgeLJ(1.0);
  TopologyPtr top_four = TopologyCreator::Create(test_ent,settings);
  temp = test_ent.Copy();
  BOOST_CHECK_THROW(top_four->Merge(temp,merge_top_four,new_ent),ost::Error); 
  TopologyPtr top_four_two = TopologyCreator::Create(test_ent,settings);
  merge_top_four->SetOBCScalings(real_vec);
  BOOST_CHECK_NO_THROW(top_four_two->Merge(temp,merge_top_four,new_ent));

  //check whether error gets thrown when gbsa radii are not set
  TopologyPtr merge_top_five(new Topology(real_vec));
  merge_top_five->SetCharges(real_vec);
  merge_top_five->SetOBCScalings(real_vec);
  merge_top_five->SetSigmas(real_vec);
  merge_top_five->SetEpsilons(real_vec);
  merge_top_five->SetFudgeQQ(1.0);
  merge_top_five->SetFudgeLJ(1.0);
  TopologyPtr top_five = TopologyCreator::Create(test_ent,settings);
  temp = test_ent.Copy();
  BOOST_CHECK_THROW(top_five->Merge(temp,merge_top_five,new_ent),ost::Error); 
  TopologyPtr top_five_two = TopologyCreator::Create(test_ent,settings);
  merge_top_five->SetGBSARadii(real_vec);
  BOOST_CHECK_NO_THROW(top_five_two->Merge(temp,merge_top_five,new_ent));

  //check whether error gets thrown when sigmas are not set
  TopologyPtr merge_top_six(new Topology(real_vec));
  merge_top_six->SetCharges(real_vec);
  merge_top_six->SetOBCScalings(real_vec);
  merge_top_six->SetGBSARadii(real_vec);  
  merge_top_six->SetEpsilons(real_vec);
  merge_top_six->SetFudgeQQ(1.0);
  merge_top_six->SetFudgeLJ(1.0);
  TopologyPtr top_six = TopologyCreator::Create(test_ent,settings);
  temp = test_ent.Copy();
  BOOST_CHECK_THROW(top_six->Merge(temp,merge_top_six,new_ent),ost::Error); 
  TopologyPtr top_six_two = TopologyCreator::Create(test_ent,settings);
  merge_top_six->SetSigmas(real_vec);
  BOOST_CHECK_NO_THROW(top_six_two->Merge(temp,merge_top_six,new_ent));

  //check whether error gets thrown when epsilons are not set
  TopologyPtr merge_top_seven(new Topology(real_vec));
  merge_top_seven->SetCharges(real_vec);
  merge_top_seven->SetOBCScalings(real_vec);
  merge_top_seven->SetGBSARadii(real_vec); 
  merge_top_seven->SetSigmas(real_vec); 
  merge_top_seven->SetFudgeQQ(1.0);
  merge_top_seven->SetFudgeLJ(1.0);
  TopologyPtr top_seven = TopologyCreator::Create(test_ent,settings);
  temp = test_ent.Copy();
  BOOST_CHECK_THROW(top_seven->Merge(temp,merge_top_seven,new_ent),ost::Error); 
  TopologyPtr top_seven_two = TopologyCreator::Create(test_ent,settings);
  merge_top_seven->SetEpsilons(real_vec);
  BOOST_CHECK_NO_THROW(top_seven_two->Merge(temp,merge_top_seven,new_ent));

  TopologyPtr merge_top_eight(new Topology(real_vec));
  merge_top_eight->SetCharges(real_vec);
  merge_top_eight->SetOBCScalings(real_vec);
  merge_top_eight->SetGBSARadii(real_vec); 
  merge_top_eight->SetSigmas(real_vec); 
  merge_top_eight->SetFudgeQQ(1.0);
  merge_top_eight->SetFudgeLJ(1.0);
  merge_top_eight->SetEpsilons(real_vec);

  //let's add every possible interaction
  merge_top_eight->AddHarmonicBond(0,1,42.0,42000.0);
  merge_top_eight->AddHarmonicAngle(0,1,2,24.0,24000.0);
  merge_top_eight->AddUreyBradleyAngle(1,2,3,1.0,2.0,3.0,4.0);
  merge_top_eight->AddPeriodicDihedral(0,1,2,3,10,2.0,42.0);
  merge_top_eight->AddPeriodicImproper(1,2,3,4,3,5.0,10.0);
  merge_top_eight->AddHarmonicImproper(0,1,2,3,8.0,12.0);
  std::vector<Real> cmap;
  cmap.push_back(0.0);
  cmap.push_back(1.0);
  cmap.push_back(2.0);
  cmap.push_back(3.0);
  merge_top_eight->AddCMap(0,1,2,3,4,2,cmap);
  merge_top_eight->AddLJPair(3,4,1.0,2.0);
  merge_top_eight->AddDistanceConstraint(2,3,10.0);
  merge_top_eight->AddExclusion(3,4);
  merge_top_eight->AddPositionConstraint(2);
  geom::Vec3 pos(1.0,2.0,3.0);
  merge_top_eight->AddHarmonicPositionRestraint(2,pos,10.0,5.0,6.0,7.0);
  merge_top_eight->AddHarmonicDistanceRestraint(1,2,10.0,1000.0);

  TopologyPtr top_eight = TopologyCreator::Create(test_ent,settings);

  uint num_harmonic_bonds = top_eight->GetNumHarmonicBonds();
  uint num_harmonic_angles = top_eight->GetNumHarmonicAngles();
  uint num_urey_bradley_angles = top_eight->GetNumUreyBradleyAngles();
  uint num_periodic_dihedrals = top_eight->GetNumPeriodicDihedrals();
  uint num_periodic_impropers = top_eight->GetNumPeriodicImpropers();
  uint num_harmonic_impropers = top_eight->GetNumHarmonicImpropers();
  uint num_cmaps = top_eight->GetNumCMaps();
  uint num_lj_pairs = top_eight->GetNumLJPairs();
  uint num_distance_constraints = top_eight->GetNumDistanceConstraints();
  uint num_exclusions = top_eight->GetNumExclusions();
  uint num_position_constraints = top_eight->GetNumPositionConstraints();
  uint num_harmonic_position_restraints = top_eight->GetNumHarmonicPositionRestraints();
  uint num_harmonic_distance_restraints = top_eight->GetNumHarmonicDistanceRestraints();

  temp = test_ent.Copy();
  top_eight->Merge(temp, merge_top_eight, new_ent);

  BOOST_CHECK(top_eight->GetNumHarmonicBonds() == num_harmonic_bonds+1);
  BOOST_CHECK(top_eight->GetNumHarmonicAngles() == num_harmonic_angles+1);
  BOOST_CHECK(top_eight->GetNumUreyBradleyAngles() == num_urey_bradley_angles+1);
  BOOST_CHECK(top_eight->GetNumPeriodicDihedrals() == num_periodic_dihedrals+1);
  BOOST_CHECK(top_eight->GetNumPeriodicImpropers() == num_periodic_impropers+1);
  BOOST_CHECK(top_eight->GetNumHarmonicImpropers() == num_harmonic_impropers+1);
  BOOST_CHECK(top_eight->GetNumCMaps() == num_cmaps+1);
  BOOST_CHECK(top_eight->GetNumLJPairs() == num_lj_pairs+1);
  BOOST_CHECK(top_eight->GetNumDistanceConstraints() == num_distance_constraints+1);
  BOOST_CHECK(top_eight->GetNumExclusions() == num_exclusions+1);
  BOOST_CHECK(top_eight->GetNumPositionConstraints() == num_position_constraints+1);
  BOOST_CHECK(top_eight->GetNumHarmonicPositionRestraints() == num_harmonic_position_restraints+1);
  BOOST_CHECK(top_eight->GetNumHarmonicDistanceRestraints() == num_harmonic_distance_restraints+1);

  uint atom_index_zero = 0 + test_ent.GetAtomCount();
  uint atom_index_one = 1 + test_ent.GetAtomCount();
  uint atom_index_two = 2 + test_ent.GetAtomCount();
  uint atom_index_three = 3 + test_ent.GetAtomCount();
  uint atom_index_four = 4 + test_ent.GetAtomCount();
  //uint atom_index_five = 5 + test_ent.GetAtomCount();

  //check whether the unique interactions are properly mapped
  BOOST_CHECK_THROW(top_eight->AddExclusion(atom_index_three,atom_index_four),ost::Error);
  BOOST_CHECK_THROW(top_eight->AddDistanceConstraint(atom_index_two,atom_index_three,10.0),ost::Error);
  BOOST_CHECK_THROW(top_eight->AddLJPair(atom_index_three,atom_index_four,1.0,2.0),ost::Error);

  //check whether the interactions are properly added (we do not check the correct parametrization at this point)

  std::vector<uint> harmonic_bond_indices = top_eight->GetHarmonicBondIndices(atom_index_zero,atom_index_one);
  std::vector<uint> harmonic_angle_indices = top_eight->GetHarmonicAngleIndices(atom_index_zero,atom_index_one,atom_index_two);
  std::vector<uint> urey_bradley_angle_indices = top_eight->GetUreyBradleyAngleIndices(atom_index_one,atom_index_two,atom_index_three);
  std::vector<uint> periodic_dihedral_indices = top_eight->GetPeriodicDihedralIndices(atom_index_zero,atom_index_one,atom_index_two,atom_index_three);
  std::vector<uint> periodic_improper_indices = top_eight->GetPeriodicImproperIndices(atom_index_one,atom_index_two,atom_index_three,atom_index_four);
  std::vector<uint> harmonic_improper_indices = top_eight->GetHarmonicImproperIndices(atom_index_zero,atom_index_one,atom_index_two,atom_index_three);
  std::vector<uint> cmap_indices = top_eight->GetCMapIndices(atom_index_zero,atom_index_one,atom_index_two,atom_index_three,atom_index_four);
  int lj_pair_index = top_eight->GetLJPairIndex(atom_index_three,atom_index_four);
  int distance_constraint_index = top_eight->GetDistanceConstraintIndex(atom_index_two,atom_index_three);
  std::vector<uint> harmonic_distance_restraint_indices = top_eight->GetHarmonicDistanceRestraintIndices(atom_index_one,atom_index_two);
  std::vector<uint> harmonic_position_restraint_indices = top_eight->GetHarmonicPositionRestraintIndices(atom_index_two);

  BOOST_CHECK(harmonic_bond_indices.size() == 1);
  BOOST_CHECK(harmonic_angle_indices.size() == 1);
  BOOST_CHECK(urey_bradley_angle_indices.size() == 1);
  BOOST_CHECK(periodic_dihedral_indices.size() == 1);
  BOOST_CHECK(periodic_improper_indices.size() == 1);
  BOOST_CHECK(harmonic_improper_indices.size() == 1);
  BOOST_CHECK(cmap_indices.size() == 1);
  BOOST_CHECK(harmonic_distance_restraint_indices.size() == 1);
  BOOST_CHECK(harmonic_position_restraint_indices.size() == 1);
  BOOST_CHECK(lj_pair_index != -1);
  BOOST_CHECK(distance_constraint_index != -1);


  //let's finally check the parameters

  uint ui0,ui1,ui2,ui3,ui4;
  int i0;
  Real r0,r1,r2,r3;
  geom::Vec3 v0;

  top_eight->GetHarmonicBondParameters(harmonic_bond_indices[0],ui0,ui1,r0,r1);
  BOOST_CHECK(r0 == 42.0 && r1 == 42000.0);

  top_eight->GetHarmonicAngleParameters(harmonic_angle_indices[0],ui0,ui1,ui2,r0,r1);
  BOOST_CHECK(r0 == 24.0 && r1 == 24000.0);

  top_eight->GetUreyBradleyAngleParameters(urey_bradley_angle_indices[0],ui0,ui1,ui2,r0,r1,r2,r3);
  BOOST_CHECK(r0 == 1.0 && r1 == 2.0 && r2 == 3.0 && r3 == 4.0);

  top_eight->GetPeriodicDihedralParameters(periodic_dihedral_indices[0],ui0,ui1,ui2,ui3,i0,r0,r1);
  BOOST_CHECK(i0 == 10 && r0 == 2.0 && r1 == 42.0);

  top_eight->GetPeriodicImproperParameters(periodic_improper_indices[0],ui0,ui1,ui2,ui3,i0,r0,r1);
  BOOST_CHECK(i0 == 3 && r0 == 5.0 && r1 == 10.0);

  top_eight->GetHarmonicImproperParameters(harmonic_improper_indices[0],ui0,ui1,ui2,ui3,r0,r1);
  BOOST_CHECK(r0 == 8.0 && r1 == 12.0);

  top_eight->GetCMapParameters(cmap_indices[0],ui0,ui1,ui2,ui3,ui4,i0,real_vec);
  BOOST_CHECK(real_vec == cmap && i0 == 2);

  top_eight->GetLJPairParameters(lj_pair_index,ui0,ui1,r0,r1);
  BOOST_CHECK(r0 == 1.0 && r1 == 2.0);

  top_eight->GetDistanceConstraintParameters(distance_constraint_index,ui0,ui1,r0);
  BOOST_CHECK(r0 == 10.0);

  top_eight->GetHarmonicPositionRestraintParameters(harmonic_position_restraint_indices[0],ui0,v0,r0,r1,r2,r3);
  BOOST_CHECK(v0 == pos && r0 == 10.0 && r1 == 5.0 && r2 == 6.0 && r3 == 7.0);

  top_eight->GetHarmonicDistanceRestraintParameters(harmonic_distance_restraint_indices[0],ui0,ui1,r0,r1);
  BOOST_CHECK(r0 == 10.0 && r1 == 1000.0);

}


BOOST_AUTO_TEST_SUITE_END();