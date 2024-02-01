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


using namespace ost::mol::mm;

BOOST_AUTO_TEST_SUITE( mol_mm );

BOOST_AUTO_TEST_CASE(test_forcefield_basics){
  Forcefield ff;

  //if a building block doesn't exist, it should throw an error
  //this is not the case for hydrogen constructors and block modifiers
  BOOST_CHECK_THROW(ff.GetBuildingBlock("IDoNotExist"),ost::Error);

  //add some atom types
  ff.AddMass("one",1.0);
  ff.AddMass("two",2.0);
  ff.AddMass("three",3.0);
  ff.AddMass("four",4.0);
  ff.AddMass("five",5.0);

  BOOST_CHECK_THROW(ff.GetBond("one","two"),ost::Error);
  BOOST_CHECK_THROW(ff.GetAngle("one","two","three"),ost::Error);
  BOOST_CHECK_THROW(ff.GetDihedrals("one","two","three","four"),ost::Error);
  BOOST_CHECK_THROW(ff.GetImpropers("one","two","three","four"),ost::Error);
  BOOST_CHECK_THROW(ff.GetCMap("one","two","three","four","five"),ost::Error);
  BOOST_CHECK_THROW(ff.GetImplicitGenborn("one"),ost::Error);
  BOOST_CHECK_THROW(ff.GetLJ("one"),ost::Error);
  BOOST_CHECK_THROW(ff.GetLJ("one","two"),ost::Error);
  BOOST_CHECK_THROW(ff.GetConstraint("one","two"),ost::Error);
  BOOST_CHECK_THROW(ff.GetMass("i_do_not_exist"),ost::Error);

  //let's add some interactions
  std::vector<String> one_string, two_string, three_string, four_string, five_string;
  std::vector<Real> one_real, two_real, three_real, four_real, five_real;

  one_string.push_back("one");
  two_string.push_back("one");
  three_string.push_back("one");
  four_string.push_back("one");
  five_string.push_back("one");

  two_string.push_back("two");
  three_string.push_back("two");
  four_string.push_back("two");
  five_string.push_back("two");

  three_string.push_back("three");
  four_string.push_back("three");
  five_string.push_back("three");

  four_string.push_back("four");
  five_string.push_back("four");

  five_string.push_back("five");


  one_real.push_back(1.0);
  two_real.push_back(1.0);
  three_real.push_back(1.0);
  four_real.push_back(1.0);
  five_real.push_back(1.0);

  two_real.push_back(2.0);
  three_real.push_back(2.0);
  four_real.push_back(2.0);
  five_real.push_back(2.0);

  three_real.push_back(3.0);
  four_real.push_back(3.0);
  five_real.push_back(3.0);

  four_real.push_back(4.0);
  five_real.push_back(4.0);

  five_real.push_back(5.0);


  InteractionPtr harmonic_bond(new Interaction(HarmonicBond));
  InteractionPtr urey_bradley_angle(new Interaction(UreyBradleyAngle));
  InteractionPtr harmonic_angle(new Interaction(HarmonicAngle));
  InteractionPtr periodic_dihedral(new Interaction(PeriodicDihedral));
  InteractionPtr periodic_improper(new Interaction(PeriodicImproper));
  InteractionPtr harmonic_improper(new Interaction(HarmonicImproper));
  InteractionPtr cmap(new Interaction(CMap));
  InteractionPtr lj(new Interaction(LJ));
  InteractionPtr lj_pair(new Interaction(LJPair));
  InteractionPtr gbsa(new Interaction(GBSA));
  InteractionPtr distance_constraint(new Interaction(DistanceConstraint));


  //check whether errors get thrown when incorrectly parametrized interactions are added
  BOOST_CHECK_THROW(ff.AddBond(harmonic_bond),ost::Error);
  harmonic_bond->SetTypes(two_string);
  BOOST_CHECK_THROW(ff.AddBond(harmonic_bond),ost::Error);
  harmonic_angle->SetParam(two_real);
  BOOST_CHECK_THROW(ff.AddAngle(harmonic_angle),ost::Error);
  harmonic_bond->SetParam(two_real);
  BOOST_CHECK_NO_THROW(ff.AddBond(harmonic_bond));
  harmonic_angle->SetTypes(three_string);

  //parametrize the rest
  urey_bradley_angle->SetTypes(three_string);
  urey_bradley_angle->SetParam(four_real);
  periodic_dihedral->SetTypes(four_string);
  periodic_dihedral->SetParam(three_real);
  periodic_improper->SetTypes(four_string);
  periodic_improper->SetParam(three_real);
  harmonic_improper->SetTypes(four_string);
  harmonic_improper->SetParam(two_real);
  std::vector<Real> cmap_param = five_real;
  cmap_param[0] = 2.0;
  cmap->SetTypes(five_string);
  cmap->SetParam(cmap_param);
  lj->SetTypes(one_string);
  lj->SetParam(two_real);
  lj_pair->SetTypes(two_string);
  lj_pair->SetParam(two_real);
  gbsa->SetTypes(one_string);
  gbsa->SetParam(two_real);
  distance_constraint->SetTypes(two_string);
  distance_constraint->SetParam(one_real);

  //check whether errors get thrown when wrong interactions get added

  BOOST_CHECK_THROW(ff.AddBond(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(ff.AddAngle(harmonic_bond),ost::Error);
  BOOST_CHECK_THROW(ff.AddDihedral(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(ff.AddImproper(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(ff.AddCMap(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(ff.AddImplicitGenborn(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(ff.AddLJ(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(ff.AddLJPair(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(ff.AddConstraint(harmonic_angle),ost::Error);

  BOOST_CHECK_NO_THROW(ff.AddBond(harmonic_bond));
  BOOST_CHECK_NO_THROW(ff.AddAngle(harmonic_angle));
  BOOST_CHECK_NO_THROW(ff.AddDihedral(periodic_dihedral));
  BOOST_CHECK_NO_THROW(ff.AddImproper(periodic_improper));
  BOOST_CHECK_NO_THROW(ff.AddImproper(harmonic_improper));
  BOOST_CHECK_NO_THROW(ff.AddCMap(cmap));
  BOOST_CHECK_NO_THROW(ff.AddImplicitGenborn(gbsa));
  BOOST_CHECK_NO_THROW(ff.AddLJ(lj));
  BOOST_CHECK_NO_THROW(ff.AddLJPair(lj_pair));
  BOOST_CHECK_NO_THROW(ff.AddConstraint(distance_constraint));


  //add a building block
  BuildingBlockPtr block(new BuildingBlock);
  block->AddAtom("a","one",0.5);
  block->AddAtom("b","two",1.0);
  ff.AddBuildingBlock("im_so_sexy",block);

  //add block modifier
  BlockModifierPtr modifier(new GromacsBlockModifier);
  ff.AddBlockModifier("el_modificator",modifier);

  //add hydrogen constructor
  HydrogenConstructorPtr hydrogen_constructor(new GromacsHydrogenConstructor);
  ff.AddHydrogenConstructor("death_metal",hydrogen_constructor);

  //add some renaming rules
  ff.AddResidueRenamingRule("im_so_sexy","even_sexier","sexy_n","sexy_c","sexy_two");
  ff.AddAtomRenamingRule("even_sexier","a","z");
  ff.AddAtomRenamingRule("sexy_two","a","z");

  //save and load it again, before we check parameters
  ff.Save("test.ff");
  ForcefieldPtr loaded_ff = Forcefield::Load("test.ff");


  //check the parameters
  BOOST_CHECK(loaded_ff->GetBond("one","two")->GetParam() == two_real);
  BOOST_CHECK(loaded_ff->GetBond("two","one")->GetParam() == two_real);

  BOOST_CHECK(loaded_ff->GetAngle("one","two","three")->GetParam() == two_real);
  BOOST_CHECK(loaded_ff->GetAngle("three","two","one")->GetParam() == two_real);

  BOOST_CHECK(loaded_ff->GetDihedrals("one","two","three","four")[0]->GetParam() == three_real);
  BOOST_CHECK(loaded_ff->GetDihedrals("four","three","two","one")[0]->GetParam() == three_real);

  BOOST_CHECK(loaded_ff->GetImpropers("one","two","three","four")[0]->GetParam() == three_real);
  BOOST_CHECK(loaded_ff->GetImpropers("four","three","two","one")[0]->GetParam() == three_real);

  BOOST_CHECK(loaded_ff->GetDihedrals("one","two","three","four")[0]->GetParam() == three_real);
  BOOST_CHECK(loaded_ff->GetDihedrals("four","three","two","one")[0]->GetParam() == three_real);

  BOOST_CHECK(loaded_ff->GetCMap("one","two","three","four","five")->GetParam() == cmap_param);
  //cmaps are directional...
  BOOST_CHECK_THROW(loaded_ff->GetCMap("five","four","three","two","one"),ost::Error);

  BOOST_CHECK(loaded_ff->GetImplicitGenborn("one")->GetParam() == two_real);

  BOOST_CHECK(loaded_ff->GetLJ("one","two",true)->GetParam() == two_real);
  BOOST_CHECK(loaded_ff->GetLJ("two","one",true)->GetParam() == two_real);

  BOOST_CHECK(loaded_ff->GetLJ("one")->GetParam() == two_real);

  BOOST_CHECK(loaded_ff->GetConstraint("one","two")->GetParam() == one_real);
  BOOST_CHECK(loaded_ff->GetConstraint("two","one")->GetParam() == one_real);

  BOOST_CHECK(loaded_ff->GetMass("two") == 2.0);

  BOOST_CHECK(loaded_ff->GetAtomType("im_so_sexy","a") == "one");

  //check whether the modifiers are still there
  BOOST_CHECK(loaded_ff->GetHydrogenConstructor("death_metal"));
  BOOST_CHECK(loaded_ff->GetBlockModifier("el_modificator"));


  //check the renaming stuff


  ost::mol::EntityHandle new_ent = ost::mol::Builder()
                                   .Chain("A")
                                     .Residue("im_so_sexy")
                                       .Atom("a",geom::Vec3(0,0,0))
                                       .Atom("b",geom::Vec3(1,1,1))
                                     .Residue("im_so_sexy")
                                       .Atom("a",geom::Vec3(2,2,2))
                                       .Atom("b",geom::Vec3(3,3,3))
                                     .Residue("im_so_sexy")
                                       .Atom("a",geom::Vec3(4,4,4))
                                       .Atom("b",geom::Vec3(5,5,5));

  ost::mol::XCSEditor ed = new_ent.EditXCS();


  loaded_ff->AssignFFSpecificNames(new_ent);
  ost::mol::ResidueHandleList res_list = new_ent.GetResidueList();
  ost::mol::AtomHandleList atom_list = new_ent.GetAtomList();

  BOOST_CHECK(res_list[0].GetName() == "sexy_n");
  BOOST_CHECK(res_list[1].GetName() == "even_sexier");
  BOOST_CHECK(res_list[2].GetName() == "sexy_c");
  BOOST_CHECK(atom_list[2].GetName() == "z");

  loaded_ff->AssignFFSpecificNames(new_ent,true);

  BOOST_CHECK(res_list[0].GetName() == "im_so_sexy");
  BOOST_CHECK(res_list[1].GetName() == "im_so_sexy");
  BOOST_CHECK(res_list[2].GetName() == "im_so_sexy");
  BOOST_CHECK(atom_list[2].GetName() == "a");

  ost::mol::EntityHandle new_ent_two = ost::mol::Builder()
                                      .Chain("A")
                                        .Residue("im_so_sexy")
                                          .Atom("a",geom::Vec3(0,0,0))
                                          .Atom("b",geom::Vec3(1,1,1));

  loaded_ff->AssignFFSpecificNames(new_ent_two);
  res_list = new_ent_two.GetResidueList();
  atom_list = new_ent_two.GetAtomList();
  BOOST_CHECK(res_list[0].GetName() == "sexy_two");
  BOOST_CHECK(atom_list[0].GetName() == "z");
  loaded_ff->AssignFFSpecificNames(new_ent_two,true);
  BOOST_CHECK(res_list[0].GetName() == "im_so_sexy");
  BOOST_CHECK(atom_list[0].GetName() == "a");

}

BOOST_AUTO_TEST_SUITE_END();