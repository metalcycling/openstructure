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
#include <ost/mol/mm/interaction.hh>
#include <ost/mol/mm/buildingblock.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/builder.hh>
#include <ost/mol/xcs_editor.hh>


#include <ost/message.hh>

using namespace ost::mol::mm;

BOOST_AUTO_TEST_SUITE( mol_mm );

BOOST_AUTO_TEST_CASE(test_block_basics) 
{
  BuildingBlock block = BuildingBlock();
  block.AddAtom("A","a",0.5);
  block.AddAtom("B","b",1.0);
  InteractionPtr bond(new Interaction(HarmonicBond));
  std::vector<String> names;
  names.push_back("A");
  names.push_back("B");
  bond->SetNames(names);
  block.AddBond(bond);
  BOOST_CHECK_THROW(block.GetCharge("C"),ost::Error);
  BOOST_CHECK_THROW(block.GetType("C"),ost::Error);
  BOOST_CHECK(block.GetCharge("A") == 0.5);
  BOOST_CHECK(block.GetType("A") == "a");
  block.ReplaceAtom("A","X","x",42.0);
  BOOST_CHECK(block.GetType("X") == "x");
  BOOST_CHECK(block.GetCharge("X")==42.0);
  BOOST_CHECK(block.GetAtoms().size() == 2);
  std::vector<InteractionPtr> bonds = block.GetBonds();
  BOOST_CHECK(bonds[0]->HasName("X"));
  block.RemoveAtom("X");
  BOOST_CHECK_THROW(block.GetType("X"),ost::Error);
  BOOST_CHECK_NO_THROW(block.GetType("B"));
  BOOST_CHECK(block.GetAtoms().size()==1);
  BOOST_CHECK(block.GetCharges().size()==1);
  BOOST_CHECK(block.GetBonds().size()==0);
  std::vector<Real> param_one,param_two;
  param_one.push_back(1.0);
  param_one.push_back(1.0);
  param_two.push_back(2.0);
  param_two.push_back(2.0);
  InteractionPtr bond_one(new Interaction(HarmonicBond));
  InteractionPtr bond_two(new Interaction(HarmonicBond));
  bond_one->SetNames(names);
  bond_one->SetParam(param_one);
  bond_two->SetNames(names);
  bond_two->SetParam(param_two);
  block.AddBond(bond_one);
  BOOST_CHECK(block.GetBonds()[0]->GetParam() == param_one);
  block.AddBond(bond_two,true);
  BOOST_CHECK(block.GetBonds()[0]->GetParam() == param_two);
  block.RemoveAtom("A");

  names[0] = "-A";
  bond_one->SetNames(names);
  names[0] = "A";
  names[1] = "+B";
  bond_two->SetNames(names);
  block.AddBond(bond_one);
  block.AddBond(bond_two);
  BOOST_CHECK(block.GetBonds().size() == 2);
  block.RemoveInteractionsToPrev();
  BOOST_CHECK(block.GetBonds().size() == 1);
  BOOST_CHECK(block.GetBonds()[0]->HasName("+B"));
  block.RemoveInteractionsToNext();
  BOOST_CHECK(block.GetBonds().size() == 0);
}

BOOST_AUTO_TEST_CASE(test_connect)
{
  ost::mol::EntityHandle e = ost::mol::Builder()
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

  ost::mol::XCSEditor ed = e.EditXCS();
  ed.DeleteBonds(e.GetBondList());

  BuildingBlock block = BuildingBlock();
  block.AddAtom("A","a",0.0);
  block.AddAtom("B","b",0.0);

  InteractionPtr bond_one(new Interaction(HarmonicBond));
  InteractionPtr bond_two(new Interaction(HarmonicBond));
  InteractionPtr bond_three(new Interaction(HarmonicBond));
  std::vector<String> two_string;
  two_string.push_back("A");
  two_string.push_back("B");
  bond_one->SetNames(two_string);
  two_string[0] = "-B";
  two_string[1] = "A";
  bond_two->SetNames(two_string);
  two_string[0] = "B";
  two_string[1] = "+A";
  bond_three->SetNames(two_string);
  block.AddBond(bond_one);
  block.AddBond(bond_two);
  block.AddBond(bond_three);

  ost::mol::ResidueHandleList res_list = e.GetResidueList();
  ost::mol::AtomHandleList atoms_one = res_list[0].GetAtomList();
  ost::mol::AtomHandleList atoms_two = res_list[1].GetAtomList();
  ost::mol::AtomHandleList atoms_three = res_list[2].GetAtomList();

  block.Connect(res_list[1],ed);
  BOOST_CHECK(ost::mol::BondExists(atoms_one[1],atoms_two[0]));
  BOOST_CHECK(ost::mol::BondExists(atoms_two[0],atoms_two[1]));
  BOOST_CHECK(ost::mol::BondExists(atoms_two[1],atoms_three[0]));

  ed.DeleteBonds(e.GetBondList());

  BOOST_CHECK_THROW(block.Connect(res_list[0],ed),ost::Error);
  BOOST_CHECK_THROW(block.Connect(res_list[2],ed),ost::Error);

  InteractionPtr bond_four(new Interaction(HarmonicBond));
  two_string[1] = "X";
  bond_four->SetNames(two_string);
  block.AddBond(bond_four);

  BOOST_CHECK_THROW(block.Connect(res_list[1],ed),ost::Error);
}

BOOST_AUTO_TEST_CASE(test_match)
{

  ost::mol::EntityHandle e = ost::mol::Builder()
                             .Chain("A")
                                .Residue("ONE")
                                  .Atom("A",geom::Vec3(0,0,0))
                                  .Atom("B",geom::Vec3(1,1,1))
                                .Residue("TWO")
                                  .Atom("A",geom::Vec3(2,2,2))
                                  .Atom("B",geom::Vec3(3,3,3))
                                .Residue("THREE")
                                  .Atom("A",geom::Vec3(4,4,4))
                                  .Atom("B",geom::Vec3(5,5,5))
                                .Residue("FOUR")
                                  .Atom("A",geom::Vec3(6,6,6))
                                  .Atom("B",geom::Vec3(7,7,7));

  ost::mol::ResidueHandleList res_list = e.GetResidueList();
  ost::mol::XCSEditor ed = e.EditXCS();

  ed.Connect(res_list[0].FindAtom("A"),res_list[0].FindAtom("B"));
  ed.Connect(res_list[0].FindAtom("B"),res_list[1].FindAtom("A"));
  ed.Connect(res_list[1].FindAtom("A"),res_list[1].FindAtom("B"));
  ed.Connect(res_list[1].FindAtom("B"),res_list[2].FindAtom("A"));
  ed.Connect(res_list[2].FindAtom("A"),res_list[2].FindAtom("B"));
  //ed.Connect(atoms_two[1],atoms_four[1]);

  BuildingBlock block = BuildingBlock();
  String info_string;
  BOOST_CHECK(!block.Match(res_list[1],true,info_string));

  block.AddAtom("A","a",0.0);
  block.AddAtom("B","b",0.0);

  BOOST_CHECK(!block.Match(res_list[1],true,info_string));

  InteractionPtr bond_one(new Interaction(HarmonicBond));
  InteractionPtr bond_two(new Interaction(HarmonicBond));
  InteractionPtr bond_three(new Interaction(HarmonicBond));

  std::vector<String> names_one;
  names_one.push_back("A");
  names_one.push_back("B");
  std::vector<String> names_two;
  names_two.push_back("-B");
  names_two.push_back("A");  
  std::vector<String> names_three;
  names_three.push_back("B");
  names_three.push_back("+A");
  
  bond_one->SetNames(names_one);
  bond_two->SetNames(names_two);
  bond_three->SetNames(names_three);

  block.AddBond(bond_one);
  BOOST_CHECK(block.Match(res_list[1],true,info_string));

  block.AddBond(bond_two);
  BOOST_CHECK(block.Match(res_list[1],true,info_string));

  block.AddBond(bond_three);
  BOOST_CHECK(block.Match(res_list[1],true,info_string));

  //let's remove the central bond and replace it with a constraint

  block.RemoveAtom("A");
  block.AddAtom("A","a",0.0);
  InteractionPtr constraint(new Interaction(DistanceConstraint));
  constraint->SetNames(names_one);
  BOOST_CHECK(!block.Match(res_list[1],true,info_string));
  block.AddConstraint(constraint);
  BOOST_CHECK(block.Match(res_list[1],true,info_string));
}

BOOST_AUTO_TEST_SUITE_END();
