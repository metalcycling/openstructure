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
#include <ost/mol/mm/gromacs_block_modifiers.hh>


#include <ost/message.hh>

using namespace ost::mol::mm;

BOOST_AUTO_TEST_SUITE( mol_mm );

BOOST_AUTO_TEST_CASE(test_gromacs_block_modifier_basics){

  GromacsBlockModifier mod;

  //let's create some interactions
  InteractionPtr harmonic_bond(new Interaction(HarmonicBond));
  InteractionPtr urey_bradley_angle(new Interaction(UreyBradleyAngle));
  InteractionPtr harmonic_angle(new Interaction(HarmonicAngle));
  InteractionPtr periodic_dihedral(new Interaction(PeriodicDihedral));
  InteractionPtr periodic_improper(new Interaction(PeriodicImproper));
  InteractionPtr harmonic_improper(new Interaction(HarmonicImproper));
  InteractionPtr cmap(new Interaction(CMap));

  //check whether error gets thrown, when names of interactions are not set
  BOOST_CHECK_THROW(mod.AddBond(harmonic_bond),ost::Error);
  BOOST_CHECK_THROW(mod.AddAngle(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(mod.AddDihedral(periodic_dihedral),ost::Error);
  BOOST_CHECK_THROW(mod.AddImproper(periodic_improper),ost::Error);
  BOOST_CHECK_THROW(mod.AddCMap(cmap),ost::Error);

  std::vector<String> two_string, three_string, four_string, five_string;

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

  harmonic_bond->SetNames(two_string);
  harmonic_angle->SetNames(three_string);
  urey_bradley_angle->SetNames(three_string);
  periodic_dihedral->SetNames(four_string);
  periodic_improper->SetNames(four_string);
  harmonic_improper->SetNames(four_string);
  cmap->SetNames(five_string);

  BOOST_CHECK_NO_THROW(mod.AddBond(harmonic_bond));
  BOOST_CHECK_NO_THROW(mod.AddAngle(harmonic_angle));
  BOOST_CHECK_NO_THROW(mod.AddAngle(urey_bradley_angle));
  BOOST_CHECK_NO_THROW(mod.AddDihedral(periodic_dihedral));
  BOOST_CHECK_NO_THROW(mod.AddImproper(periodic_improper));
  BOOST_CHECK_NO_THROW(mod.AddImproper(harmonic_improper));
  BOOST_CHECK_NO_THROW(mod.AddCMap(cmap));

  //check whether error get thrown when wrong interactions get added

  BOOST_CHECK_THROW(mod.AddBond(harmonic_angle),ost::Error);
  BOOST_CHECK_THROW(mod.AddAngle(harmonic_bond),ost::Error);
  BOOST_CHECK_THROW(mod.AddDihedral(harmonic_bond),ost::Error);
  BOOST_CHECK_THROW(mod.AddImproper(periodic_dihedral),ost::Error);
  BOOST_CHECK_THROW(mod.AddCMap(periodic_improper),ost::Error);
  
}

BOOST_AUTO_TEST_SUITE_END();
