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
#include <ost/mol/mol.hh>
#include <ost/mol/builder.hh>

#include <ost/message.hh>

using namespace ost::mol::mm;



BOOST_AUTO_TEST_SUITE( mol_mm );

BOOST_AUTO_TEST_CASE(test_set_stuff) 
{
  Interaction harmonic_bond = Interaction(HarmonicBond);
  Interaction harmonic_angle = Interaction(HarmonicAngle);
  Interaction urey_bradley = Interaction(UreyBradleyAngle);
  Interaction periodic_dihedral = Interaction(PeriodicDihedral);
  Interaction periodic_improper = Interaction(PeriodicImproper);
  Interaction harmonic_improper = Interaction(HarmonicImproper);
  Interaction cmap = Interaction(CMap);
  Interaction lj = Interaction(LJ);
  Interaction lj_pair = Interaction(LJPair);
  Interaction gbsa = Interaction(GBSA);
  Interaction distance_constraint = Interaction(DistanceConstraint);
  Interaction exclusion = Interaction(Exclusion);
  Interaction harmonic_distance_restraint = Interaction(HarmonicDistanceRestraint);
  Interaction harmonic_position_restraint = Interaction(HarmonicPositionRestraint);

  //std::vector<String> one_string, two_string, three_string, four_string, five_string;
  //std::vector<Real> one_real, two_real, three_real, four_real, five_real, zero_real;

  std::vector<Real> one_real(1);
  std::vector<Real> two_real(2);
  std::vector<Real> three_real(3);
  std::vector<Real> four_real(4);
  std::vector<Real> five_real(5);
  std::vector<Real> seven_real(7);
  std::vector<Real> zero_real(0);

  std::vector<String> one_string(1);
  std::vector<String> two_string(2);
  std::vector<String> three_string(3);
  std::vector<String> four_string(4);
  std::vector<String> five_string(5);


  BOOST_CHECK(!harmonic_bond.IsParametrized());

  BOOST_CHECK_THROW(harmonic_bond.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(harmonic_angle.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(urey_bradley.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(periodic_dihedral.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(periodic_improper.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(harmonic_improper.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(cmap.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(lj.SetNames(two_string),ost::Error);
  BOOST_CHECK_THROW(lj_pair.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(gbsa.SetNames(two_string),ost::Error);
  BOOST_CHECK_THROW(distance_constraint.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(exclusion.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(harmonic_distance_restraint.SetNames(one_string),ost::Error);
  BOOST_CHECK_THROW(harmonic_position_restraint.SetNames(two_string),ost::Error);


  BOOST_CHECK_THROW(harmonic_bond.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(harmonic_angle.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(urey_bradley.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(periodic_dihedral.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(periodic_improper.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(harmonic_improper.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(cmap.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(lj.SetTypes(two_string),ost::Error);
  BOOST_CHECK_THROW(lj_pair.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(gbsa.SetTypes(two_string),ost::Error);
  BOOST_CHECK_THROW(distance_constraint.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(exclusion.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(harmonic_distance_restraint.SetTypes(one_string),ost::Error);
  BOOST_CHECK_THROW(harmonic_position_restraint.SetTypes(two_string),ost::Error);

  BOOST_CHECK_THROW(harmonic_bond.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(harmonic_angle.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(urey_bradley.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(periodic_dihedral.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(periodic_improper.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(harmonic_improper.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(cmap.SetParam(two_real),ost::Error);
  BOOST_CHECK_THROW(lj.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(lj_pair.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(gbsa.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(distance_constraint.SetParam(two_real),ost::Error);
  BOOST_CHECK_THROW(exclusion.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(harmonic_distance_restraint.SetParam(one_real),ost::Error);
  BOOST_CHECK_THROW(harmonic_position_restraint.SetParam(one_real),ost::Error);



  BOOST_CHECK_NO_THROW(harmonic_bond.SetNames(two_string));
  BOOST_CHECK_NO_THROW(harmonic_angle.SetNames(three_string));
  BOOST_CHECK_NO_THROW(urey_bradley.SetNames(three_string));
  BOOST_CHECK_NO_THROW(periodic_dihedral.SetNames(four_string));
  BOOST_CHECK_NO_THROW(periodic_improper.SetNames(four_string));
  BOOST_CHECK_NO_THROW(harmonic_improper.SetNames(four_string));
  BOOST_CHECK_NO_THROW(cmap.SetNames(five_string));
  BOOST_CHECK_NO_THROW(lj.SetNames(one_string));
  BOOST_CHECK_NO_THROW(lj_pair.SetNames(two_string));
  BOOST_CHECK_NO_THROW(gbsa.SetNames(one_string));
  BOOST_CHECK_NO_THROW(distance_constraint.SetNames(two_string));
  BOOST_CHECK_NO_THROW(exclusion.SetNames(two_string));
  BOOST_CHECK_NO_THROW(harmonic_distance_restraint.SetNames(two_string));
  BOOST_CHECK_NO_THROW(harmonic_position_restraint.SetNames(one_string));
  

  BOOST_CHECK_NO_THROW(harmonic_bond.SetTypes(two_string));
  BOOST_CHECK_NO_THROW(harmonic_angle.SetTypes(three_string));
  BOOST_CHECK_NO_THROW(urey_bradley.SetTypes(three_string));
  BOOST_CHECK_NO_THROW(periodic_dihedral.SetTypes(four_string));
  BOOST_CHECK_NO_THROW(periodic_improper.SetTypes(four_string));
  BOOST_CHECK_NO_THROW(harmonic_improper.SetTypes(four_string));
  BOOST_CHECK_NO_THROW(cmap.SetTypes(five_string));
  BOOST_CHECK_NO_THROW(lj.SetTypes(one_string));
  BOOST_CHECK_NO_THROW(lj_pair.SetTypes(two_string));
  BOOST_CHECK_NO_THROW(gbsa.SetTypes(one_string));
  BOOST_CHECK_NO_THROW(distance_constraint.SetTypes(two_string));
  BOOST_CHECK_NO_THROW(exclusion.SetTypes(two_string));
  BOOST_CHECK_NO_THROW(harmonic_distance_restraint.SetTypes(two_string));
  BOOST_CHECK_NO_THROW(harmonic_position_restraint.SetTypes(one_string));

  BOOST_CHECK_NO_THROW(harmonic_bond.SetParam(two_real));
  BOOST_CHECK_NO_THROW(harmonic_angle.SetParam(two_real));
  BOOST_CHECK_NO_THROW(urey_bradley.SetParam(four_real));
  BOOST_CHECK_NO_THROW(periodic_dihedral.SetParam(three_real));
  BOOST_CHECK_NO_THROW(periodic_improper.SetParam(three_real));
  BOOST_CHECK_NO_THROW(harmonic_improper.SetParam(two_real));
  BOOST_CHECK_NO_THROW(lj.SetParam(two_real));
  BOOST_CHECK_NO_THROW(lj_pair.SetParam(two_real));
  BOOST_CHECK_NO_THROW(gbsa.SetParam(two_real));
  BOOST_CHECK_NO_THROW(distance_constraint.SetParam(one_real));
  BOOST_CHECK_NO_THROW(exclusion.SetParam(zero_real));
  BOOST_CHECK_NO_THROW(harmonic_distance_restraint.SetParam(two_real));
  BOOST_CHECK_NO_THROW(harmonic_position_restraint.SetParam(seven_real));

  std::vector<Real> valid_cmap_parameters;
  valid_cmap_parameters.push_back(2.0);
  valid_cmap_parameters.push_back(0.0);
  valid_cmap_parameters.push_back(0.0);
  valid_cmap_parameters.push_back(0.0);
  valid_cmap_parameters.push_back(0.0);

  BOOST_CHECK_NO_THROW(cmap.SetParam(valid_cmap_parameters));
  BOOST_CHECK(harmonic_bond.IsParametrized());

  Interaction new_bond = Interaction(HarmonicBond);
  std::vector<String> with_wildcard, without_wildcard;
  with_wildcard.push_back("A");
  without_wildcard.push_back("A");
  with_wildcard.push_back("X");
  without_wildcard.push_back("A");

  //should be initialized with false
  BOOST_CHECK(!new_bond.HasNameWildcard());
  BOOST_CHECK(!new_bond.HasTypeWildcard());

  new_bond.SetNames(without_wildcard);
  new_bond.SetTypes(without_wildcard);
  BOOST_CHECK(!new_bond.HasNameWildcard());
  BOOST_CHECK(!new_bond.HasTypeWildcard());

  new_bond.SetNames(with_wildcard);
  new_bond.SetTypes(with_wildcard);
  BOOST_CHECK(new_bond.HasNameWildcard());
  BOOST_CHECK(new_bond.HasTypeWildcard());

  //Check hasname/hastype
  BOOST_CHECK(new_bond.HasType("A"));
  BOOST_CHECK(!new_bond.HasType("B"));
  BOOST_CHECK(new_bond.HasName("A"));
  BOOST_CHECK(!new_bond.HasType("B"));

}

BOOST_AUTO_TEST_CASE(test_match_stuff) 
{
  Interaction dihedral = Interaction(PeriodicDihedral);

  std::vector<String> with_wildcard, without_wildcard;
  with_wildcard.push_back("X");
  with_wildcard.push_back("A");
  with_wildcard.push_back("A");
  with_wildcard.push_back("A");
  without_wildcard.push_back("A");
  without_wildcard.push_back("A");
  without_wildcard.push_back("A");
  without_wildcard.push_back("A");

  dihedral.SetTypes(with_wildcard);
  BOOST_CHECK(dihedral.MatchTypes(with_wildcard));
  BOOST_CHECK(dihedral.MatchTypes(without_wildcard));
  dihedral.SetTypes(without_wildcard);
  BOOST_CHECK(dihedral.MatchTypes(without_wildcard));
  BOOST_CHECK(!dihedral.MatchTypes(with_wildcard));

  dihedral.SetNames(with_wildcard);
  BOOST_CHECK(dihedral.MatchNames(with_wildcard));
  BOOST_CHECK(dihedral.MatchNames(without_wildcard));
  dihedral.SetNames(without_wildcard);
  BOOST_CHECK(dihedral.MatchNames(without_wildcard));
  BOOST_CHECK(!dihedral.MatchNames(with_wildcard));

}

BOOST_AUTO_TEST_CASE(test_replace_atom){
  Interaction bond = Interaction(HarmonicBond);
  std::vector<String> strings;
  strings.push_back("A");
  strings.push_back("B");
  bond.SetNames(strings);
  bond.ReplaceAtom("A","X","X");
  std::vector<String> new_names = bond.GetNames();
  BOOST_CHECK(new_names != strings);
  bond.SetTypes(strings);
  bond.ReplaceAtom("X","A","X");
  std::vector<String> new_types = bond.GetTypes();
  BOOST_CHECK(new_types == new_names);

}

BOOST_AUTO_TEST_CASE(test_getatoms){
  ost::mol::EntityHandle e = ost::mol::Builder()
                             .Chain("A")
                                .Residue("ONE")
                                  .Atom("A",geom::Vec3(-5,-5,-5))
                                  .Atom("B",geom::Vec3(-5, 5,-5))
                                .Residue("TWO")
                                  .Atom("A",geom::Vec3(1,1,1))
                                  .Atom("B",geom::Vec3(2,2,2))
                                .Residue("THREE")
                                  .Atom("A",geom::Vec3(3,3,3))
                                  .Atom("B",geom::Vec3(4,4,4));

  ost::mol::ResidueHandleList res_list = e.GetResidueList();
  ost::mol::AtomHandleList atom_list_one  = res_list[0].GetAtomList();
  ost::mol::AtomHandleList atom_list_two  = res_list[1].GetAtomList();
  ost::mol::AtomHandleList atom_list_three  = res_list[2].GetAtomList();

  Interaction bond = Interaction(HarmonicBond);
  std::vector<String> strings;
  strings.push_back("A");
  strings.push_back("B");
  bond.SetNames(strings);
  ost::mol::AtomHandleList bond_atom_list = bond.GetAtoms(res_list[0]);
  BOOST_CHECK(atom_list_one[0] == bond_atom_list[0]);
  BOOST_CHECK(atom_list_one[1] == bond_atom_list[1]);
  strings[1] = "X";
  bond.SetNames(strings);
  BOOST_CHECK_THROW(bond.GetAtoms(res_list[0]),ost::Error);
  strings[0] = "-A";
  strings[1] = "B";
  bond.SetNames(strings);
  BOOST_CHECK_THROW(bond.GetAtoms(res_list[0]),ost::Error);
  bond_atom_list = bond.GetAtoms(res_list[1]);  
  BOOST_CHECK(bond_atom_list[0] == atom_list_one[0]);
  BOOST_CHECK(bond_atom_list[1] == atom_list_two[1]);
  strings[1]="+B";
  bond.SetNames(strings);
  BOOST_CHECK_THROW(bond.GetAtoms(res_list[2]),ost::Error);
  bond_atom_list = bond.GetAtoms(res_list[1]);
  BOOST_CHECK(bond_atom_list[0] == atom_list_one[0]);
  BOOST_CHECK(bond_atom_list[1] == atom_list_three[1]); 
}
BOOST_AUTO_TEST_SUITE_END();
