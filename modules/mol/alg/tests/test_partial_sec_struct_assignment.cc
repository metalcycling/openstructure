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


#include <ost/mol/alg/sec_struct.hh>
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/mol/mol.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/mol/builder.hh>
#include <ost/message.hh>



BOOST_AUTO_TEST_SUITE( mol_alg );

BOOST_AUTO_TEST_CASE(partial_sec_struct_assignment) {

  // load a test structure
  String filepath = "testfiles/1a0s.pdb";
  ost::io::PDBReader reader(filepath, ost::io::IOProfile());
  ost::mol::EntityHandle test_ent = ost::mol::CreateEntity();
  reader.Import(test_ent);

  ost::mol::EntityView view = test_ent.CreateFullView();
  ost::mol::ResidueViewList res_list = view.GetResidueList();

  std::vector<int> res_indices;
  std::vector<geom::Vec3> ca_positions;
  std::vector<int> donor_for_one;
  std::vector<int> donor_for_two;
  std::vector<int> connected_to_next;

  ost::mol::alg::PrepareSSData(res_list, res_indices, ca_positions,
                               donor_for_one, donor_for_two,
                               connected_to_next);

  String full_pred = ost::mol::alg::RawEstimateSS(ca_positions,
                                                  0, ca_positions.size(),
                                                  donor_for_one, donor_for_two,
                                                  connected_to_next);

  // select some start indices to test
  int max_num = ca_positions.size();
  std::vector<int> test_locations;
  test_locations.push_back(max_num/10);
  test_locations.push_back((max_num/10)*2);
  test_locations.push_back((max_num/10)*3);
  test_locations.push_back((max_num/10)*4);
  test_locations.push_back((max_num/10)*5);
  test_locations.push_back((max_num/10)*6);
  test_locations.push_back((max_num/10)*7);
  test_locations.push_back((max_num/10)*8);
  test_locations.push_back((max_num/10)*9);
  int stretch_length = std::min(12, max_num - test_locations.back());

  for(uint i = 0; i < test_locations.size(); ++i) {
    String expected_ss = full_pred.substr(test_locations[i], stretch_length);
    String calculated_ss = ost::mol::alg::RawEstimateSS(ca_positions,
                                                        test_locations[i],
                                                        stretch_length,
                                                        donor_for_one,
                                                        donor_for_two,
                                                        connected_to_next);

    BOOST_CHECK(expected_ss == calculated_ss);
  }

}

BOOST_AUTO_TEST_SUITE_END();
