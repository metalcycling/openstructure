//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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

/*
  Author: Ansgar Philippsen
*/

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/framework.hpp>
#include <ost/base.hh>

#include "test_data.hh"
#include "test_extent.hh"
#include "test_image.hh"
#include "test_image_state.hh"
#include "test_image_state_visitor.hh"
#include "test_index.hh"
#include "test_point.hh"
#include "test_sampling.hh"
#include "test_size.hh"

#include "test_value_holder.hh"

using namespace boost::unit_test;

test_suite* CreateEntityTest();
test_suite* CreateQueryTest();
test_suite* CreateConnTest();
test_suite* CreateICSTest();
test_suite* InsertChainTest();
test_suite* CreateDeleteTest();
test_suite* CreateViewTest();
test_suite* CreateAtomGroupsTest();
test_suite* CreateIteratorsTest();
test_suite* CreateResidueTest();
test_suite* CreateCoordGroupTest();
test_suite* CreateDomainsTest();

bool init_ost_img_base_unit_tests() {
  try {
    framework::master_test_suite().add(CreateDataTest());
    framework::master_test_suite().add(CreateExtentTest());
    framework::master_test_suite().add(CreateImageTest());
    framework::master_test_suite().add(CreateImageStateTest());
    framework::master_test_suite().add(CreateImageStateVisitorTest());
    framework::master_test_suite().add(CreateIndexTest());
    framework::master_test_suite().add(CreatePointTest());
    framework::master_test_suite().add(CreateSamplingTest());
    framework::master_test_suite().add(CreateSizeTest());
    framework::master_test_suite().add(CreateValueHolderTest());
    framework::master_test_suite().add(CreateDomainsTest());
  } catch(std::exception& e) {
    return false;
  }
  return true;
}


int main(int argc, char** argv)
{
  return boost::unit_test::unit_test_main(init_ost_img_base_unit_tests,
                                          argc, argv);
}
