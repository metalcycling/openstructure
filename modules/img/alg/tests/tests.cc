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

#include <iostream>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/test/framework.hpp>

using namespace boost::unit_test;


extern test_suite* CreateShiftTest();
extern test_suite* CreateMirrorTest();
extern test_suite* CreateTransformTest();
extern test_suite* CreateDiscreteShrinkTest();
extern test_suite* CreateFillTest();
extern test_suite* CreateStatTest();
extern test_suite* CreateClearTest();
extern test_suite* CreateFFTTest();
extern test_suite* CreateNegateTest();
extern test_suite* CreateConjugateTest();
extern test_suite* CreateNormalizerTest();

bool init_ost_img_alg_unit_tests() {
  try {
    framework::master_test_suite().add(CreateShiftTest()); 
    framework::master_test_suite().add(CreateMirrorTest()); 
    framework::master_test_suite().add(CreateTransformTest()); 
    framework::master_test_suite().add(CreateDiscreteShrinkTest()); 
    framework::master_test_suite().add(CreateFillTest());
    framework::master_test_suite().add(CreateStatTest()); 
    framework::master_test_suite().add(CreateConjugateTest()); 
    framework::master_test_suite().add(CreateClearTest());      
    framework::master_test_suite().add(CreateNegateTest());
    framework::master_test_suite().add(CreateFFTTest());          
    framework::master_test_suite().add(CreateNormalizerTest());          
  } catch(std::exception& e) {
    return false;
  }
  return true;
}


int main(int argc, char** argv)
{
  return boost::unit_test::unit_test_main(init_ost_img_alg_unit_tests, 
                                          argc, argv);
}
