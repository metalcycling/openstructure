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

#include <ost/geom/geom.hh>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "helper.hh"
using namespace geom;

BOOST_AUTO_TEST_SUITE( geom );


BOOST_AUTO_TEST_CASE(init_mat2)
{
  // default
  Mat2 m1;
  BOOST_CHECK(match(m1,1.0,0.0, 0.0,1.0));

  //with values
  Mat2 m2(1.1,2.2,3.3,4.4);
  BOOST_CHECK(match(m2,1.1,2.2,3.3,4.4));

  // copy ctor
  Mat2 m3(m2);
  BOOST_CHECK(match(m3,1.1,2.2,3.3,4.4));

  // array
  Real arr[]={9.9,8.8,7.7,6.6};
  Mat2 m4(arr);
  BOOST_CHECK(match(m4,arr[0],arr[1],arr[2],arr[3]));

  // assignement op
  m1=m4;
  BOOST_CHECK(match(m1,m4(0,0),m4(0,1),m4(1,0),m4(1,1)));
}

BOOST_AUTO_TEST_SUITE_END();
