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

BOOST_AUTO_TEST_CASE(init_mat4)
{
  // default
  Mat4 m1;
  BOOST_CHECK(match(m1,
        1.0,0.0,0.0,0.0, 
        0.0,1.0,0.0,0.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1.0));

  //with values
  Mat4 m2(1.1,2.2,3.3,4.4,
    5.5,6.6,7.7,8.8,
    9.9,1.1,2.2,3.3,
    4.4,5.5,6.6,7.7);
  BOOST_CHECK(match(m2,
        1.1,2.2,3.3,4.4,
        5.5,6.6,7.7,8.8,
        9.9,1.1,2.2,3.3,
        4.4,5.5,6.6,7.7));

  // copy ctor
  Mat4 m3(m2);
  BOOST_CHECK(match(m2,
        1.1,2.2,3.3,4.4,
        5.5,6.6,7.7,8.8,
        9.9,1.1,2.2,3.3,
        4.4,5.5,6.6,7.7));

  // array
  double arr[]={9.9,8.8,7.7,6.6,5.5,4.4,3.3,2.2,1.1,9.9,8.8,7.7,6.6,5.5,4.4,3.3};
  Mat4 m4(arr);
  BOOST_CHECK(match(m4,
        arr[0],arr[1],arr[2],arr[3],
        arr[4],arr[5],arr[6],arr[7],
        arr[8],arr[9],arr[10],arr[11],
        arr[12],arr[13],arr[14],arr[15]));

  // assignement op
  m1=m4;
  BOOST_CHECK(match(m1,
        m4(0,0),m4(0,1),m4(0,2),m4(0,3),
        m4(1,0),m4(1,1),m4(1,2),m4(1,3),
        m4(2,0),m4(2,1),m4(2,2),m4(2,3),
        m4(3,0),m4(3,1),m4(3,2),m4(3,3)));

  // conversion
  Mat4 m5(Mat2(9.0,8.0,7.0,6.0));
  BOOST_CHECK(match(m5,
        9.0,8.0,0.0,0.0,
        7.0,6.0,0.0,0.0,
        0.0,0.0,1.0,0.0,
        0.0,0.0,0.0,1.0));

  Mat4 m6(Mat3(1.0,2.0,3.0, 4.0,5.0,6.0, 7.0,8.0,9.0));
  BOOST_CHECK(match(m6,
        1.0,2.0,3.0,0.0,
        4.0,5.0,6.0,0.0,
        7.0,8.0,9.0,0.0,
        0.0,0.0,0.0,1.0));
        
}

BOOST_AUTO_TEST_SUITE_END();
