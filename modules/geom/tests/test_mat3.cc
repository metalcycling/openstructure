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

BOOST_AUTO_TEST_CASE(init_mat3)
{
  // default
  Mat3 m1;
  BOOST_CHECK(match(m1,1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0));

  //with values
  Mat3 m2(1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9);
  BOOST_CHECK(match(m2,1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9));

  // copy ctor
  Mat3 m3(m2);
  BOOST_CHECK(match(m3,1.1,2.2,3.3,4.4,5.5,6.6,7.7,8.8,9.9));

  // array
  Real arr[]={9.9,8.8,7.7,6.6,5.5,4.4,3.3,2.2,1.1};
  Mat3 m4(arr);
  BOOST_CHECK(match(m4,arr[0],arr[1],arr[2],arr[3],arr[4],arr[5],arr[6],arr[7],arr[8]));

  // assignement op
  m1=m4;
  BOOST_CHECK(match(m1,m4(0,0),m4(0,1),m4(0,2),m4(1,0),m4(1,1),m4(1,2),m4(2,0),m4(2,1),m4(2,2)));

  Mat3 m5(Mat2(2.0,3.0,4.0,5.0));
  BOOST_CHECK(match(m5,2.0,3.0,0.0,4.0,5.0,0.0,0.0,0.0,1.0));
}

BOOST_AUTO_TEST_CASE(mult_mat3)
{
  Mat3 m1=rnd_mat3();
  Mat3 m2=rnd_mat3();
  Mat3 m3=m1*m2;
  Mat3 m4;
  for(int ir=0;ir<3;++ir) {
    for(int ic=0;ic<3;++ic) {
      double tmp=0.0;
      for(int tt=0;tt<3;++tt) {
        tmp+=m1(ir,tt)*m2(tt,ic);
      }
      m4(ir,ic)=tmp;
    }
  }
  BOOST_CHECK(Equal(m3,m4));

  m4=rnd_mat3();
  Mat3 m5=m2*m4;

  BOOST_CHECK(Equal(m1*m5,m3*m4));  
}

BOOST_AUTO_TEST_SUITE_END();
