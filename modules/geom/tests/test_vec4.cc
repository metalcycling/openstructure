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

BOOST_AUTO_TEST_CASE(init_vec4)
{
  // default
  Vec4 v1;
  BOOST_CHECK(match(v1,0.0,0.0,0.0,1.0));

  //with values
  Vec4 v2(5.9,1.5,2.3,4.5);
  BOOST_CHECK(match(v2,5.9,1.5,2.3,4.5));

  // copy ctor
  Vec4 v3(v2);
  BOOST_CHECK(match(v3,5.9,1.5,2.3,4.5));

  // array
  double arr[]={3.3,6.5,4.3,2.1};
  Vec4 v4(arr);
  BOOST_CHECK(match(v4,arr[0],arr[1],arr[2],arr[3]));

  // assignement op
  v1=v4;
  BOOST_CHECK(match(v1,v4[0],v4[1],v4[2],v4[3]));

  // conversion from vec2
  Vec4 v5(Vec2(2.8,1.5));
  BOOST_CHECK(match(v5,2.8,1.5,0.0,1.0));

  // conversion from vec4
  Vec4 v6(Vec3(0.4,1.2,4.0));
  BOOST_CHECK(match(v6,0.4,1.2,4.0,1.0));
}


BOOST_AUTO_TEST_CASE(operators_vec4)
{
  Vec4 v1(0.8,1.2,2.3,3.4);
  Vec4 v2(0.8,1.2,2.3,3.4);
  Vec4 v3(0.2,0.1,0.2,0.3);

  // equality
  BOOST_CHECK(v1==v2);
  BOOST_CHECK(v1!=v3);

  // additive
  v2=v1+v3;
  BOOST_CHECK(match(v2,1.0,1.3,2.5,3.7));
  v2=v1-v3;
  BOOST_CHECK(match(v2,0.6,1.1,2.1,3.1));

  // multipliable with scalar
  v2=v1*2.0;
  BOOST_CHECK(match(v2,1.6,2.4,4.6,6.8));
  
  // dividable with scalar
  v2=v1/2.0;
  BOOST_CHECK(match(v2,0.4,0.6,1.15,1.7));

  // negative
  v2=-v1;
  BOOST_CHECK(match(v2,-0.8,-1.2,-2.3,-3.4));

  // reciprocal
  v2=1.0/v1;
  BOOST_CHECK(match(v2,1.0/0.8,1.0/1.2,1.0/2.3,1.0/3.4));

  // chaining
  v2=-(3.0*v3-v1/2.0);
  BOOST_CHECK(match(v2,-0.2,0.3,0.55,0.8));
}

BOOST_AUTO_TEST_SUITE_END();
