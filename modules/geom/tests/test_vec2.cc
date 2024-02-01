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
  
BOOST_AUTO_TEST_CASE(init_vec2)
{
  // default
  Vec2 v1;
  BOOST_CHECK(match(v1,0.0,0.0));

  //with values
  Vec2 v2(1.5,2.3);
  BOOST_CHECK(match(v2,1.5,2.3));

  // copy ctor
  Vec2 v3(v2);
  BOOST_CHECK(match(v3,1.5,2.3));

  // array
  double arr[]={6.5,4.3};
  Vec2 v4(arr);
  BOOST_CHECK(match(v4,arr[0],arr[1]));

  // assignement op
  v1=v4;
  BOOST_CHECK(match(v1,v4[0],v4[1]));

  // conversion from vec3
  Vec2 v5(Vec3(1.2,2.8,1.5));
  BOOST_CHECK(match(v5,1.2,2.8));

  // conversion from vec4
  Vec2 v6(Vec4(0.5,2.0,4.0,2.0));
  BOOST_CHECK(match(v6, 0.25, 1.0));
  BOOST_CHECK_THROW( Vec2(Vec4(1.0,1.0,1.0,0.0)), DivideByZeroException);
}

BOOST_AUTO_TEST_CASE(operators_vec2)
{
  Vec2 v1(1.2,2.3);
  Vec2 v2(1.2,2.3);
  Vec2 v3(0.1,0.2);

  // equality
  BOOST_CHECK(v1==v2);
  BOOST_CHECK(v1!=v3);

  // additive
  v2=v1+v3;
  BOOST_CHECK(match(v2,1.3,2.5));
  v2=v1-v3;
  BOOST_CHECK(match(v2,1.1,2.1));

  // multipliable with scalar
  v2=v1*2.0;
  BOOST_CHECK(match(v2,2.4,4.6));
  
  // dividable with scalar
  v2=v1/2.0;
  BOOST_CHECK(match(v2,0.6,1.15));

  // negative
  v2=-v1;
  BOOST_CHECK(match(v2,-1.2,-2.3));

  // reciprocal
  v2=1.0/v1;
  BOOST_CHECK(match(v2,1.0/1.2,1.0/2.3));

  // chaining
  v2=-(3.0*v3-v1/2.0);
  BOOST_CHECK(match(v2,0.3,0.55));
}

BOOST_AUTO_TEST_SUITE_END();
