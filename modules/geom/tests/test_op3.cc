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

BOOST_AUTO_TEST_CASE(vecmat_mult3)
{
  Vec3 v1(1.0,2.0,3.0);
  Mat3 m(1.0,2.0,3.0, 
         4.0,5.0,6.0, 
         7.0,8.0,9.0);

  BOOST_CHECK(match(v1*m,30.0,36.0,42.0));
  BOOST_CHECK(match(m*v1,14.0,32.0,50.0));

  Vec3 v2(7.0,8.0,9.0);
  BOOST_CHECK(match(CompMultiply(v1,v2),7.0,16.0,27.0));
  BOOST_CHECK(match(CompDivide(v2,v1),7.0,4.0,3.0));

  Mat3 m1=m;
  Mat3 m2(-4,0,1,0,8,-10,8,3,0);
  Mat3 m3=m1*m2;
  Vec3 v3=v2*m1;
  BOOST_CHECK(Equal(v3*m2,v2*m3));
}

BOOST_AUTO_TEST_CASE(vec3list_op)
{
  Vec3List l1,l2;
  BOOST_CHECK_THROW(MinDistance(l1,l2),GeomException);
  l1.push_back(Vec3(0.0,0.0,0.0));
  l1.push_back(Vec3(1.0,0.0,0.0));
  BOOST_CHECK_THROW(MinDistance(l1,l2),GeomException);
  BOOST_CHECK_THROW(MinDistance(l2,l1),GeomException);
  l2.push_back(Vec3(4.0,4.0,0.0));
  l2.push_back(Vec3(2.0,4.0,5.0));
  l2.push_back(Vec3(0.0,3.0,6.0));
  BOOST_CHECK(MinDistance(l1,l2)==5.0);
}

BOOST_AUTO_TEST_SUITE_END();
