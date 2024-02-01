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

#include "tests.hh"

#include <ost/img/image.hh>
#include <ost/img/alg/randomize.hh>
#include <ost/img/alg/transform.hh>

namespace {

using namespace ost::img;

void test() 
{
  alg::Randomize rnd;
  alg::Transform transform;
}

void test_point()
{
  Point p0(3,-2,1);
  Point p1a(9,-6,3);
  Point p1b(alg::Scale(3.0).Apply(p0));
  Point p2a(Point(2,3,1));
  Point p2b=alg::Rotate2D(M_PI_2).Apply(p0);

  std::ostringstream msg;
  msg << p1a << " != " << p1b;
  BOOST_CHECK_MESSAGE(p1a==p1b,msg.str());
  msg.str("");
  msg << p2a << " != " << p2b;
  BOOST_CHECK_MESSAGE(p2a==p2b,msg.str());
}

void test_extent()
{
  Extent e0(Point(0,0),Point(2,3));
  Extent e1a(Point(-1,-1),Point(7,10));
  Extent e1b=alg::Scale(3.0).Apply(e0);
  Extent e2a(Point(-3,0),Point(0,2));
  Extent e2b=alg::Rotate2D(M_PI_2).Apply(e0);

  std::ostringstream msg;
  msg << e1a << " != " << e1b;
  BOOST_CHECK_MESSAGE(e1a==e1b,msg.str());
  msg.str("");
  msg << e2a << " != " << e2b;
  BOOST_CHECK_MESSAGE(e2a==e2b,msg.str());

  e0=Extent(Point(-3,-2),Size(8,10));
  e1a=Extent(Point(-10,-7),Size(24,30));
  e1b=alg::Scale(3.0).Apply(e0);

  e2a=Extent(Point(-7,-3),Size(10,8));
  e2b=alg::Rotate2D(M_PI_2).Apply(e0);

  msg.str("");
  msg << e1a << " != " << e1b;
  BOOST_CHECK_MESSAGE(e1a==e1b,msg.str());
  msg.str("");
  msg << e2a << " != " << e2b;
  BOOST_CHECK_MESSAGE(e2a==e2b,msg.str());
}

} // ns

test_suite* CreateTransformTest()
{
  test_suite* ts=BOOST_TEST_SUITE("tf transform Test");

  ts->add(BOOST_TEST_CASE(&test));
  ts->add(BOOST_TEST_CASE(&test_point));
  ts->add(BOOST_TEST_CASE(&test_extent));

  return ts;
}
