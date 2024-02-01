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

#include "test_point.hh"

#include <ost/img/point.hh>
using namespace ost::img;

namespace test_point {

void Construction()
{
  // default
  Point p0;
  BOOST_CHECK(p0[0]==0);
  BOOST_CHECK(p0[1]==0);
  BOOST_CHECK(p0[2]==0);

  // 1D
  Point p1(3);
  BOOST_CHECK(p1[0]==3);
  BOOST_CHECK(p1[1]==0);
  BOOST_CHECK(p1[2]==0);

  // 2D
  Point p2(-2,8);
  BOOST_CHECK(p2[0]==-2);
  BOOST_CHECK(p2[1]==8);
  BOOST_CHECK(p2[2]==0);

  // 3D
  Point p3(4,-1,5);
  BOOST_CHECK(p3[0]==4);
  BOOST_CHECK(p3[1]==-1);
  BOOST_CHECK(p3[2]==5);
}

} // namespace 

test_suite* CreatePointTest()
{
  using namespace test_point;
  test_suite* ts=BOOST_TEST_SUITE("Point Test");

  ts->add(BOOST_TEST_CASE(&Construction));

  return ts;
}
