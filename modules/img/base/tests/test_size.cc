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

#include "test_size.hh"

#include <ost/img/size.hh>
using namespace ost::img;

namespace test_size {

void Construction()
{
  // default
  Size s0;
  BOOST_CHECK(s0[0]==1);
  BOOST_CHECK(s0[1]==1);
  BOOST_CHECK(s0[2]==1);
  BOOST_CHECK(s0.GetDim()==1);

  // 1D
  Size s1(4);
  BOOST_CHECK(s1[0]==4);
  BOOST_CHECK(s1[1]==1);
  BOOST_CHECK(s1[2]==1);
  BOOST_CHECK(s1.GetDim()==1);
  
  // 2D
  Size s2(7,2);
  BOOST_CHECK(s2[0]==7);
  BOOST_CHECK(s2[1]==2);
  BOOST_CHECK(s2[2]==1);
  BOOST_CHECK(s2.GetDim()==2);
  
  // 3D
  Size s3(3,4,2);
  BOOST_CHECK(s3[0]==3);
  BOOST_CHECK(s3[1]==4);
  BOOST_CHECK(s3[2]==2);
  BOOST_CHECK(s3.GetDim()==3);
}

} // namespace 

test_suite* CreateSizeTest()
{
  using namespace test_size;
  test_suite* ts=BOOST_TEST_SUITE("Size Test");

  ts->add(BOOST_TEST_CASE(&Construction));

  return ts;
}
