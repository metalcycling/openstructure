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
#include <ost/img/alg/discrete_shrink.hh>

namespace {

using namespace ost::img;
using namespace ost::img::alg;
void test() 
{
  alg::Randomize rnd;
  DiscreteShrink discrete_shrink(Size(2,2));

  ImageHandle h=CreateImage(Extent(Point(0,0),Point(5,5)));
  h.ApplyIP(rnd);

  ImageHandle rh=h.Apply(discrete_shrink);

  BOOST_REQUIRE(rh.GetExtent()==Extent(Point(0,0),Point(2,2)));

  for(ExtentIterator it(rh.GetExtent());!it.AtEnd();++it) {
    Point p(it);
    Point p2(p[0]*2,p[1]*2);
    Real sm=0.25*(h.GetReal(p2+Point(0,0))+
      h.GetReal(p2+Point(1,0))+
      h.GetReal(p2+Point(0,1))+
      h.GetReal(p2+Point(1,1)));
    BOOST_CHECK_SMALL(Real(std::fabs(sm-rh.GetReal(it))),Real(1e-7));
  }
  
}

} // ns

test_suite* CreateDiscreteShrinkTest()
{
  test_suite* ts=BOOST_TEST_SUITE("DiscreteShrink Alg Test");

  ts->add(BOOST_TEST_CASE(&test));

  return ts;
}
