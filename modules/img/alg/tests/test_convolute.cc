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
#include <ost/img/alg/convolute.hh>
#include <ost/img/alg/fft.hh>
#include <ost/img/alg/dft.hh>

namespace {

using namespace ost::img;
using namespace ost::img::alg;

void test() 
{
  ImageHandle im1=CreateImage(Extent(Point(0),Size(8)),REAL);
  ImageHandle im2=CreateImage(Extent(Point(0),Size(8)),REAL);
  im1.ApplyIP(alg::Randomize());
  im2.ApplyIP(alg::Randomize());
#if 0
  ImageHandle im3=im1.Apply(alg::ExplicitConvolute(im2));

  BOOST_REQUIRE(im3.GetExtent()==im1.GetExtent());

  im1.ApplyIP(alg::DFT());
  im2.ApplyIP(alg::DFT());

  im1*=im2;
  im1.ApplyIP(alg::DFT());

  for(ExtentIterator it(im3.GetExtent()); !it.AtEnd(); ++it) {
    //std::cerr << it << " " << im1.GetReal(it) << " " << im3.GetReal(it) << std::endl;
    BOOST_REQUIRE(std::fabs(im3.GetReal(it)-im1.GetReal(it))<1e-10);
  }
#endif
}

} // ns

test_suite* CreateConvoluteTest()
{
  test_suite* ts=BOOST_TEST_SUITE("Convolute Alg Test");

  ts->add(BOOST_TEST_CASE(&test));

  return ts;
}
