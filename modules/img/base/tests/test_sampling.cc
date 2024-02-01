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

#include "test_sampling.hh"

#include <ost/img/pixel_sampling.hh>
#include <ost/img/extent.hh>
#include <ost/img/vecmat.hh>

using namespace ost::img;

namespace test_sampling {

void test()
{
  PixelSampling ps(geom::Vec3(2.0,4.0,1.0),SPATIAL,Extent(Size(4,10,20)));
  BOOST_CHECK(geom::Length(ps.GetPixelSampling()-geom::Vec3(2.0,4.0,1.0))<1e-5);
  BOOST_CHECK(geom::Length(ps.GetSpatialSampling()-geom::Vec3(2.0,4.0,1.0))<1e-5);
  ps.SetDomain(FREQUENCY);
  BOOST_CHECK(geom::Length(ps.GetPixelSampling()-geom::Vec3(0.125,0.025,0.05))<1e-5);
  BOOST_CHECK(geom::Length(ps.GetSpatialSampling()-geom::Vec3(2.0,4.0,1.0))<1e-5);

  PixelSampling ps2(geom::Vec3(4.0,2.0,8.0),SPATIAL,Extent(Size(4,10,20)));
  BOOST_CHECK(geom::Length(ps2.Point2Coord(Point(2,3,5))-geom::Vec3(8.0,6.0,40.0))<1e-5);
  ps2.SetDomain(FREQUENCY);
  BOOST_CHECK(geom::Length(ps2.Point2Coord(Point(1,4,7))-geom::Vec3(0.0625,0.2,0.04375))<1e-5);
}

} // namespace 

test_suite* CreateSamplingTest()
{
  using namespace test_sampling;
  test_suite* ts=BOOST_TEST_SUITE("Sampling Test");

  ts->add(BOOST_TEST_CASE(&test));

  return ts;
}
