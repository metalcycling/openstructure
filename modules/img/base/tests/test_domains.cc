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

#include <sstream>

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
using boost::unit_test_framework::test_suite;

#include <ost/img/image_state/image_state_spatial_domain.hh>
#include <ost/img/image_state/image_state_frequency_domain.hh>
#include <ost/img/image_state/image_state_half_frequency_domain.hh>

namespace {

using namespace ost::img;
using namespace ost::img::image_state;

void test_spatial()
{
  std::ostringstream mesg;

  SpatialDomain sd(Extent(Size(5,6)));
  Extent le(Point(0,0),Point(4,5));
  Extent pe(Point(0,0),Point(4,5));

  BOOST_CHECK(sd.GetDomain()==SPATIAL);
  mesg.str("");
  mesg << " expected " << le << " but got " << sd.GetExtent();
  BOOST_CHECK_MESSAGE(sd.GetExtent()==le,mesg.str());
  mesg.str("");
  mesg << " expected " << le << " but got " << sd.GetLogicalExtent();
  BOOST_CHECK_MESSAGE(sd.GetLogicalExtent()==le,mesg.str());
  mesg.str("");
  mesg << " expected " << pe << " but got " << sd.GetPhysicalExtent();
  BOOST_CHECK_MESSAGE(sd.GetPhysicalExtent()==pe,mesg.str());

  sd=SpatialDomain(Extent(Size(3,4,5)));
  le=Extent(Point(0,0,0),Point(2,3,4));
  pe=Extent(Point(0,0,0),Point(2,3,4));
  mesg.str("");
  mesg << " expected " << le << " but got " << sd.GetExtent();
  BOOST_CHECK_MESSAGE(sd.GetExtent()==le,mesg.str());
  mesg.str("");
  mesg << " expected " << le << " but got " << sd.GetLogicalExtent();
  BOOST_CHECK_MESSAGE(sd.GetLogicalExtent()==le,mesg.str());
  mesg.str("");
  mesg << " expected " << pe << " but got " << sd.GetPhysicalExtent();
  BOOST_CHECK_MESSAGE(sd.GetPhysicalExtent()==pe,mesg.str());
}

void test_half_freq()
{
  std::ostringstream mesg;

  HalfFrequencyDomain hd(Extent(Size(5,6)));
  Extent le(Point(-2,-2),Point(2,3));
  Extent pe(Point(-2,0),Point(2,3));

  BOOST_CHECK(hd.GetDomain()==HALF_FREQUENCY);
  mesg.str("");
  mesg << " expected " << pe << " but got " << hd.GetExtent();
  BOOST_CHECK_MESSAGE(hd.GetExtent()==pe,mesg.str());
  mesg.str("");
  mesg << " expected " << le << " but got " << hd.GetLogicalExtent();
  BOOST_CHECK_MESSAGE(hd.GetLogicalExtent()==le,mesg.str());
  mesg.str("");
  mesg << " expected " << pe << " but got " << hd.GetPhysicalExtent();
  BOOST_CHECK_MESSAGE(hd.GetPhysicalExtent()==pe,mesg.str());

  // 3d
  hd=HalfFrequencyDomain(Extent(Size(3,4,5)));
  le=Extent(Point(-1,-1,-2),Point(1,2,2));
  pe=Extent(Point(-1,-1,0),Point(1,2,2));

  BOOST_CHECK(hd.GetDomain()==HALF_FREQUENCY);
  mesg.str("");
  mesg << " expected " << pe << " but got " << hd.GetExtent();
  BOOST_CHECK_MESSAGE(hd.GetExtent()==pe,mesg.str());
  mesg.str("");
  mesg << " expected " << le << " but got " << hd.GetLogicalExtent();
  BOOST_CHECK_MESSAGE(hd.GetLogicalExtent()==le,mesg.str());
  mesg.str("");
  mesg << " expected " << pe << " but got " << hd.GetPhysicalExtent();
  BOOST_CHECK_MESSAGE(hd.GetPhysicalExtent()==pe,mesg.str());

}


} // ns

test_suite* CreateDomainsTest()
{
  test_suite* ts=BOOST_TEST_SUITE("Domains Test");

  ts->add(BOOST_TEST_CASE(test_spatial));
  ts->add(BOOST_TEST_CASE(test_half_freq));

  return ts;
}
