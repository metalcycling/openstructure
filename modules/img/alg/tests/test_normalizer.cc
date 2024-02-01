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
  Author: Andreas Schenk
*/

#include <iostream>

#include "tests.hh"
#include <boost/version.hpp>
#if BOOST_VERSION<105900
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#include <ost/img/image.hh>
#include  <ost/img/alg/normalizer_factory.hh>
#include  <ost/img/alg/stat.hh>



namespace {

using namespace ost::img;
using namespace ost::img::alg;

void test() 
{
#if BOOST_VERSION<105900
  boost::test_tools::close_at_tolerance<Real> close_test(boost::test_tools::percent_tolerance(0.001));
#else
  boost::math::fpc::close_at_tolerance<Real> close_test(boost::math::fpc::percent_tolerance(0.001));
#endif
  ost::img::ImageHandle testimage=ost::img::CreateImage(ost::img::Extent(ost::img::Point(0,0),ost::img::Point(3,3)));
  int counter=0;
  for (ost::img::ExtentIterator i(testimage.GetExtent()); !i.AtEnd(); ++i, ++counter) {
   testimage.SetReal(i, counter);
  }
  ost::img::alg::Normalizer norm=ost::img::alg::CreateLinearRangeNormalizer(testimage,0.0,65535.0);
  ost::img::ImageHandle scaled_image=testimage.Apply(norm);
  scaled_image+=0.01; //if all values are > 0.0 we can use close_at_tolerance
  bool failed=false;
  ost::img::ExtentIterator eit(testimage.GetExtent());
  for(;!eit.AtEnd();++eit) {
    if( ! close_test(scaled_image.GetReal(eit),testimage.GetReal(eit)/15.0*65535.0+0.01)){
      failed=true;
      break;
    }
  }
   if(failed){
    BOOST_ERROR("Normalizer failed at point " 
                << ost::img::Point(eit)<< ". Should be " 
                << testimage.GetReal(eit)/15.0*65535.0+0.01 << ", but "
                << scaled_image.GetReal(eit) << " found.");
    }
}

} // ns

test_suite* CreateNormalizerTest()
{
  test_suite* ts=BOOST_TEST_SUITE("img alg Normalizer Test");

  ts->add(BOOST_TEST_CASE(&test));

  return ts;
}
