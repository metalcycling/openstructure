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



#include <ost/img/alg/histogram.hh>

#include "tests.hh"
#include <ost/img/image.hh>

using namespace ost::img;
using namespace ost::img::alg;

namespace test_histogram {

void test()
{
  ImageHandle i = CreateImage(Size(10));
  for (int k=0;k<10;++k) {
    i.SetReal(Point(k),Real(k));
  }

  Histogram h(10,0.0,9.0);

  Histogram::Bins bins = h.GetBins();
  for (int k=0;k<10;++k) {
    BOOST_CHECK(bins[k]==0);
  }

  i.Apply(h);

  bins = h.GetBins();

  std::ostringstream msg;

  for (int k=0;k<10;++k) {
    msg.str("");
    msg << " expected bins[" << k << "] == 1 but got " << bins[k];
    BOOST_CHECK_MESSAGE(bins[k]==1,msg.str());
  }
  try {
      BOOST_CHECK_THROW(Histogram(1,0.0,1.0),HistogramError);
  }catch(...){
    BOOST_ERROR( "Failed to catch HisotgramException." );
  }


}

} // namespace 

test_suite* CreateHistogramTest()
{
  using namespace test_histogram;
  test_suite* ts=BOOST_TEST_SUITE("Histogram Test");

  ts->add(BOOST_TEST_CASE(&test));

  return ts;
}
