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

#include "test_index.hh"

#include <ost/img/image_state/index.hh>

using namespace ost::img;
using namespace ost::img::image_state;

namespace test_index {

// test both Index and IndexIterator 

void TestIterator()
{  
  IndexIterator it(5,6,7);
  std::ostringstream msg;

  BOOST_CHECK(it.AtStart());
  for(int w=0; w<7;++w) {
    for(int v=0; v<6;++v) {
      for(int u=0; u<5;++u) {
	msg.str("");
	msg << Index(u,v,w) << " != " << it;
	BOOST_REQUIRE_MESSAGE(it == Index(u,v,w),msg.str());
	++it;
      }
    }
  }
  BOOST_CHECK(it.AtEnd());

  --it.ToEnd();

  for(int w=6; w>=0;--w) {
    for(int v=5; v>=0;--v) {
      for(int u=4; u>=0;--u) {
	msg.str("");
	msg << Index(u,v,w) << " != " << it;
	BOOST_REQUIRE_MESSAGE(it == Index(u,v,w),msg.str());
	--it;
      }
    }
  }
  BOOST_CHECK(it.AtStart());

}

} // namespace 

test_suite* CreateIndexTest()
{
  using namespace test_index;
  test_suite* ts=BOOST_TEST_SUITE("Index Test");

  ts->add(BOOST_TEST_CASE(&TestIterator));

  return ts;
}
