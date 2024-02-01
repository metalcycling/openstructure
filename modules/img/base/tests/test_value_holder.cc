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

#include <boost/mpl/vector.hpp>

#include "test_value_holder.hh"

#include <ost/img/image_state/value_holder.hh>
#include <ost/img/image_state/index.hh>
#include <ost/img/value_util.hh>

using namespace ost::img;
using namespace ost::img::image_state;

namespace value_holder_test {


template <typename V>
void ReadWrite()
{
  int wi=4;
  int he=3;
  int de=5;
  ValueHolder<V> vh(Size(wi,he,de));

  for(int u=0;u<wi;++u) {
    for(int v=0;v<he;++v) {
      for(int w=0;w<de;++w) {
	V val = Random<V>();
	int indx = w+de*(v+he*u);
	vh.Value(Index(u,v,w))=val;
	BOOST_REQUIRE(vh.Value(Index(u,v,w)) == val);
	BOOST_REQUIRE(vh.GetData()[indx] == val);
      }
    }
  }
}

template <typename V>
void CopyCtor()
{
  int wi=4;
  int he=3;
  int de=5;

  ValueHolder<V> vh1(Size(wi,he,de));

  for(int u=0;u<wi;++u) {
    for(int v=0;v<he;++v) {
      for(int w=0;w<de;++w) {
	V val = Random<V>();
	vh1.Value(Index(u,v,w))=val;
      }
    }
  }

  ValueHolder<V> vh2(vh1);
  for(int u=0;u<wi;++u) {
    for(int v=0;v<he;++v) {
      for(int w=0;w<de;++w) {
	BOOST_REQUIRE(vh2.Value(Index(u,v,w))==vh1.Value(Index(u,v,w)));
      }
    }
  }
}

}

test_suite* CreateValueHolderTest()
{
  using namespace value_holder_test;
  
  test_suite* ts=BOOST_TEST_SUITE("ValueHolder Test");

  ts->add( BOOST_TEST_CASE(ReadWrite<Real>) );
  ts->add( BOOST_TEST_CASE(ReadWrite<Complex>) );
  ts->add( BOOST_TEST_CASE(CopyCtor<Real>) );
  ts->add( BOOST_TEST_CASE(CopyCtor<Complex>) );
  

  return ts;
}
