//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
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


#include <ost/geom/geom.hh>


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "helper.hh"
using namespace geom;

BOOST_AUTO_TEST_SUITE( geom );


BOOST_AUTO_TEST_CASE(line_init2)
{
  Line3 line(geom::Vec2(0,0), geom::Vec2(2,0));
  BOOST_CHECK_EQUAL(geom::Length(line.GetDirection()), 1.0);
}

BOOST_AUTO_TEST_CASE(rectangle2)
{
  Rectangle2 rect;
  try{
    BOOST_CHECK_THROW(rect[2]=geom::Vec2(0,0), OutOfRangeException);
  }catch(...)
  {
    BOOST_ERROR("Failed to catch OutOfRangeException");
  }
}

BOOST_AUTO_TEST_CASE(is_on_line2)
{
  Line2 line(geom::Vec2(0,0), geom::Vec2(1,0));
  BOOST_CHECK(IsOnLine(line, geom::Vec2(0.5,0.0)));
  BOOST_CHECK(IsOnLine(line, geom::Vec2(1.0,0.0)));
  BOOST_CHECK(IsOnLine(line, geom::Vec2(0.0,0.0)));
  BOOST_CHECK(IsOnLine(line, geom::Vec2(-5,0.0)));
  BOOST_CHECK(IsOnLine(line, geom::Vec2(10.0,0.0)));
  BOOST_CHECK(!IsOnLine(line, geom::Vec2(0.5,0.1)));
  BOOST_CHECK(!IsOnLine(line, geom::Vec2(1.0,0.1)));
  
  line=Line2(geom::Vec2(1,0), geom::Vec2(3,2));
  for (int i=-10; i<10; ++i) {
    BOOST_CHECK(IsOnLine(line, line.At(i)));
  }
  BOOST_CHECK(!IsOnLine(line, geom::Vec2(3,2.1)));
}

BOOST_AUTO_TEST_SUITE_END();
