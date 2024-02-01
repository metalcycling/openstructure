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

BOOST_AUTO_TEST_CASE(vecmat_mult2)
{
  Vec2 v1(1.0,2.0);
  Mat2 m(1.0,2.0, 
         3.0,4.0);

  BOOST_CHECK(match(v1*m,7.0,10.0));
  BOOST_CHECK(match(m*v1,5.0,11.0));

  Vec2 v2(5.0,6.0);
  BOOST_CHECK(match(CompMultiply(v1,v2),5.0,12.0));
  BOOST_CHECK(match(CompDivide(v2,v1),5.0,3.0));

  Mat2 m1=rnd_mat2();
  Mat2 m2=rnd_mat2();
  Mat2 m3=m1*m2;
  Vec2 v3=v2*m1;
  BOOST_CHECK(Equal(v3*m2,v2*m3));
}

BOOST_AUTO_TEST_SUITE_END();
