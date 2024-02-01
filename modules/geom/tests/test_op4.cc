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

BOOST_AUTO_TEST_CASE(vecmat_mult4)
{
  Vec4 v1(5.0,4.0,3.0,2.0);
  Mat4 m(0.5,1.0,1.5,2.0, 
	 2.5,3.0,3.5,4.0,
	 4.5,5.0,5.5,6.0,
	 6.5,7.0,7.5,8.0);

  BOOST_CHECK(match(v1*m,39.0,46.0,53.0,60.0));
  BOOST_CHECK(match(m*v1,15.0,43.0,71.0,99.0));

  Vec4 v2(10.0,12.0,9.0,4.0);
  BOOST_CHECK(match(CompMultiply(v1,v2),50.0,48.0,27.0,8.0));
  BOOST_CHECK(match(CompDivide(v2,v1),2.0,3.0,3.0,2.0));
}

BOOST_AUTO_TEST_SUITE_END();
