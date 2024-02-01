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
 
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <ost/mol/surface_handle.hh>
#include <cmath>

using namespace ost;
using namespace ost::mol;

BOOST_AUTO_TEST_SUITE( mol_base );


BOOST_AUTO_TEST_CASE(test_throw_invalid_surface)
{
  SurfaceHandle surf;
  BOOST_CHECK_THROW(CheckHandleValidity(surf), InvalidHandle);
  
  surf=CreateSurface();
  BOOST_CHECK_NO_THROW(CheckHandleValidity(surf));
}

BOOST_AUTO_TEST_SUITE_END();
