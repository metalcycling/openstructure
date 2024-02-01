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

/*
  Author: Andreas Schenk
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <ost/io/img/load_map.hh>
#include <ost/img/image.hh>
#include <ost/io/io_exception.hh>

using namespace ost;
using namespace ost::io;


BOOST_AUTO_TEST_SUITE( io );


BOOST_AUTO_TEST_CASE(exception) 
{
  try {
      BOOST_CHECK_THROW(ost::img::ImageHandle loadedimage=LoadImage("nonexistent.mrc"),IOException);
  }catch(...){
    BOOST_ERROR( "Failed to catch IOException." );
  }
  ost::img::ImageHandle ih=ost::img::CreateImage();
  try {
      BOOST_CHECK_THROW(SaveImage(ih,"nonexistent.ABC"),IOUnknownFormatException);
  }catch(...){
    BOOST_ERROR( "Failed to catch IOUnknownFormatException." );
  }
}

BOOST_AUTO_TEST_SUITE_END();
