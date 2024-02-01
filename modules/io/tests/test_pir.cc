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
  Author: Marco Biasini
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/seq/invalid_sequence.hh>
#include <ost/io/seq/pir_io_handler.hh>
#include <ost/io/seq/load.hh>
#include <ost/io/io_exception.hh>


using namespace ost;
using namespace ost::io;


BOOST_AUTO_TEST_SUITE( io );


BOOST_AUTO_TEST_CASE(pir_filetypes) 
{
  BOOST_CHECK(PirIOHandler::ProvidesImport("","pir"));
  BOOST_CHECK(PirIOHandler::ProvidesImport("alignment.pir"));
}

BOOST_AUTO_TEST_CASE(pir_simple) 
{
  seq::SequenceHandle s=LoadSequence("testfiles/pir/simple.pir");
  BOOST_CHECK_EQUAL(s.GetName(), "name");
  BOOST_CHECK_EQUAL(s.GetString(), "ABCDEFGHIJKLMNOP");
}

BOOST_AUTO_TEST_CASE(pir_seq_types) 
{
  seq::SequenceList s=LoadSequenceList("testfiles/pir/seq-types.pir");
}

BOOST_AUTO_TEST_CASE(pir_errors)
{
  BOOST_CHECK_THROW(SequenceFromString(">", "pir"), IOException);
  BOOST_CHECK_THROW(SequenceFromString(">P1", "pir"), IOException);
  BOOST_CHECK_THROW(SequenceFromString(">P1;", "pir"), IOException);
  BOOST_CHECK_THROW(SequenceFromString(">P1;name\n", "pir"), IOException);
  BOOST_CHECK_THROW(SequenceFromString(">P1;name\ndescription\n", "pir"), 
                    IOException);  
}

BOOST_AUTO_TEST_CASE(pir_no_star) 
{
  BOOST_CHECK_THROW(LoadSequence("testfiles/pir/no-star.pir"),  
                    seq::InvalidSequence);
}

BOOST_AUTO_TEST_SUITE_END();
