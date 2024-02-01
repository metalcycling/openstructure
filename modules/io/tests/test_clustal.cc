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

#include <ost/io/seq/clustal_io_handler.hh>
#include <ost/io/seq/load.hh>

using namespace ost;
using namespace ost::io;
const static String MOUSE="MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTAITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS";
const static String HUMAN="MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTAITTSQDLQWLVQPTLISSMAQSQGQPLASQPPVVDPYDMPGTSYSTPGMSGYSSGGASGS";


BOOST_AUTO_TEST_SUITE( io );


BOOST_AUTO_TEST_CASE(clustal) 
{
  BOOST_CHECK(ClustalIOHandler::ProvidesImport("","clustal"));
  BOOST_CHECK(ClustalIOHandler::ProvidesImport("alignment.aln"));
  seq::AlignmentHandle aln=LoadAlignment("testfiles/alignment.aln");
  BOOST_CHECK_EQUAL(aln.GetCount(), 2);
  BOOST_CHECK_EQUAL(aln.FindSequence("FOSB_MOUSE").GetString(), 
                    aln.GetSequence(0).GetString());
  BOOST_CHECK_EQUAL(aln.FindSequence("FOSB_HUMAN").GetString(), 
                    aln.GetSequence(1).GetString());  
  BOOST_CHECK_EQUAL(aln.GetSequence(0).GetString(), MOUSE);
  BOOST_CHECK_EQUAL(aln.GetSequence(1).GetString(), HUMAN);
}

BOOST_AUTO_TEST_SUITE_END();
