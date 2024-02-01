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
 *  Authors: Marco Biasini
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <ost/seq/alg/sequence_identity.hh>



using namespace ost;
using namespace ost::seq;

BOOST_AUTO_TEST_SUITE(ost_seq_alg);

BOOST_AUTO_TEST_CASE(seqid_one) 
{
  SequenceHandle s1=CreateSequence("S1", "abcdefgh");
  SequenceHandle s2=CreateSequence("S2", "abcd----");
  
  AlignmentHandle aln=CreateAlignment();
  aln.AddSequence(s1);
  aln.AddSequence(s2);
  
  Real seqid=alg::SequenceIdentity(aln, alg::RefMode::LONGER_SEQUENCE, 0, 1);
  BOOST_CHECK_CLOSE(seqid, Real(50.0), Real(1e-6));

  // check the default sequence indices
  seqid=alg::SequenceIdentity(aln, alg::RefMode::LONGER_SEQUENCE);
  BOOST_CHECK_CLOSE(seqid, Real(50.0), Real(1e-6));

  // check the default value for reference mode
  seqid=alg::SequenceIdentity(aln);
  BOOST_CHECK_CLOSE(seqid, Real(50.0), Real(1e-6));  
  
  seqid=alg::SequenceIdentity(aln, alg::RefMode::ALIGNMENT);
  BOOST_CHECK_CLOSE(seqid, Real(100.0), Real(1e-6));  
}

BOOST_AUTO_TEST_CASE(seqid_two) 
{
  SequenceHandle s1=CreateSequence("S1", "abcd");
  SequenceHandle s2=CreateSequence("S2", "xbcx");
  
  AlignmentHandle aln=CreateAlignment();
  aln.AddSequence(s1);
  aln.AddSequence(s2);
  
  Real seqid=alg::SequenceIdentity(aln, alg::RefMode::LONGER_SEQUENCE, 0, 1);
  BOOST_CHECK_CLOSE(seqid, Real(50.0), Real(1e-6));

  // check the default sequence indices
  seqid=alg::SequenceIdentity(aln, alg::RefMode::LONGER_SEQUENCE);
  BOOST_CHECK_CLOSE(seqid, Real(50.0), Real(1e-6));

  // check the default value for reference mode
  seqid=alg::SequenceIdentity(aln);
  BOOST_CHECK_CLOSE(seqid, Real(50.0), Real(1e-6));  
  
  seqid=alg::SequenceIdentity(aln, alg::RefMode::ALIGNMENT);
  BOOST_CHECK_CLOSE(seqid, Real(50.0), Real(1e-6));  
}

BOOST_AUTO_TEST_CASE(seqid_three) 
{
  SequenceHandle s1=CreateSequence("S1", "a--d");
  SequenceHandle s2=CreateSequence("S2", "x--x");
  
  AlignmentHandle aln=CreateAlignment();
  aln.AddSequence(s1);
  aln.AddSequence(s2);
  
  Real seqid=alg::SequenceIdentity(aln, alg::RefMode::LONGER_SEQUENCE, 0, 1);
  BOOST_CHECK_CLOSE(seqid, Real(0.0), Real(1e-6));

  // check the default sequence indices
  seqid=alg::SequenceIdentity(aln, alg::RefMode::LONGER_SEQUENCE);
  BOOST_CHECK_CLOSE(seqid, Real(0.0), Real(1e-6));

  // check the default value for reference mode
  seqid=alg::SequenceIdentity(aln);
  BOOST_CHECK_CLOSE(seqid, Real(0.0), Real(1e-6));  
  
  seqid=alg::SequenceIdentity(aln, alg::RefMode::ALIGNMENT);
  BOOST_CHECK_CLOSE(seqid, Real(0.0), Real(1e-6));  
}

BOOST_AUTO_TEST_CASE(seqid_empty) 
{
  SequenceHandle s1=CreateSequence("S1", "---");
  SequenceHandle s2=CreateSequence("S2", "---");
  
  AlignmentHandle aln=CreateAlignment();
  aln.AddSequence(s1);
  aln.AddSequence(s2);
  
  Real seqid=alg::SequenceIdentity(aln, alg::RefMode::LONGER_SEQUENCE, 0, 1);
  BOOST_CHECK_CLOSE(seqid, Real(0.0), Real(1e-6));

  // check the default sequence indices
  seqid=alg::SequenceIdentity(aln, alg::RefMode::LONGER_SEQUENCE);
  BOOST_CHECK_CLOSE(seqid, Real(0.0), Real(1e-6));

  // check the default value for reference mode
  seqid=alg::SequenceIdentity(aln);
  BOOST_CHECK_CLOSE(seqid, Real(0.0), Real(1e-6));  
  
  seqid=alg::SequenceIdentity(aln, alg::RefMode::ALIGNMENT);
  BOOST_CHECK_CLOSE(seqid, Real(0.0), Real(1e-6));  
}


BOOST_AUTO_TEST_SUITE_END();
