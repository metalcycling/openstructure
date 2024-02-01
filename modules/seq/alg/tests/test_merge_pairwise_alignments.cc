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
 *  Authors: Marco Biasini, Juergen Haas
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <ost/seq/alg/merge_pairwise_alignments.hh>
#include <ost/integrity_error.hh>


using namespace ost;
using namespace ost::seq;

BOOST_AUTO_TEST_SUITE(ost_seq_alg);

BOOST_AUTO_TEST_CASE(merge_pairwise_alignments_one) 
{
  SequenceHandle ref=CreateSequence("REF", "abcdefghijklmn", "R");
  SequenceHandle s1=CreateSequence("S1", "abcd---efghijklmn");
  SequenceHandle s2=CreateSequence("S2", "abcdxyzefghijklmn", "X");
  AlignmentHandle aln1=CreateAlignment();
  aln1.AddSequence(s1);
  aln1.AddSequence(s2);

  SequenceHandle s3=CreateSequence("S1", "abcdefghij---klmn");
  SequenceHandle s4=CreateSequence("S2", "abcdefghijxyzklmn", "Y");

  AlignmentHandle aln2=CreateAlignment();
  aln2.AddSequence(s3);
  aln2.AddSequence(s4);
  AlignmentList l;
  l.push_back(aln1);
  l.push_back(aln2);
  AlignmentHandle m=alg::MergePairwiseAlignments(l, ref);
  ConstSequenceList seqs=m.GetSequences();
  BOOST_CHECK_EQUAL(seqs[0].GetString(), "abcd---efghij---klmn");
  BOOST_CHECK_EQUAL(seqs[1].GetString(), "abcdxyzefghij---klmn");
  BOOST_CHECK_EQUAL(seqs[2].GetString(), "abcd---efghijxyzklmn");
  BOOST_CHECK_EQUAL(seqs[0].GetRole(), "R");
  BOOST_CHECK_EQUAL(seqs[1].GetRole(), "X");
  BOOST_CHECK_EQUAL(seqs[2].GetRole(), "Y");
}
  
BOOST_AUTO_TEST_CASE(merge_pairwise_alignments_two) 
{
  SequenceHandle ref=CreateSequence("REF", "abcdefghijklmn");
  SequenceHandle s1=CreateSequence("S1",   "abcd---efghijklmn");
  SequenceHandle s2=CreateSequence("S2",   "abcdxyzefghijklmn");
  AlignmentHandle aln1=CreateAlignment();
  aln1.AddSequence(s1);
  aln1.AddSequence(s2);
  
  SequenceHandle s3=CreateSequence("S1", "abcd-efghijklmn");
  SequenceHandle s4=CreateSequence("S2", "abcdyefghijklmn");
  
  AlignmentHandle aln2=CreateAlignment();
  aln2.AddSequence(s3);
  aln2.AddSequence(s4);
  AlignmentList l;
  l.push_back(aln1);
  l.push_back(aln2);
  AlignmentHandle m=alg::MergePairwiseAlignments(l, ref);
  ConstSequenceList seqs=m.GetSequences();
  BOOST_CHECK_EQUAL(seqs[0].GetString(), "abcd---efghijklmn");
  BOOST_CHECK_EQUAL(seqs[1].GetString(), "abcdxyzefghijklmn");
  BOOST_CHECK_EQUAL(seqs[2].GetString(), "abcd--yefghijklmn");
}


BOOST_AUTO_TEST_CASE(merge_pairwise_alignments_three)
{
  SequenceHandle ref=CreateSequence("REF", "abcdefghijk");
  SequenceHandle s1=CreateSequence("S1",   "--abcdefghijk");
  SequenceHandle s2=CreateSequence("S2",   "xyabcdefghijk");
  AlignmentHandle aln1=CreateAlignment();
  aln1.AddSequence(s1);
  aln1.AddSequence(s2);

  SequenceHandle s3=CreateSequence("S1", "-abcdefghijk");
  SequenceHandle s4=CreateSequence("S2", "zabcdefghijk");

  AlignmentHandle aln2=CreateAlignment();
  aln2.AddSequence(s3);
  aln2.AddSequence(s4);
  AlignmentList l;
  l.push_back(aln1);
  l.push_back(aln2);
  AlignmentHandle m=alg::MergePairwiseAlignments(l, ref);
  ConstSequenceList seqs=m.GetSequences();
  BOOST_CHECK_EQUAL(seqs[0].GetString(), "--abcdefghijk");
  BOOST_CHECK_EQUAL(seqs[1].GetString(), "xyabcdefghijk");
  BOOST_CHECK_EQUAL(seqs[2].GetString(), "-zabcdefghijk");
}

BOOST_AUTO_TEST_CASE(merge_pairwise_alignments_four)
{
  SequenceHandle ref=CreateSequence("REF", "abcdefghijk");
  SequenceHandle s1=CreateSequence("S1",   "abcdefghijk--");
  SequenceHandle s2=CreateSequence("S2",   "abcdefghijkxy");
  AlignmentHandle aln1=CreateAlignment();
  aln1.AddSequence(s1);
  aln1.AddSequence(s2);

  SequenceHandle s3=CreateSequence("S1", "-abcdefghijk");
  SequenceHandle s4=CreateSequence("S2", "zabcdefghijk");

  AlignmentHandle aln2=CreateAlignment();
  aln2.AddSequence(s3);
  aln2.AddSequence(s4);
  AlignmentList l;
  l.push_back(aln1);
  l.push_back(aln2);
  AlignmentHandle m=alg::MergePairwiseAlignments(l, ref);
  ConstSequenceList seqs=m.GetSequences();
  BOOST_CHECK_EQUAL(seqs[0].GetString(), "-abcdefghijk--");
  BOOST_CHECK_EQUAL(seqs[1].GetString(), "-abcdefghijkxy");
  BOOST_CHECK_EQUAL(seqs[2].GetString(), "zabcdefghijk--");
}

BOOST_AUTO_TEST_CASE(merge_pairwise_alignments_five)
{
  SequenceHandle ref=CreateSequence("REF", "abcdefghijk");
  ref.SetOffset(4);
  SequenceHandle s1=CreateSequence("S1",   "abcdefghijk--");
  SequenceHandle s2=CreateSequence("S2",   "abcdefghijkxy");
  s1.SetOffset(4);
  AlignmentHandle aln1=CreateAlignment();
  aln1.AddSequence(s1); 
  aln1.AddSequence(s2);

  SequenceHandle s3=CreateSequence("S1", "-abcdefghijk");
  SequenceHandle s4=CreateSequence("S2", "zabcdefghijk");
  s3.SetOffset(4);
  
  AlignmentHandle aln2=CreateAlignment();
  aln2.AddSequence(s3);
  aln2.AddSequence(s4);
  AlignmentList l;
  l.push_back(aln1);
  l.push_back(aln2);
  AlignmentHandle m=alg::MergePairwiseAlignments(l, ref);
  ConstSequenceList seqs=m.GetSequences();  
  
  BOOST_CHECK_EQUAL(seqs[0].GetString(), "-abcdefghijk--");
  BOOST_CHECK_EQUAL(seqs[0].GetOffset(), 4);
  BOOST_CHECK_EQUAL(seqs[1].GetString(), "-abcdefghijkxy");
  BOOST_CHECK_EQUAL(seqs[2].GetString(), "zabcdefghijk--");
}

BOOST_AUTO_TEST_CASE(merge_pairwise_alignments_six)
{
  SequenceHandle ref=CreateSequence("REF", "abcdefghijk");
  ref.SetOffset(4);
  SequenceHandle s1=CreateSequence("S1",   "abcdefghijk--");
  SequenceHandle s2=CreateSequence("S2",   "abcdefghijkxy");
  s1.SetOffset(5);
  AlignmentHandle aln1=CreateAlignment();
  aln1.AddSequence(s1); 
  aln1.AddSequence(s2);

  SequenceHandle s3=CreateSequence("S1", "-abcdefghijk");
  SequenceHandle s4=CreateSequence("S2", "zabcdefghijk");
  s3.SetOffset(4);
  
  AlignmentHandle aln2=CreateAlignment();
  aln2.AddSequence(s3);
  aln2.AddSequence(s4);
  AlignmentList l;
  l.push_back(aln1);
  l.push_back(aln2);
  BOOST_CHECK_THROW (alg::MergePairwiseAlignments(l, ref),IntegrityError);
  
}


BOOST_AUTO_TEST_SUITE_END();
