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

#include <ost/mol/mol.hh>

#include <ost/invalid_handle.hh>
#include <ost/integrity_error.hh>

#include <ost/mol/mol.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/aligned_region.hh>
#include <ost/seq/sequence_handle.hh>
#include <ost/seq/sequence_list.hh>
#include <ost/seq/impl/sequence_list_impl.hh>
#include <ost/seq/impl/sequence_impl.hh>
#include <ost/seq/sequence_list.hh>

using namespace ost;
using namespace ost::seq;
using namespace ost::mol;


BOOST_AUTO_TEST_SUITE( aligned_region );

BOOST_AUTO_TEST_CASE(aligned_region_triv)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","ab-cdef");
  SequenceHandle seq2 = CreateSequence("S2","ghijk-l");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  AlignedRegion ar = AlignedRegion(a,0,2,1);
  BOOST_CHECK(ar.GetAlignmentHandle()==a);
  BOOST_CHECK_EQUAL(ar.GetStart(),0);
  BOOST_CHECK_EQUAL(ar.GetEnd(),2);
  BOOST_CHECK_EQUAL(ar.GetLength(),2);
  BOOST_CHECK_EQUAL(ar.GetMaster(),1);
  ar.SetMaster(0);
  BOOST_CHECK_EQUAL(ar.GetMaster(),0);
  ar.SetMaster(-1);
  BOOST_CHECK_EQUAL(ar.GetMaster(),-1);
  BOOST_CHECK_THROW(ar.SetMaster(2),IntegrityError);
  BOOST_CHECK_THROW(ar.SetMaster(-2),IntegrityError);

  AlignedRegion ar_same = a.MakeRegion(0,2);
  BOOST_CHECK(ar==ar_same);
}

BOOST_AUTO_TEST_CASE(aligned_region_delete)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","ab-cdef");
  SequenceHandle seq2 = CreateSequence("S2","ghijk-l");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  AlignedRegion ar = a.MakeRegion(1,2);
  ar.Delete();
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"acdef");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"gjk-l");

  a = CreateAlignment();
  seq1 = CreateSequence("S1","ab-cdef");
  seq2 = CreateSequence("S2","ghijk-l");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  ar = a.MakeRegion(0,7);
  ar.Delete();
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"");

  a = CreateAlignment();
  seq1 = CreateSequence("S1","ab-cdef");
  seq2 = CreateSequence("S2","ghijk-l");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  ar = a.MakeRegion(0,0);
  ar.Delete();
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"ab-cdef");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"ghijk-l");
}

BOOST_AUTO_TEST_CASE(aligned_region_replace)
{
  AlignmentHandle a1 = CreateAlignment();
  SequenceHandle seq11 = CreateSequence("S1","ab-cdef");
  SequenceHandle seq12 = CreateSequence("S2","ghijk-l");
  a1.AddSequence(seq11);
  a1.AddSequence(seq12);
  AlignedRegion ar1 = a1.MakeRegion(2,2);
  AlignmentHandle a2 = CreateAlignment();
  SequenceHandle seq21 = CreateSequence("S1","zzzzzzz");
  SequenceHandle seq22 = CreateSequence("S2","ee-----");
  a2.AddSequence(seq21);
  a2.AddSequence(seq22);
  AlignedRegion ar2 = a2.MakeRegion(2,2);
  ar1.Replace(ar2);
  BOOST_CHECK_EQUAL(a1.GetSequence(0).GetString(),"abzzdef");
  BOOST_CHECK_EQUAL(a1.GetSequence(1).GetString(),"gh--k-l");

  a1 = CreateAlignment();
  seq11 = CreateSequence("S1","ab-cdef");
  seq12 = CreateSequence("S2","ghijk-l");
  a1.AddSequence(seq11);
  a1.AddSequence(seq12);
  ar1 = a1.MakeRegion(1,2);
  a2 = CreateAlignment();
  seq21 = CreateSequence("S1","zzzzzzz");
  seq22 = CreateSequence("S2","ee-----");
  a2.AddSequence(seq21);
  a2.AddSequence(seq22);
  ar2 = a2.MakeRegion(1,0);
  ar1.Replace(ar2);
  BOOST_CHECK_EQUAL(a1.GetSequence(0).GetString(),"acdef");
  BOOST_CHECK_EQUAL(a1.GetSequence(1).GetString(),"gjk-l");

  a1 = CreateAlignment();
  seq11 = CreateSequence("S1","ab-cdef");
  seq12 = CreateSequence("S2","ghijk-l");
  a1.AddSequence(seq11);
  a1.AddSequence(seq12);
  ar1 = a1.MakeRegion(1,2);
  a2 = CreateAlignment();
  seq21 = CreateSequence("S1","zzzzzzz");
  seq22 = CreateSequence("S2","ee-----");
  a2.AddSequence(seq21);
  a2.AddSequence(seq22);
  ar2 = a2.MakeRegion(1,3);
  ar1.Replace(ar2);
  BOOST_CHECK_EQUAL(a1.GetSequence(0).GetString(),"azzzcdef");
  BOOST_CHECK_EQUAL(a1.GetSequence(1).GetString(),"ge--jk-l");
}

BOOST_AUTO_TEST_CASE(aligned_region_shift)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","ab-c");
  SequenceHandle seq2 = CreateSequence("S2","ghij");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  AlignedRegion ar = a.MakeRegion(1,2);

  BOOST_CHECK_NO_THROW(ar.ShiftLeft(1));
  BOOST_CHECK_EQUAL(ar.GetStart(),0);
  BOOST_CHECK_EQUAL(ar.GetEnd(),2);
  BOOST_CHECK_EQUAL(ar.GetLength(),2);
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"b-ac");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"higj");
  BOOST_CHECK_NO_THROW(ar.ShiftLeft(0));
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"b-ac");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"higj");

  BOOST_CHECK_THROW(ar.ShiftLeft(1),IntegrityError);

  BOOST_CHECK_NO_THROW(ar.ShiftRight(2));
  BOOST_CHECK_EQUAL(ar.GetStart(),2);
  BOOST_CHECK_EQUAL(ar.GetEnd(),4);
  BOOST_CHECK_EQUAL(ar.GetLength(),2);
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"acb-");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"gjhi");

  BOOST_CHECK_NO_THROW(ar.ShiftRight(0));
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"acb-");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"gjhi");

  BOOST_CHECK_THROW(ar.ShiftRight(1),IntegrityError);
}

BOOST_AUTO_TEST_SUITE_END();
