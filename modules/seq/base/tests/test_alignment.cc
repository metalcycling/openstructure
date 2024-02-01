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

#include <ost/invalid_handle.hh>
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


struct Fixture {
  Fixture() {
    eh = CreateEntity();
    eh.SetName("TestEntity");
    XCSEditor e=eh.EditXCS();
    chain = e.InsertChain("A");
    res1 = e.AppendResidue(chain, "ARG");
    e.InsertAtom(res1, "CA",geom::Vec3(1,0,0));
    res2 = e.AppendResidue(chain, "ARG");
    e.InsertAtom(res2, "CA",geom::Vec3(0,1,0));
    res3 = e.AppendResidue(chain, "ARG");
    e.InsertAtom(res3, "CA",geom::Vec3(0,0,1));
  }
  EntityHandle eh;
  ChainHandle  chain;
  ResidueHandle res1;
  ResidueHandle res2;
  ResidueHandle res3;
};

BOOST_AUTO_TEST_SUITE( alignment );

BOOST_AUTO_TEST_CASE(alignment_triv)
{
  AlignmentHandle a;
  BOOST_CHECK_THROW(a.AddSequence(CreateSequence("S1","-asdf-")),InvalidHandle);
  a = CreateAlignment();
  BOOST_CHECK_NO_THROW(a.AddSequence(CreateSequence("S1","-asdf-")));
  BOOST_CHECK_THROW(a.AddSequence(SequenceHandle()),InvalidSequence);
  BOOST_CHECK_THROW(a.AddSequence(CreateSequence("S2","-asdf-f")),
                    Error);
  SequenceList list = CreateSequenceList();
  list.AddSequence(CreateSequence("S1", "-asdf-"));
  list.AddSequence(CreateSequence("S2", "fasdfas"));
  BOOST_CHECK_THROW(AlignmentFromSequenceList(list), Error);
  list = CreateSequenceList();
  list.AddSequence(CreateSequence("S1", "-asdf-"));
  list.AddSequence(CreateSequence("S2", "fasdfa"));
  BOOST_CHECK_NO_THROW(AlignmentFromSequenceList(list));
}

BOOST_AUTO_TEST_CASE(aln_throw_invalid)
{
  AlignmentHandle s;
  BOOST_CHECK_THROW(CheckHandleValidity(s), InvalidHandle);
  s=CreateAlignment();
  BOOST_CHECK_NO_THROW(CheckHandleValidity(s));
}


BOOST_AUTO_TEST_CASE(alignment_count_and_length)
{
  AlignmentHandle a = CreateAlignment();
  BOOST_CHECK_EQUAL(a.GetCount(),0);
  BOOST_CHECK_EQUAL(a.GetLength(),0);
  a.AddSequence(CreateSequence("S1","asdf"));
  BOOST_CHECK_EQUAL(a.GetCount(),1);
  BOOST_CHECK_EQUAL(a.GetLength(),4);
  a.AddSequence(CreateSequence("S2","qwer"));
  BOOST_CHECK_EQUAL(a.GetCount(),2);
  BOOST_CHECK_EQUAL(a.GetLength(),4);
}

BOOST_AUTO_TEST_CASE(alignment_onelettercode)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","-a-sdf");
  SequenceHandle seq2 = CreateSequence("S2","q--wer");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(0,0),'-');
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(0,0),seq1.GetOneLetterCode(0));
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(1,0),'q');
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(1,0),seq2.GetOneLetterCode(0));
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(0,3),'s');
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(0,3),seq1.GetOneLetterCode(3));
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(1,3),'w');
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(1,3),seq2.GetOneLetterCode(3));
  BOOST_CHECK_THROW(a.GetOneLetterCode(0,-1),std::out_of_range);
  BOOST_CHECK_THROW(a.GetOneLetterCode(1,6),std::out_of_range);
  BOOST_CHECK_THROW(a.GetOneLetterCode(2,0),std::out_of_range);
  BOOST_CHECK_THROW(a.GetOneLetterCode(-1,0), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(alignment_pos_and_index)
{
  AlignmentHandle a = CreateAlignment();
  a.AddSequence(CreateSequence("S1","-a-sdf"));
  a.AddSequence(CreateSequence("S2","q--wer"));
  BOOST_CHECK_EQUAL(a.GetPos(0,1), 3);
  BOOST_CHECK_EQUAL(a.GetResidueIndex(0,3), 1);
  BOOST_CHECK_EQUAL(a.GetPos(1,2), 4);
  BOOST_CHECK_EQUAL(a.GetResidueIndex(1,4), 2);
  BOOST_CHECK_THROW(a.GetPos(0,-1),Error);
  BOOST_CHECK_THROW(a.GetPos(1,4),Error);
  BOOST_CHECK_THROW(a.GetPos(2,0),std::out_of_range);
  BOOST_CHECK_THROW(a.GetPos(-1,0),std::out_of_range);
}

BOOST_AUTO_TEST_CASE(alignment_getseq)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","-a-sdf");
  SequenceHandle seq2 = CreateSequence("S2","q--wer");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  BOOST_CHECK_EQUAL(a.GetSequence(0), seq1);
  BOOST_CHECK_EQUAL(a.GetSequence(1), seq2);
  BOOST_CHECK_THROW(a.GetSequence(-1), std::out_of_range);
  BOOST_CHECK_THROW(a.GetSequence(2), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(alignment_attach_view)
{
  Fixture f;
  EntityView v = f.eh.CreateFullView();
  EntityView v2 = f.eh.CreateFullView();

  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","r-r");
  SequenceHandle seq2 = CreateSequence("S2","aa-");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  a.AttachView(0,v);
  a.AttachView(1,v2);
  BOOST_CHECK_EQUAL(a.GetResidue(0,0), f.res1);
  BOOST_CHECK_EQUAL(a.GetResidue(0,2), f.res3);
  BOOST_CHECK_EQUAL(a.GetResidue(1,0), f.res1);
  BOOST_CHECK_EQUAL(a.GetResidue(1,1), f.res2);
  BOOST_CHECK_NO_THROW(a.AttachView(0,EntityView()));
  BOOST_CHECK_EQUAL(a.GetResidue(0,1), EntityView());
  BOOST_CHECK_EQUAL(a.GetResidue(0,3), EntityView());
  BOOST_CHECK_THROW(a.GetResidue(-1,0),std::out_of_range);
}

BOOST_AUTO_TEST_CASE(alignment_cut)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","r-r");
  SequenceHandle seq2 = CreateSequence("S2","aa-");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  a.Cut(1,2);
  BOOST_CHECK_EQUAL(a.GetLength(),2);
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(0,0),'r');
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(0,1),'r');
  BOOST_CHECK_THROW(a.GetOneLetterCode(0,2),std::out_of_range);
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(1,0),'a');
  BOOST_CHECK_EQUAL(a.GetOneLetterCode(1,1),'-');
  BOOST_CHECK_THROW(a.GetOneLetterCode(1,2),std::out_of_range);
}

BOOST_AUTO_TEST_CASE(alignment_shift_region)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","abcdef");
  SequenceHandle seq2 = CreateSequence("S2","ghijkl");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  a.ShiftRegion(2, 4, 2);
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"abefcd");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"ghklij");

  a.ShiftRegion(2, 2, 2);//Shift nothing
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"abefcd");

  a.ShiftRegion(2,3,1,1);
  BOOST_CHECK_EQUAL(a.GetSequence(0).GetString(),"abefcd");
  BOOST_CHECK_EQUAL(a.GetSequence(1).GetString(),"ghlkij");

  BOOST_CHECK_THROW(a.ShiftRegion(0,2,1,2),std::out_of_range);
  BOOST_CHECK_THROW(a.ShiftRegion(0,2,1,-2),std::out_of_range);

  BOOST_CHECK_THROW(a.ShiftRegion(0,2,-1),std::out_of_range);
  BOOST_CHECK_THROW(a.ShiftRegion(0,5,2),std::out_of_range);
}

BOOST_AUTO_TEST_CASE(alignment_aligned_region)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","abcdef");
  SequenceHandle seq2 = CreateSequence("S2","ghijkl");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  AlignedRegion ar = a.MakeRegion(1,2);
  BOOST_CHECK_EQUAL(ar.GetLength(),2);
  BOOST_CHECK_EQUAL(ar.GetStart(),1);
  BOOST_CHECK_EQUAL(ar.GetEnd(),3);

  BOOST_CHECK_THROW(ar = a.MakeRegion(-1,2,0),std::out_of_range);
  BOOST_CHECK_THROW(ar = a.MakeRegion(0,7,0),std::out_of_range);
  BOOST_CHECK_NO_THROW(ar = a.MakeRegion(2,4,0));
  BOOST_CHECK_NO_THROW(ar = a.MakeRegion(3,2,0));
  BOOST_CHECK_NO_THROW(ar = a.MakeRegion(3,3,0));
}

BOOST_AUTO_TEST_SUITE_END();
