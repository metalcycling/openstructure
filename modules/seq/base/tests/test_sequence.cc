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

#include <ost/seq/sequence_handle.hh>
#include <ost/seq/invalid_sequence.hh>
#include <ost/seq/sequence_op.hh>

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

BOOST_AUTO_TEST_SUITE( seq );

BOOST_AUTO_TEST_CASE(seq_triv) 
{
  BOOST_CHECK_NO_THROW(CreateSequence("S1", "-afbcdefghijkLMNOPQRSTUV-"));
  BOOST_CHECK_THROW(CreateSequence("S1", "1"), InvalidSequence);
  BOOST_CHECK_THROW(CreateSequence("S1", "*"), InvalidSequence);
  BOOST_CHECK_THROW(CreateSequence("S1", " "), InvalidSequence);
  SequenceHandle s=CreateSequence("S1","-afbcdefghijkLMNOPQRSTUV-");
  BOOST_CHECK_EQUAL(s.GetString(),"-afbcdefghijkLMNOPQRSTUV-");
  BOOST_CHECK_NO_THROW(s.SetString("-afc--de-f"));
  BOOST_CHECK_EQUAL(s.GetString(),"-afc--de-f");
  BOOST_CHECK_THROW(s.SetString("1"), InvalidSequence);
}

BOOST_AUTO_TEST_CASE(match)
{
  BOOST_CHECK(Match("abcdefghijkl", "ABcDeFgHiJkL"));
  BOOST_CHECK(Match("abcxXxxxxjXl", "ABcDeFgHiJkL"));
  BOOST_CHECK(Match("ABcDeFgHiJkL", "ABcDeFXxiJxL"));
  BOOST_CHECK(!Match("abc", "abcd"));
  BOOST_CHECK(!Match("abc", "aby"));
}

BOOST_AUTO_TEST_CASE(seq_throw_invalid)
{
  SequenceHandle s;
  BOOST_CHECK_THROW(CheckHandleValidity(s), InvalidHandle);
  s=CreateSequence("A", "abcd");
  BOOST_CHECK_NO_THROW(CheckHandleValidity(s));
}

BOOST_AUTO_TEST_CASE(seq_length)
{
  SequenceHandle s=CreateSequence("S1", "abfcdadeaf");
  BOOST_CHECK_EQUAL(s.GetLength(),10);
  s=CreateSequence("S1", "-afc--de-f");
  BOOST_CHECK_EQUAL(s.GetLength(),10);
  s=CreateSequence("S1", "");
  BOOST_CHECK_EQUAL(s.GetLength(),0);
}

BOOST_AUTO_TEST_CASE(seq_string)
{
  SequenceHandle s=CreateSequence("S1", "abfcdadeaf");
  BOOST_CHECK_EQUAL(s.GetString(),"abfcdadeaf");
  BOOST_CHECK_EQUAL(s.GetGaplessString(),"abfcdadeaf");
  s=CreateSequence("S1", "-afc--de-f");
  BOOST_CHECK_EQUAL(s.GetString(),"-afc--de-f");
  BOOST_CHECK_EQUAL(s.GetGaplessString(),"afcdef");
  s=CreateSequence("S1", "");
  BOOST_CHECK_EQUAL(s.GetString(),"");
  BOOST_CHECK_EQUAL(s.GetGaplessString(),"");
}

BOOST_AUTO_TEST_CASE(seq_normalise)
{
  SequenceHandle s=CreateSequence("S1", "abfcdadeaf");
  BOOST_CHECK_EQUAL(s.GetString(),"abfcdadeaf");
  s.Normalise();
  BOOST_CHECK_EQUAL(s.GetString(),"ABFCDADEAF");
  s=CreateSequence("S1", ".afc..de.f");
  BOOST_CHECK_EQUAL(s.GetString(),".afc..de.f");
  s.Normalise();
  BOOST_CHECK_EQUAL(s.GetString(),"AFCDEF");
  s=CreateSequence("S1", "-afc--de-f");
  BOOST_CHECK_EQUAL(s.GetString(),"-afc--de-f");
  s.Normalise();
  BOOST_CHECK_EQUAL(s.GetString(),"AFCDEF");
  s=CreateSequence("S1", ".afc-.de-f");
  BOOST_CHECK_EQUAL(s.GetString(),".afc-.de-f");
  s.Normalise();
  BOOST_CHECK_EQUAL(s.GetString(),"AFCDEF");
}

BOOST_AUTO_TEST_CASE(seq_onelettercode)
{
  SequenceHandle s=CreateSequence("S1", "abfcdadeaf");
  BOOST_CHECK_EQUAL(s.GetOneLetterCode(0),'a');
  BOOST_CHECK_EQUAL(s.GetOneLetterCode(3),'c');
  BOOST_CHECK_EQUAL(s.GetOneLetterCode(9),'f');
  BOOST_CHECK_THROW(s.GetOneLetterCode(-1),std::out_of_range);
  BOOST_CHECK_THROW(s.GetOneLetterCode(10),std::out_of_range);
  s=CreateSequence("S1", "-afc--de-f");
  BOOST_CHECK_EQUAL(s.GetOneLetterCode(0),'-');
  BOOST_CHECK_EQUAL(s.GetOneLetterCode(1),'a');
  BOOST_CHECK_EQUAL(s.GetOneLetterCode(9),'f');
  BOOST_CHECK_THROW(s.GetOneLetterCode(-1),std::out_of_range);
  BOOST_CHECK_THROW(s.GetOneLetterCode(10),std::out_of_range);
  s=CreateSequence("S1", "");
  BOOST_CHECK_THROW(s.GetOneLetterCode(0), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(seq_getnum) 
{
  SequenceHandle s=CreateSequence("S1", "-afc--de-f");
  BOOST_CHECK_THROW(s.GetResidueIndex(0), Error);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(1),1-1);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(2),2-1);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(3),3-1);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(6),4-1);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(7),5-1);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(9),6-1);
  BOOST_CHECK_THROW(s.GetResidueIndex(-1), std::out_of_range);
  BOOST_CHECK_THROW(s.GetResidueIndex(10), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(seq_getpos) 
{
  SequenceHandle s=CreateSequence("S1", "-afc--de-f");

  BOOST_CHECK(s.GetPos(0)==1);
  BOOST_CHECK(s.GetPos(1)==2);
  BOOST_CHECK(s.GetPos(2)==3);
  BOOST_CHECK(s.GetPos(3)==6);
  BOOST_CHECK(s.GetPos(4)==7);
  BOOST_CHECK(s.GetPos(5)==9);
}

BOOST_AUTO_TEST_CASE(seq_offset)
{
  SequenceHandle s=CreateSequence("S1", "-afc--de-f");
  SequenceHandle s2=CreateSequence("S1", "-afc--de-f");
  s.SetOffset(2);
  BOOST_CHECK_THROW(s.GetPos(-1), Error);
  BOOST_CHECK_THROW(s.GetPos(0), Error);
  BOOST_CHECK_THROW(s.GetPos(1), Error);
  BOOST_CHECK_EQUAL(s.GetPos(2),1);
  BOOST_CHECK_EQUAL(s.GetPos(2),s2.GetPos(2-2));
  BOOST_CHECK_EQUAL(s.GetPos(5),6);
  BOOST_CHECK_EQUAL(s.GetPos(5),s2.GetPos(5-2));
  BOOST_CHECK_EQUAL(s.GetPos(7),9);
  BOOST_CHECK_EQUAL(s.GetPos(7),s2.GetPos(7-2));
  BOOST_CHECK_THROW(s.GetPos(8), Error);
  BOOST_CHECK_THROW(s.GetResidueIndex(0), Error);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(1),2);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(1),s2.GetResidueIndex(1)+2);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(6), 5);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(6),s2.GetResidueIndex(6)+2);
  BOOST_CHECK_THROW(s.GetResidueIndex(10), std::out_of_range);
  BOOST_CHECK_THROW(s.GetResidueIndex(-1), std::out_of_range);

  s.SetOffset(-1);
  BOOST_CHECK_THROW(s.GetPos(-2), Error);
  BOOST_CHECK_EQUAL(s.GetPos(-1), 1);
  BOOST_CHECK_EQUAL(s.GetPos(0), 2);
  BOOST_CHECK_EQUAL(s.GetPos(0),s2.GetPos(0+1));
  BOOST_CHECK_EQUAL(s.GetPos(2), 6);
  BOOST_CHECK_EQUAL(s.GetPos(2),s2.GetPos(2+1));
  BOOST_CHECK_EQUAL(s.GetPos(4),9);
  BOOST_CHECK_EQUAL(s.GetPos(4),s2.GetPos(4+1));
  BOOST_CHECK_THROW(s.GetPos(5), Error);
  BOOST_CHECK_THROW(s.GetResidueIndex(0), Error);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(1),-1);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(1),s2.GetResidueIndex(1)-1);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(6), 2);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(6),s2.GetResidueIndex(6)-1);
  BOOST_CHECK_THROW(s.GetResidueIndex(10), std::out_of_range);
}

BOOST_AUTO_TEST_CASE(seq_gaps)
{
  SequenceHandle s=CreateSequence("S1", "");
  BOOST_CHECK_EQUAL(s.GetFirstNonGap(),-1);
  BOOST_CHECK_EQUAL(s.GetLastNonGap(),-1);

  s=CreateSequence("S1", "afcdef");
  BOOST_CHECK_EQUAL(s.GetFirstNonGap(),0);
  BOOST_CHECK_EQUAL(s.GetLastNonGap(),5);

  s=CreateSequence("S1", "-afc--de-f-");
  BOOST_CHECK_EQUAL(s.GetFirstNonGap(),1);
  BOOST_CHECK_EQUAL(s.GetLastNonGap(),9);
}

BOOST_AUTO_TEST_CASE(seq_append_olc)
{
  SequenceHandle s=CreateSequence("S1", "");
  // check if the shift-table gets setup properly
  s.Append('-');
  s.Append('-');
  s.Append('a');
  s.Append('b');
  s.Append('-');
  s.Append('c');
  BOOST_CHECK_EQUAL(s.GetResidueIndex(2), 0);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(3), 1);
  BOOST_CHECK_EQUAL(s.GetResidueIndex(5), 2);
}

BOOST_AUTO_TEST_CASE(seq_attach_view)
{
  Fixture f;
  SequenceHandle s=CreateSequence("S1", "r-r");
  BOOST_CHECK_EQUAL(s.HasAttachedView(),false);
  BOOST_CHECK_EQUAL(s.GetAttachedView(),EntityView());
  BOOST_CHECK_EQUAL(s.GetResidue(0),ResidueHandle());

  EntityView v = f.eh.CreateFullView();
  s.AttachView(v);
  BOOST_CHECK_EQUAL(s.HasAttachedView(),true);
  BOOST_CHECK_EQUAL(s.GetAttachedView(),v);

  BOOST_CHECK_EQUAL(s.GetResidue(0),f.res1);
  BOOST_CHECK_EQUAL(s.GetResidue(s.GetPos(0)),f.res1);
  BOOST_CHECK_EQUAL(s.GetResidue(1),ResidueHandle());
  BOOST_CHECK_EQUAL(s.GetResidue(2),f.res3);
  BOOST_CHECK_EQUAL(s.GetResidue(s.GetPos(1)),f.res3);

  v = EntityView();
  s.AttachView(v);
  BOOST_CHECK_EQUAL(s.HasAttachedView(),false);
  BOOST_CHECK_EQUAL(s.GetAttachedView(),EntityView());
  BOOST_CHECK_EQUAL(s.GetResidue(0),ResidueHandle());
}

BOOST_AUTO_TEST_SUITE_END();
