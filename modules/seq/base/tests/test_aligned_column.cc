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
#include <ost/seq/aligned_column.hh>
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

BOOST_AUTO_TEST_SUITE( aligned_column );

BOOST_AUTO_TEST_CASE(aligned_column_triv)
{
  AlignedColumn ac =AlignedColumn();
  BOOST_CHECK_THROW(ac[0], InvalidHandle);

  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","ab-cdef");
  SequenceHandle seq2 = CreateSequence("S2","ghijk-l");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  AlignedRegion ar = AlignedRegion(a,1,3,1);
  ac = ar[0];

  BOOST_CHECK_EQUAL(ac.GetIndex(), 1);
  ac = ar[1];
  BOOST_CHECK_EQUAL(ac.GetIndex(), 2);
  BOOST_CHECK_THROW(ar[2], std::out_of_range);
  BOOST_CHECK_THROW(ar[-1], std::out_of_range);

  AlignedColumn ac_test = AlignedColumn(a,2);
  char c = ac[0];
  BOOST_CHECK_EQUAL(c, ac_test[0]);
  BOOST_CHECK_EQUAL(c, '-');
  c = ac[1];
  BOOST_CHECK_EQUAL(c, ac_test[1]);
  BOOST_CHECK_EQUAL(c, 'i');
}

BOOST_AUTO_TEST_CASE(aligned_get_residue)
{
  Fixture f;
  EntityView v1 = f.eh.CreateFullView();
  EntityView v2 = f.eh.CreateFullView();

  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","ab-cdef");
  SequenceHandle seq2 = CreateSequence("S2","ghijk-l");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  a.AttachView(0,v1);
  a.AttachView(1,v2);
  AlignedRegion ar = AlignedRegion(a,1,3,1);
  AlignedColumn ac = ar[0];
  BOOST_CHECK_EQUAL(f.res2,ac.GetResidue(0));
  ac = ar[1];
  BOOST_CHECK_EQUAL(f.res2,ac.GetResidue(1));
  BOOST_CHECK_THROW(ac.GetResidue(2),std::out_of_range);
  BOOST_CHECK_THROW(ac.GetResidue(-1),std::out_of_range);
}

BOOST_AUTO_TEST_CASE(aligned_get_residue_index)
{
  AlignmentHandle a = CreateAlignment();
  SequenceHandle seq1 = CreateSequence("S1","ab-cdef");
  SequenceHandle seq2 = CreateSequence("S2","ghijk-l");
  a.AddSequence(seq1);
  a.AddSequence(seq2);
  AlignedRegion ar = AlignedRegion(a,1,3,1);
  AlignedColumn ac = ar[0];
  BOOST_CHECK_EQUAL(ac.GetResidueIndex(0),1);
  BOOST_CHECK_EQUAL(ac.GetResidueIndex(1),1);

  ac = ar[1];
  BOOST_CHECK_THROW(ac.GetResidueIndex(0),Error);
  BOOST_CHECK_EQUAL(ac.GetResidueIndex(1),2);
  BOOST_CHECK_THROW(ac.GetResidueIndex(2),std::out_of_range);
  BOOST_CHECK_THROW(ac.GetResidueIndex(-1),std::out_of_range);
}

BOOST_AUTO_TEST_SUITE_END();
