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
  Author: Gerardo Tauriello, Gabriel Studer
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <iostream>

#include <ost/seq/profile_handle.hh>
#include <ost/io/seq/hhm_io_handler.hh>
#include <ost/io/seq/pssm_io_handler.hh>

using namespace ost::seq;

BOOST_AUTO_TEST_SUITE( profile );

namespace {

void CheckNullModel(const ProfileColumn& nm) {
  BOOST_CHECK_GE(nm.GetEntropy(), Real(0));
  Real sum_freq = 0;
  const Real* freq_ptr = nm.freqs_begin();
  for (uint i = 0; i < 20; ++i) {
    BOOST_CHECK_GE(freq_ptr[i], Real(0));
    sum_freq += freq_ptr[i];
  }
  BOOST_CHECK_CLOSE(sum_freq, Real(1), Real(1e-2));
}

} // anon ns

BOOST_AUTO_TEST_CASE(profile_column)
{
  const char * olc = "ACDEFGHIKLMNPQRSTVWY";
  
  // check empty column
  ProfileColumn pc;
  BOOST_CHECK_EQUAL(pc.GetEntropy(), Real(0));
  Real* freq_ptr = pc.freqs_begin();
  for (uint i = 0; i < 20; ++i) {
    BOOST_CHECK_EQUAL(freq_ptr[i], Real(0));
    BOOST_CHECK_EQUAL(pc.GetFreq(olc[i]), Real(0));
    // check overwrite
    const Real tst = Real(i)/Real(20);
    pc.SetFreq(olc[i], tst);
    BOOST_CHECK_EQUAL(pc.GetFreq(olc[i]), tst);
  }

  // check null models
  ProfileColumn nm = ProfileColumn::BLOSUMNullModel();
  CheckNullModel(nm);
  CheckNullModel(ProfileColumn::HHblitsNullModel());

  // check assign/copy
  ProfileColumn nm2(nm);
  BOOST_CHECK(nm == nm2);
  freq_ptr = nm.freqs_begin();
  Real* freq_ptr2 = nm2.freqs_begin();
  for (uint i = 0; i < 20; ++i) {
    BOOST_CHECK_EQUAL(freq_ptr[i], freq_ptr2[i]);
  }
  pc = nm;
  BOOST_CHECK(nm == pc);
}

BOOST_AUTO_TEST_CASE(profile_handle)
{
  // generate dummy
  ProfileHandle prof;
  BOOST_CHECK_EQUAL(prof.size(), size_t(0));
  BOOST_CHECK_EQUAL(prof.GetSequence(), "");
  CheckNullModel(prof.GetNullModel());
  ProfileColumn pc0 = ProfileColumn::BLOSUMNullModel();
  prof.SetNullModel(pc0);
  ProfileColumn pc1, pc2, pc3, pc4;
  pc1.SetFreq('H', 0.74);
  pc1.SetFreq('E', 0.01);
  pc1.SetFreq('L', 0.25);
  pc2.SetFreq('L', 1);
  pc3.SetFreq('Y', 0.5);
  pc3.SetFreq('E', 0.5);
  pc4.SetFreq('A', 0.5);
  pc4.SetFreq('H', 0.5);
  prof.AddColumn(pc0);
  prof.AddColumn(pc1);
  prof.AddColumn(pc2);
  BOOST_CHECK_EQUAL(prof.size(), prof.GetSequence().length());
  prof.AddColumn(pc3);
  prof.AddColumn(pc4);
  prof.SetSequence("RRSPP");
  BOOST_CHECK_EQUAL(prof.size(), prof.GetSequence().length());

  // check invalid sequence settings
  BOOST_CHECK_THROW(prof.SetSequence("RRSP"), ost::Error);
  BOOST_CHECK_THROW(prof.SetSequence("RRSPPT"), ost::Error);

  // check components
  BOOST_CHECK_EQUAL(prof.size(), size_t(5));
  BOOST_CHECK(prof.GetNullModel() == pc0);
  BOOST_CHECK_EQUAL(prof.GetSequence(), "RRSPP");

  // check columns with different methods of getting the columns
  BOOST_CHECK(prof[0] == pc0);
  BOOST_CHECK(*(prof.columns_begin()+1) == pc1);
  BOOST_CHECK(prof.GetColumns()[2] == pc2);
  BOOST_CHECK(prof.at(3) == pc3);
  BOOST_CHECK(*(prof.columns_end()-1) == pc4);

  // check entropies
  Real entropy1 = -0.74*log(0.74) - 0.01*log(0.01) - 0.25*log(0.25);
  Real entropy2 = 0;
  Real entropy3 = -log(0.5);
  Real avg_entropy = (pc0.GetEntropy() + entropy1 + entropy2 + 2*entropy3)/5;
  BOOST_CHECK_CLOSE(prof[1].GetEntropy(), entropy1, Real(1e-3));
  BOOST_CHECK_CLOSE(prof[2].GetEntropy(), entropy2, Real(1e-3));
  BOOST_CHECK_CLOSE(prof[3].GetEntropy(), entropy3, Real(1e-3));
  BOOST_CHECK_CLOSE(prof.GetAverageEntropy(), avg_entropy, Real(1e-3));

  // check scoring (scores precomputed to make sure they don't randomly change!)
  ProfileHandlePtr pp = prof.Extract(1, 4);
  BOOST_CHECK_CLOSE(prof.GetAverageScore(*pp, 0), -1.81746, Real(1e-3));
  BOOST_CHECK_CLOSE(prof.GetAverageScore(*pp, 1), 2.8012, Real(1e-3));
  BOOST_CHECK_CLOSE(prof.GetAverageScore(*pp, 2), -4.12038, Real(1e-3));
  BOOST_CHECK_THROW(prof.GetAverageScore(*pp, -1), ost::Error);
  BOOST_CHECK_THROW(prof.GetAverageScore(*pp, 3), ost::Error);

  // extract parts
  avg_entropy = (entropy1 + entropy2 + entropy3)/3;
  BOOST_CHECK_EQUAL(pp->size(), size_t(3));
  BOOST_CHECK_EQUAL(pp->GetSequence(), "RSP");
  BOOST_CHECK_CLOSE(pp->GetAverageEntropy(), avg_entropy, Real(1e-3));

  pp = prof.Extract(0, 3);
  avg_entropy = (pc0.GetEntropy() + entropy1 + entropy2)/3;
  BOOST_CHECK_EQUAL(pp->size(), size_t(3));
  BOOST_CHECK_EQUAL(pp->GetSequence(), "RRS");
  BOOST_CHECK_CLOSE(pp->GetAverageEntropy(), avg_entropy, Real(1e-3));

  pp = prof.Extract(3, 5);
  BOOST_CHECK_EQUAL(pp->size(), size_t(2));
  BOOST_CHECK_EQUAL(pp->GetSequence(), "PP");
  BOOST_CHECK_CLOSE(pp->GetAverageEntropy(), entropy3, Real(1e-3));

  BOOST_CHECK_THROW(prof.Extract(3, 3), ost::Error);
  BOOST_CHECK_THROW(prof.Extract(-1, 3), ost::Error);
  BOOST_CHECK_THROW(prof.Extract(3, 6), ost::Error);
}


BOOST_AUTO_TEST_SUITE_END();
