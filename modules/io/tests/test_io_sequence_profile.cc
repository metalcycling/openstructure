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
#include <ost/io/seq/load.hh>
#include <ost/io/seq/hhm_io_handler.hh>
#include <ost/io/seq/pssm_io_handler.hh>
#include <ost/io/io_exception.hh>

using namespace ost::seq;
using namespace ost::io;

BOOST_AUTO_TEST_SUITE( io );

// helper to compare profiles element by element (with delta)
// delta given in % !
template <typename T, typename T2>
void compareProfiles(T& t, T2& t2, Real delta = 0) {
  BOOST_REQUIRE_EQUAL(t.size(), t2.size());
  BOOST_CHECK_EQUAL(t.GetSequence(), t2.GetSequence());
  BOOST_CHECK_CLOSE(t.GetAverageEntropy(), t2.GetAverageEntropy(), delta*20);
  const Real* f = t.GetNullModel().freqs_begin();
  const Real* f2 = t2.GetNullModel().freqs_begin();
  for (uint j = 0; j < 20; ++j) {
    BOOST_CHECK_CLOSE(f[j], f2[j], delta);
  }
  uint N = t.size();
  for (uint i = 0; i < N; ++i) {
    f = t[i].freqs_begin();
    f2 = t2[i].freqs_begin();
    for (uint j = 0; j < 20; ++j) {
      BOOST_CHECK_CLOSE(f[j], f2[j], delta);
    }
  }
}

BOOST_AUTO_TEST_CASE(hhm_loading)
{
  // check hhm factory
  ProfileIOHandlerFactoryBasePtr fac(new HhmIOHandlerFactory);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test.hhm"),     true);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test.hhm.gz"),  true);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test.pssm"),    false);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test.pssm.gz"), false);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test", "hhm"),  true);
  BOOST_CHECK_EQUAL(fac->ProvidesExport("test.hhm"),     false);
  BOOST_CHECK_EQUAL(fac->ProvidesExport("test.hhm.gz"),  false);
  BOOST_CHECK_EQUAL(fac->ProvidesExport("test", "hhm"),  false);
  ProfileIOHandlerPtr io_handler = fac->Create();

  // load file
  ProfileHandlePtr prof(new ProfileHandle);
  io_handler->Import(*prof, "testfiles/test_hmm.hhm");
  BOOST_CHECK_EQUAL(prof->size(), size_t(5));
  BOOST_CHECK_EQUAL(prof->GetSequence(), "ETGSP");

  // we first check, whether the null model has been read correctly
  // and then randomly pick a column to check

  ProfileColumn null_model = prof->GetNullModel();

  String olc = "ACDEFGHIKLMNPQRSTVWY";

  Real correct_aa_freqs[] = {static_cast<Real>(pow(2.0, -0.001*3706)),
                             static_cast<Real>(pow(2.0, -0.001*5728)),
                             static_cast<Real>(pow(2.0, -0.001*4211)),
                             static_cast<Real>(pow(2.0, -0.001*4064)),
                             static_cast<Real>(pow(2.0, -0.001*4839)),
                             static_cast<Real>(pow(2.0, -0.001*3729)),
                             static_cast<Real>(pow(2.0, -0.001*4763)),
                             static_cast<Real>(pow(2.0, -0.001*4308)),
                             static_cast<Real>(pow(2.0, -0.001*4069)),
                             static_cast<Real>(pow(2.0, -0.001*3323)),
                             static_cast<Real>(pow(2.0, -0.001*5509)),
                             static_cast<Real>(pow(2.0, -0.001*4640)),
                             static_cast<Real>(pow(2.0, -0.001*4464)),
                             static_cast<Real>(pow(2.0, -0.001*4937)),
                             static_cast<Real>(pow(2.0, -0.001*4285)),
                             static_cast<Real>(pow(2.0, -0.001*4423)),
                             static_cast<Real>(pow(2.0, -0.001*3815)),
                             static_cast<Real>(pow(2.0, -0.001*3783)),
                             static_cast<Real>(pow(2.0, -0.001*6325)),
                             static_cast<Real>(pow(2.0, -0.001*4665))};

  for (int i = 0; i < 20; ++i) {
    BOOST_CHECK_CLOSE(null_model.GetFreq(olc[i]), correct_aa_freqs[i],
                      Real(1e-5));
  }

  // let's check, whether a particular column has been read correctly  
  ProfileColumn col = (*prof)[2];

  memset(correct_aa_freqs, 0, 20*sizeof(Real));
  correct_aa_freqs[0] = static_cast<Real>(pow(2.0, -0.001*3676));
  correct_aa_freqs[1] = static_cast<Real>(pow(2.0, -0.001*2597));
  correct_aa_freqs[5] = static_cast<Real>(pow(2.0, -0.001*5359));
  correct_aa_freqs[6] = static_cast<Real>(pow(2.0, -0.001*3275));
  correct_aa_freqs[11] = static_cast<Real>(pow(2.0, -0.001*3292));
  correct_aa_freqs[12] = static_cast<Real>(pow(2.0, -0.001*5077));
  correct_aa_freqs[13] = static_cast<Real>(pow(2.0, -0.001*3826));
  correct_aa_freqs[15] = static_cast<Real>(pow(2.0, -0.001*2409));
  correct_aa_freqs[16] = static_cast<Real>(pow(2.0, -0.001*3733));
  correct_aa_freqs[17] = static_cast<Real>(pow(2.0, -0.001*4503));
  correct_aa_freqs[19] = static_cast<Real>(pow(2.0, -0.001*3070));

  for(int i = 0; i < 20; ++i){
    BOOST_CHECK_CLOSE(col.GetFreq(olc[i]), correct_aa_freqs[i], Real(1e-5));
  }

  // check, whether the entropy of the column above gets correctly calculated
  Real entropy = 0.0;
  entropy -= correct_aa_freqs[0] * log(correct_aa_freqs[0]);
  entropy -= correct_aa_freqs[1] * log(correct_aa_freqs[1]);
  entropy -= correct_aa_freqs[5] * log(correct_aa_freqs[5]);
  entropy -= correct_aa_freqs[6] * log(correct_aa_freqs[6]);
  entropy -= correct_aa_freqs[11] * log(correct_aa_freqs[11]);
  entropy -= correct_aa_freqs[12] * log(correct_aa_freqs[12]);
  entropy -= correct_aa_freqs[13] * log(correct_aa_freqs[13]);
  entropy -= correct_aa_freqs[15] * log(correct_aa_freqs[15]);
  entropy -= correct_aa_freqs[16] * log(correct_aa_freqs[16]);
  entropy -= correct_aa_freqs[17] * log(correct_aa_freqs[17]);
  entropy -= correct_aa_freqs[19] * log(correct_aa_freqs[19]);
  BOOST_CHECK_CLOSE((*prof)[2].GetEntropy(), entropy, Real(1e-5));

  // look at avg entropy...
  BOOST_CHECK_CLOSE(prof->GetAverageEntropy(), 1.48704, Real(1e-2));

  // reload with LoadSequenceProfile
  ProfileHandlePtr prof2 = LoadSequenceProfile("testfiles/test_hmm.hhm");
  BOOST_CHECK(*prof == *prof2);

  // reload
  io_handler->Import(*prof2, "testfiles/test_hmm.hhm");
  BOOST_CHECK(*prof == *prof2);
}

BOOST_AUTO_TEST_CASE(pssm_loading)
{
  // check pssm factory
  ProfileIOHandlerFactoryBasePtr fac(new PssmIOHandlerFactory);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test.hhm"),     false);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test.hhm.gz"),  false);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test.pssm"),    true);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test.pssm.gz"), true);
  BOOST_CHECK_EQUAL(fac->ProvidesImport("test", "pssm"), true);
  BOOST_CHECK_EQUAL(fac->ProvidesExport("test.pssm"),    false);
  BOOST_CHECK_EQUAL(fac->ProvidesExport("test.pssm.gz"), false);
  BOOST_CHECK_EQUAL(fac->ProvidesExport("test", "pssm"), false);
  ProfileIOHandlerPtr io_handler = fac->Create();

  // load file
  ProfileHandlePtr prof(new ProfileHandle);
  io_handler->Import(*prof, "testfiles/test_pssm.pssm");

  // compare manually
  ProfileHandle prof_pssm_ref;
  prof_pssm_ref.SetNullModel(ProfileColumn::BLOSUMNullModel());
  ProfileColumn pc1, pc2, pc3, pc4, pc5;
  pc1.SetFreq('R', 0.74);
  pc1.SetFreq('Q', 0.01);
  pc1.SetFreq('K', 0.25);

  pc2.SetFreq('A', 0.02/0.99);
  pc2.SetFreq('R', 0.61/0.99);
  pc2.SetFreq('N', 0.03/0.99);
  pc2.SetFreq('Q', 0.02/0.99);
  pc2.SetFreq('E', 0.02/0.99);
  pc2.SetFreq('G', 0.02/0.99);
  pc2.SetFreq('H', 0.01/0.99);
  pc2.SetFreq('I', 0.01/0.99);
  pc2.SetFreq('L', 0.01/0.99);
  pc2.SetFreq('K', 0.13/0.99);
  pc2.SetFreq('P', 0.03/0.99);
  pc2.SetFreq('S', 0.06/0.99);
  pc2.SetFreq('T', 0.01/0.99);
  pc2.SetFreq('V', 0.01/0.99);
  pc3.SetFreq('A', 0.17/0.99);
  pc3.SetFreq('R', 0.07/0.99);
  pc3.SetFreq('N', 0.03/0.99);

  pc3.SetFreq('A', 0.17);
  pc3.SetFreq('R', 0.07);
  pc3.SetFreq('N', 0.03);
  pc3.SetFreq('D', 0.02);
  pc3.SetFreq('C', 0.01);
  pc3.SetFreq('Q', 0.03);
  pc3.SetFreq('E', 0.02);
  pc3.SetFreq('G', 0.03);
  pc3.SetFreq('L', 0.02);
  pc3.SetFreq('K', 0.04);
  pc3.SetFreq('M', 0.01);
  pc3.SetFreq('F', 0.01);
  pc3.SetFreq('P', 0.06);
  pc3.SetFreq('S', 0.37);
  pc3.SetFreq('T', 0.08);
  pc3.SetFreq('V', 0.03);

  pc4.SetFreq('A', 0.11/0.98);
  pc4.SetFreq('R', 0.06/0.98);
  pc4.SetFreq('D', 0.01/0.98);
  pc4.SetFreq('Q', 0.02/0.98);
  pc4.SetFreq('E', 0.02/0.98);
  pc4.SetFreq('I', 0.01/0.98);
  pc4.SetFreq('L', 0.03/0.98);
  pc4.SetFreq('K', 0.04/0.98);
  pc4.SetFreq('F', 0.01/0.98);
  pc4.SetFreq('P', 0.58/0.98);
  pc4.SetFreq('S', 0.03/0.98);
  pc4.SetFreq('T', 0.05/0.98);
  pc4.SetFreq('V', 0.01/0.98);

  pc5.SetFreq('A', 0.04);
  pc5.SetFreq('R', 0.01);
  pc5.SetFreq('N', 0.02);
  pc5.SetFreq('D', 0.03);
  pc5.SetFreq('Q', 0.02);
  pc5.SetFreq('E', 0.03);
  pc5.SetFreq('H', 0.01);
  pc5.SetFreq('I', 0.01);
  pc5.SetFreq('L', 0.04);
  pc5.SetFreq('K', 0.05);
  pc5.SetFreq('M', 0.01);
  pc5.SetFreq('F', 0.01);
  pc5.SetFreq('P', 0.63);
  pc5.SetFreq('S', 0.03);
  pc5.SetFreq('T', 0.02);
  pc5.SetFreq('V', 0.04);

  prof_pssm_ref.AddColumn(pc1);
  prof_pssm_ref.AddColumn(pc2);
  prof_pssm_ref.AddColumn(pc3);
  prof_pssm_ref.AddColumn(pc4);
  prof_pssm_ref.AddColumn(pc5);
  prof_pssm_ref.SetSequence("RRSPP");

  compareProfiles(*prof, prof_pssm_ref, 1e-3);

  // reload with LoadSequenceProfile
  ProfileHandlePtr prof2 = LoadSequenceProfile("testfiles/test_pssm.pssm");
  BOOST_CHECK(*prof == *prof2);

  // reload
  io_handler->Import(*prof2, "testfiles/test_pssm.pssm");
  BOOST_CHECK(*prof == *prof2);
}

BOOST_AUTO_TEST_CASE(hhm_loading_file_vs_str)
{
  // load from file
  ProfileHandlePtr prof_from_file(new ProfileHandle);
  ProfileIOHandlerPtr io_handler(new HhmIOHandler);
  io_handler->Import(*prof_from_file, "testfiles/test_hmm.hhm");

  // load from string
  std::ifstream in("testfiles/test_hmm.hhm", std::ios::in | std::ios::binary);
  std::string str;
  if (in) {
    in.seekg(0, std::ios::end);
    str.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&str[0], str.size());
    in.close();
  } else {
    throw IOException("Cannot read hhm file in unit test.");
  }
  ProfileHandlePtr prof_from_string(new ProfileHandle);
  io_handler->ImportFromString(*prof_from_string, str);

  compareProfiles(*prof_from_file, *prof_from_string, 1e-3);
}

BOOST_AUTO_TEST_CASE(pssm_loading_file_vs_str)
{
  // load from file
  ProfileHandlePtr prof_from_file(new ProfileHandle);
  ProfileIOHandlerPtr io_handler(new PssmIOHandler);
  io_handler->Import(*prof_from_file, "testfiles/test_pssm.pssm");

  // load from string
  std::ifstream in("testfiles/test_pssm.pssm", std::ios::in | std::ios::binary);
  std::string str;
  if (in) {
    in.seekg(0, std::ios::end);
    str.resize(in.tellg());
    in.seekg(0, std::ios::beg);
    in.read(&str[0], str.size());
    in.close();
  } else {
    throw IOException("Cannot read pssm file in unit test.");
  }
  ProfileHandlePtr prof_from_string(new ProfileHandle);
  io_handler->ImportFromString(*prof_from_string, str);

  compareProfiles(*prof_from_file, *prof_from_string, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END();
