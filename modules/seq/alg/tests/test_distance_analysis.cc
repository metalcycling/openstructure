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
  Author: Gerardo Tauriello
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <ost/mol/atom_view.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/io/mol/load_entity.hh>
#include <ost/seq/alg/clip_alignment.hh>
#include <ost/seq/alg/distance_map.hh>
#include <ost/seq/alg/variance_map.hh>
#include <ost/integrity_error.hh>
#include <ost/mol/alg/local_dist_diff_test.hh>

using namespace ost;
using namespace ost::seq;

BOOST_AUTO_TEST_SUITE(ost_seq_alg);

// HELPERS
AlignmentHandle GetClippingAln() {
  AlignmentHandle aln = CreateAlignment();
  aln.AddSequence(CreateSequence("ref", "DDFAGTHN"));
  aln.AddSequence(CreateSequence("A",   "-DFAG---"));
  aln.AddSequence(CreateSequence("B",   "-----THN"));
  aln.AddSequence(CreateSequence("C",   "--FAG---"));
  aln.AddSequence(CreateSequence("D",   "DDFAG---"));
  // attach view
  mol::EntityHandle ent = io::LoadEntity("testfiles/dummy.pdb");
  for (int i_s = 1; i_s < aln.GetCount(); ++i_s) {
    aln.AttachView(i_s, ent.CreateFullView());
    aln.SetSequenceOffset(i_s, aln.GetSequence(i_s).GetFirstNonGap());
  }
  return aln;
}

AlignmentHandle GetDistanceMapAln() {
  // prepare alignment
  AlignmentHandle aln = CreateAlignment();
  aln.AddSequence(CreateSequence("ref", "DDFAGTHN"));
  aln.AddSequence(CreateSequence("A",   "DDFAGT--"));
  aln.AddSequence(CreateSequence("B",   "--FAGTH-"));
  // attach views
  mol::EntityHandle ent = io::LoadEntity("testfiles/dummy.pdb");
  // move CA of res. at index 3 (res. num = 4) for ent
  mol::AtomHandle atom = ent.FindAtom("A", 4, "CA");
  ent.EditXCS(mol::UNBUFFERED_EDIT).SetAtomPos(atom, geom::Vec3(0,0,0));
  BOOST_CHECK_EQUAL(ent.FindAtom("A", 4, "CA").GetPos(), geom::Vec3(0,0,0));
  mol::EntityHandle ent2 = io::LoadEntity("testfiles/dummy.pdb");
  aln.AttachView(1, ent.CreateFullView());
  aln.AttachView(2, ent2.CreateFullView());
  aln.SetSequenceOffset(2, 2);
  return aln;
}

alg::DistanceMapPtr GetDistanceMap() {
  AlignmentHandle aln = GetDistanceMapAln();
  return alg::CreateDistanceMap(aln);
}

// TESTS
BOOST_AUTO_TEST_CASE(test_clip_alignment_no_change) {
  // clip it at level 2 (no change!)
  AlignmentHandle aln_cut = GetClippingAln();
  int clip_start = alg::ClipAlignment(aln_cut, 2);
  BOOST_CHECK_EQUAL(clip_start, 0);
  BOOST_CHECK_EQUAL(aln_cut.GetLength(), 8);
  BOOST_CHECK_EQUAL(aln_cut.GetCount(), 5);
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(0).GetString(), "DDFAGTHN");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(1).GetString(), "-DFAG---");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(2).GetString(), "-----THN");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(3).GetString(), "--FAG---");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(4).GetString(), "DDFAG---");
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(0), 0);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(1), 1);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(2), 5);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(3), 2);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(4), 0);
}

BOOST_AUTO_TEST_CASE(test_clip_alignment_no_offset) {
  // clip it at level 3 without updating offsets and deleting
  AlignmentHandle aln_cut = GetClippingAln();
  int clip_start = alg::ClipAlignment(aln_cut, 3, false, false);
  BOOST_CHECK_EQUAL(clip_start, 1);
  BOOST_CHECK_EQUAL(aln_cut.GetLength(), 4);
  BOOST_CHECK_EQUAL(aln_cut.GetCount(), 5);
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(0).GetString(), "DFAG");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(1).GetString(), "DFAG");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(2).GetString(), "----");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(3).GetString(), "-FAG");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(4).GetString(), "DFAG");
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(0), 0);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(1), 1);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(2), 5);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(3), 2);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(4), 0);
}

BOOST_AUTO_TEST_CASE(test_clip_alignment_def) {
  // clip it at level 3 with updating offsets and deleting
  AlignmentHandle aln_cut = GetClippingAln();
  int clip_start = alg::ClipAlignment(aln_cut, 3, true, true);
  BOOST_CHECK_EQUAL(clip_start, 1);
  BOOST_CHECK_EQUAL(aln_cut.GetLength(), 4);
  BOOST_CHECK_EQUAL(aln_cut.GetCount(), 4);
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(0).GetString(), "DFAG");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(1).GetString(), "DFAG");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(2).GetString(), "-FAG");
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(3).GetString(), "DFAG");
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(0), 0); // REF!
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(1), 1);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(2), 2);
  BOOST_CHECK_EQUAL(aln_cut.GetSequenceOffset(3), 1);
}

BOOST_AUTO_TEST_CASE(test_clip_alignment_remove_all) {
  // clip it at level 5 to remove all
  AlignmentHandle aln_cut = GetClippingAln();
  int clip_start = alg::ClipAlignment(aln_cut, 5, true, true);
  BOOST_CHECK_EQUAL(clip_start, -1);
  BOOST_CHECK_EQUAL(aln_cut.GetCount(), 0);
}

BOOST_AUTO_TEST_CASE(test_clip_alignment_generic) {
  // generic tests
  AlignmentHandle aln_cut = CreateAlignment();
  int clip_start = alg::ClipAlignment(aln_cut);
  BOOST_CHECK_EQUAL(clip_start, -1);
  BOOST_CHECK_EQUAL(aln_cut.GetCount(), 0);
  aln_cut.AddSequence(CreateSequence("ref", "DDFAG"));
  aln_cut.AddSequence(CreateSequence("A",   "-----"));
  clip_start = alg::ClipAlignment(aln_cut, 1, true, true);
  BOOST_CHECK_EQUAL(clip_start, 0);
  BOOST_CHECK_EQUAL(aln_cut.GetCount(), 1);
  BOOST_CHECK_EQUAL(aln_cut.GetSequence(0).GetString(), "DDFAG");
  aln_cut.AddSequence(CreateSequence("A",   "-----"));
  clip_start = alg::ClipAlignment(aln_cut, 2, true, true);
  BOOST_CHECK_EQUAL(clip_start, -1);
  BOOST_CHECK_EQUAL(aln_cut.GetCount(), 0);
}

BOOST_AUTO_TEST_CASE(test_distance_map) {
  // setup
  AlignmentHandle aln = GetDistanceMapAln();
  alg::DistanceMapPtr d_map = alg::CreateDistanceMap(aln);
  geom::Vec3 pos2 = aln.GetResidue(2, 2).FindAtom("CA").GetPos();
  geom::Vec3 pos3 = aln.GetResidue(2, 3).FindAtom("CA").GetPos();
  geom::Vec3 pos4 = aln.GetResidue(2, 4).FindAtom("CA").GetPos();
  // check map size
  BOOST_CHECK_EQUAL(d_map->GetSize(), 8u);
  BOOST_CHECK_EQUAL(d_map->GetNumStructures(), 2u);
  // check single entry
  const alg::Distances& di = d_map->GetDistances(2,4);
  BOOST_CHECK(di == d_map->GetDistances(4,2));
  BOOST_CHECK_EQUAL(di.GetDataSize(), 2u);
  std::pair<Real,int> data0 = di.GetDataElement(0);
  std::pair<Real,int> data1 = di.GetDataElement(1);
  BOOST_CHECK_EQUAL(data0.second + data1.second, 1+2);
  BOOST_CHECK_EQUAL(data0.first, geom::Length(pos4 - pos2));
  BOOST_CHECK_EQUAL(data1.first, data0.first);
  // check modified entry
  const alg::Distances& di2 = d_map->GetDistances(3,4);
  BOOST_CHECK_EQUAL(di2.GetDataSize(), 2u);
  data0 = di2.GetDataElement(0);
  data1 = di2.GetDataElement(1);
  if (data0.second == 2) std::swap(data0, data1);
  BOOST_CHECK_EQUAL(data0.second, 1);
  BOOST_CHECK_EQUAL(data1.second, 2);
  BOOST_CHECK_EQUAL(data0.first, geom::Length(pos4));
  BOOST_CHECK_EQUAL(data1.first, geom::Length(pos4 - pos3));
  BOOST_CHECK(data0.first != data1.first);
  // min/max and other features
  std::pair<Real,int> min_data = di2.GetMin();
  std::pair<Real,int> max_data = di2.GetMax();
  if (data0.first > data1.first) {
    BOOST_CHECK(max_data == data0);
    BOOST_CHECK(min_data == data1);
  } else {
    BOOST_CHECK(max_data == data1);
    BOOST_CHECK(min_data == data0);
  }
  const Real exp_avg = (data0.first + data1.first)/2;
  BOOST_CHECK_CLOSE(di2.GetAverage(), exp_avg, 1e-3);
  const Real exp_std = std::abs(data0.first - exp_avg);
  BOOST_CHECK_CLOSE(di2.GetStdDev(), exp_std, 1e-3);
  BOOST_CHECK_CLOSE(di2.GetWeightedStdDev(2), exp_std * std::exp(-exp_avg/4),
                    1e-3);
  BOOST_CHECK_CLOSE(di2.GetNormStdDev(), exp_std/exp_avg, 1e-3);
  // check other entries
  BOOST_CHECK_EQUAL(d_map->GetDistances(1,2).GetDataSize(), 1u);
  BOOST_CHECK_EQUAL(d_map->GetDistances(1,7).GetDataSize(), 0u);
  BOOST_CHECK_EQUAL(d_map->GetDistances(4,7).GetDataSize(), 0u);
}

BOOST_AUTO_TEST_CASE(test_variance_map) {
  alg::DistanceMapPtr d_map = GetDistanceMap();
  // check variance map
  alg::VarianceMapPtr std_mat = alg::CreateVarianceMap(d_map);
  BOOST_CHECK_EQUAL(std_mat->GetSize(), 8);
  for (uint row = 0; row < 8; ++row) {
    for (uint col = 0; col < row; ++col) {
      if ((col == 3 && row <= 5) || (row == 3 && col >= 2)) {
        BOOST_CHECK(std_mat->Get(row, col) > 0);
      } else {
        BOOST_CHECK_EQUAL(std_mat->Get(row, col), Real(0));
      }
      BOOST_CHECK_EQUAL(std_mat->Get(row, col), std_mat->Get(col, row));
    }
    BOOST_CHECK_EQUAL(std_mat->Get(row, row), Real(0));
  }
  Real max_std = std::max(std_mat->Get(3,4), std_mat->Get(2,3));
  max_std = std::max(std_mat->Get(3,5), max_std);
  BOOST_CHECK_EQUAL(std_mat->Min(), Real(0));
  BOOST_CHECK_EQUAL(std_mat->Max(), max_std);
  // check json export
  String json = std_mat->GetJsonString();
  BOOST_CHECK_EQUAL(std::count(json.begin(), json.end(), '['), 8+1);
  BOOST_CHECK_EQUAL(std::count(json.begin(), json.end(), ']'), 8+1);
  BOOST_CHECK_EQUAL(std::count(json.begin(), json.end(), ','), 7*8 + 7);
}

BOOST_AUTO_TEST_CASE(test_dist_to_mean) {
  alg::DistanceMapPtr d_map = GetDistanceMap();
  // check distance to mean
  alg::Dist2MeanPtr dist_to_mean = alg::CreateDist2Mean(d_map);
  uint num_res = dist_to_mean->GetNumResidues();
  uint num_str = dist_to_mean->GetNumStructures();
  BOOST_CHECK_EQUAL(num_res, 8u);
  BOOST_CHECK_EQUAL(num_str, 2u);
  for (uint row = 0; row < num_res; ++row) {
    for (uint col = 0; col < num_str; ++col) {
      if (row >= 2 && row <= 5) {
        BOOST_CHECK(dist_to_mean->Get(row, col) != 0);
      } else {
        if (dist_to_mean->Get(row, col) != 0) std::cout << row << ", " << col;
        BOOST_CHECK_EQUAL(dist_to_mean->Get(row, col), Real(0));
      }
    }
  }
  // check json export
  String json = dist_to_mean->GetJsonString();
  BOOST_CHECK_EQUAL(std::count(json.begin(), json.end(), '['), 8+1);
  BOOST_CHECK_EQUAL(std::count(json.begin(), json.end(), ']'), 8+1);
  BOOST_CHECK_EQUAL(std::count(json.begin(), json.end(), ','), 1*8 + 7);
}

BOOST_AUTO_TEST_CASE(test_avg_lddt) {

  // deliberately load different structures than in other tests
  // as we need a bit a bigger structure to also test things like the 15A 
  // threshold in lDDT. We compare the ouput of MeanlDDTHA with the classic
  // lDDT. As there are only two structures in the alignment, there
  // is no real average and the outputs from the two algorithms should be 
  // consistent.

  mol::EntityView calmodulin_1 = 
  io::LoadEntity("testfiles/calmodulin_01.pdb").Select("aname=CA");
  mol::EntityView calmodulin_2 = 
  io::LoadEntity("testfiles/calmodulin_02.pdb").Select("aname=CA");

  std::stringstream ss;
  ost::mol::ResidueViewList r_list = calmodulin_1.GetResidueList();
  for(auto r_it = r_list.begin(); r_it != r_list.end(); ++r_it) {
    ss << r_it->GetOneLetterCode();
  }

  AlignmentHandle aln = CreateAlignment();
  aln.AddSequence(CreateSequence("ref", ss.str()));
  aln.AddSequence(CreateSequence("A",   ss.str()));
  aln.AddSequence(CreateSequence("B",   ss.str()));

  aln.AttachView(1, calmodulin_1);
  aln.AttachView(2, calmodulin_2);

  ost::mol::alg::GlobalRDMap dl_1 = mol::alg::CreateDistanceList(calmodulin_1, 15.0);
  ost::mol::alg::GlobalRDMap dl_2 = mol::alg::CreateDistanceList(calmodulin_2, 15.0);

  ost::mol::alg::LDDTHA(calmodulin_1, dl_2);
  ost::mol::alg::LDDTHA(calmodulin_2, dl_1);

  alg::DistanceMapPtr dmap = alg::CreateDistanceMap(aln);
  alg::MeanlDDTPtr mean_lddt = alg::CreateMeanlDDTHA(dmap);

  r_list = calmodulin_1.GetResidueList();
  int r_idx = 0;
  for(auto r_it = r_list.begin(); r_it != r_list.end(); ++r_it, ++r_idx) {
    BOOST_CHECK_CLOSE(r_it->GetFloatProp("locallddt"), mean_lddt->Get(r_idx, 0), 1e-4);
  }

  r_list = calmodulin_2.GetResidueList();
  r_idx = 0;
  for(auto r_it = r_list.begin(); r_it != r_list.end(); ++r_it, ++r_idx) {
    BOOST_CHECK_CLOSE(r_it->GetFloatProp("locallddt"), mean_lddt->Get(r_idx, 1), 1e-4);
  }
}

BOOST_AUTO_TEST_SUITE_END();
