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
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <ost/seq/alg/merge_pairwise_alignments.hh>
#include <ost/seq/alg/sequence_identity.hh>
#include <ost/seq/alg/sequence_similarity.hh>
#include <ost/seq/alg/ins_del.hh>
#include <ost/seq/alg/conservation.hh>
#include <ost/seq/alg/subst_weight_matrix.hh>
#include <ost/seq/alg/local_align.hh>
#include <ost/seq/alg/global_align.hh>
#include <ost/seq/alg/semiglobal_align.hh>
#include <ost/seq/alg/entropy.hh>
#include <ost/seq/alg/pair_subst_weight_matrix.hh>
#include <ost/seq/alg/contact_weight_matrix.hh>
#include <ost/seq/alg/contact_prediction_score.hh>
#include <ost/seq/alg/clip_alignment.hh>
#include <ost/seq/alg/distance_map.hh>
#include <ost/seq/alg/variance_map.hh>
#include <ost/seq/alg/hmm_pseudo_counts.hh>
#include <ost/seq/alg/hmm_score.hh>

using namespace boost::python;
using namespace ost::seq;
using namespace ost::seq::alg;

////////////////////////////////////////////////////////////////////
// wrappers
namespace {
tuple WrapDistancesGetMin(const Distances& distance) {
  const std::pair<Real, int> dist = distance.GetMin();
  return boost::python::make_tuple(dist.first, dist.second);
}
tuple WrapDistancesGetMax(const Distances& distance) {
  const std::pair<Real, int> dist = distance.GetMax();
  return boost::python::make_tuple(dist.first, dist.second);
}
tuple WrapDistancesGetDataElement(const Distances& distance, uint index) {
  const std::pair<Real, int> dist = distance.GetDataElement(index);
  return boost::python::make_tuple(dist.first, dist.second);
}

template <typename T>
list GetList(const T& data, uint num_rows, uint num_cols) {
  list ret;
  for (uint row = 0; row < num_rows; ++row) {
    list my_row;
    for (uint col = 0; col < num_cols; ++col) {
      my_row.append(data(row, col));
    }
    ret.append(my_row);
  }
  return ret;
}

list VarMapGetData(const VarianceMapPtr v_map) {
  return GetList(*v_map, v_map->GetSize(), v_map->GetSize());
}

list DistToMeanGetData(const Dist2MeanPtr d2m) {
  return GetList(*d2m, d2m->GetNumResidues(), d2m->GetNumStructures());
}

void AAPseudoCountsSimple(ProfileHandle& profile, Real a, Real b, Real c) {
  AddAAPseudoCounts(profile, a, b, c);
}

void AAPseudoCountsAngermueller(ProfileHandle& profile, const ContextProfileDB& db,
                                Real a, Real b, Real c) {
  AddAAPseudoCounts(profile, db, a, b, c);
}

} // anon ns
////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////
// Work on alignments
void export_aln_alg()
{
  enum_<RefMode::Type>("RefMode")
    .value("ALIGNMENT", RefMode::ALIGNMENT)
    .value("LONGER_SEQUENCE", RefMode::LONGER_SEQUENCE)
    .export_values()
  ;

  def("SequenceIdentity", &SequenceIdentity, 
      (arg("ref_mode")=RefMode::ALIGNMENT, arg("seq_a")=0, arg("seq_b")=1));

  def("SequenceSimilarity", &SequenceSimilarity,
      (arg("aln"),arg("subst_weight"),arg("normalize")=false,
       arg("seq_a")=0,arg("seq_b")=1));

  class_<AlignedRegionList>("AlignedRegionList", init<>())
    .def(vector_indexing_suite<AlignedRegionList>())
  ;

  class_<InsDel>("InsDel", init<const AlignmentHandle&, int, int>())
    .def(init<const AlignmentHandle&>())
    .def("GetDeletions", &InsDel::GetDeletions)
    .def("GetInsertions", &InsDel::GetInsertions)
  ;
  
  def("MergePairwiseAlignments", &MergePairwiseAlignments);
  def("Conservation", &Conservation, (arg("assign")=true, arg("prop_name")="cons", arg("ignore_gap")=false));
  def("LocalAlign", &LocalAlign, (arg("seq1"), arg("seq2"),arg("subst_weight"), 
      arg("gap_open")=-5, arg("gap_ext")=-2));
  def("GlobalAlign", &GlobalAlign,(arg("seq1"),arg("seq2"),arg("subst_weight"), 
      arg("gap_open")=-5, arg("gap_ext")=-2));
  def("SemiGlobalAlign", &SemiGlobalAlign,(arg("seq1"),arg("seq2"),arg("subst_weight"), 
      arg("gap_open")=-5, arg("gap_ext")=-2));
  def("ShannonEntropy", &ShannonEntropy, (arg("aln"), arg("ignore_gaps")=true));
}

////////////////////////////////////////////////////////////////////
// Contact Prediction
void export_contact_prediction()
{
  class_<PairSubstWeightMatrix>("PairSubstWeightMatrix",init< std::vector <std::vector <std::vector <std::vector <Real> > > >,std::vector <char> >())
    .def(init<>())
    .def_readonly("weights",&PairSubstWeightMatrix::weights)
    .def_readonly("aa_list",&PairSubstWeightMatrix::aa_list)
  ;

  class_<ContactWeightMatrix>("ContactWeightMatrix",init< std::vector <std::vector <Real> >,std::vector <char> >())
    .def(init<>())
    .def_readonly("weights",&ContactWeightMatrix::weights)
    .def_readonly("aa_list",&ContactWeightMatrix::aa_list)
  ;

  class_<ContactPredictionScoreResult>("ContactPredictionScoreResult",no_init)
    .def_readonly("matrix",&ContactPredictionScoreResult::matrix)
    .def_readonly("sorted_indices",&ContactPredictionScoreResult::sorted_indices)
    .def("RefreshSortedIndices",&ContactPredictionScoreResult::RefreshSortedIndices)
    .def("GetScore",&ContactPredictionScoreResult::GetScore,(arg("i"),arg("j")))
    .def("SetScore",&ContactPredictionScoreResult::SetScore,(arg("i"),arg("j"),arg("score")))
  ;
  def("CalculateMutualInformation", &CalculateMutualInformation,(arg("aln"),arg("w")=LoadConstantContactWeightMatrix(),
                                                        arg("apc_correction")=true,arg("zpx_transformation")=true,arg("small_number_correction")=0.05));
  def("CalculateContactScore", &CalculateContactScore,(arg("aln"), arg("w")=LoadDefaultContactWeightMatrix()));
  def("CalculateContactSubstitutionScore", &CalculateContactSubstitutionScore,(arg("aln"), arg("ref_seq_index")=0, arg("w")=LoadDefaultPairSubstWeightMatrix()));
  def("LoadDefaultContactWeightMatrix",LoadDefaultContactWeightMatrix);
  def("LoadConstantContactWeightMatrix",LoadConstantContactWeightMatrix);
  def("LoadDefaultPairSubstWeightMatrix",LoadDefaultPairSubstWeightMatrix);

  // NOTE: anything after this is within SubstWeightMatrix scope
  scope mat_scope = class_<SubstWeightMatrix, SubstWeightMatrixPtr>
                      ("SubstWeightMatrix", init<>())
                      .def("GetWeight", &SubstWeightMatrix::GetWeight)
                      .def("SetWeight", &SubstWeightMatrix::SetWeight)
                      .def("GetMinWeight", &SubstWeightMatrix::GetMinWeight)
                      .def("GetMaxWeight", &SubstWeightMatrix::GetMaxWeight)
                      .def("AssignPreset", &SubstWeightMatrix::AssignPreset)
  ;

  enum_<SubstWeightMatrix::Preset>("Preset")
    .value("BLOSUM45", SubstWeightMatrix::BLOSUM45)
    .value("BLOSUM62", SubstWeightMatrix::BLOSUM62)
    .value("BLOSUM80", SubstWeightMatrix::BLOSUM80)
    .value("BLOSUM100", SubstWeightMatrix::BLOSUM100)
  ;
}

////////////////////////////////////////////////////////////////////
// getting/analyzing distance matrices from alignments
void export_distance_analysis()
{
  def("ClipAlignment", &ClipAlignment, (arg("aln"), arg("n_seq_thresh")=2, 
      arg("set_offset")=true, arg("remove_empty")=true));
  def("CreateDistanceMap", &CreateDistanceMap, (arg("aln")));
  def("CreateVarianceMap", &CreateVarianceMap, (arg("d_map"), arg("sigma")=25));
  def("CreateDist2Mean", &CreateDist2Mean, (arg("d_map")));

  class_<Distances>("Distances", no_init)
    .def("GetDataSize", &Distances::GetDataSize)
    .def("GetAverage", &Distances::GetAverage)
    .def("GetMin", &WrapDistancesGetMin)
    .def("GetMax", &WrapDistancesGetMax)
    .def("GetDataElement", &WrapDistancesGetDataElement, (arg("index")))
    .def("GetStdDev", &Distances::GetStdDev)
    .def("GetWeightedStdDev", &Distances::GetWeightedStdDev, (arg("sigma")))
    .def("GetNormStdDev", &Distances::GetNormStdDev)
  ;

  class_<DistanceMap, DistanceMapPtr,
         boost::noncopyable>("DistanceMap", no_init)
    .def("GetDistances", &DistanceMap::GetDistances,
         return_value_policy<reference_existing_object>(),
         (arg("i_res1"), arg("i_res2")))
    .def("GetSize", &DistanceMap::GetSize)
    .def("GetNumStructures", &DistanceMap::GetNumStructures)
  ;

  class_<VarianceMap, VarianceMapPtr,
         boost::noncopyable>("VarianceMap", no_init)
    .def("Get", &VarianceMap::Get, (arg("i_res1"), arg("i_res2")),
         return_value_policy<return_by_value>())
    .def("GetSize", &VarianceMap::GetSize)
    .def("Min", &VarianceMap::Min)
    .def("Max", &VarianceMap::Max)
    .def("ExportDat", &VarianceMap::ExportDat, (arg("file_name")))
    .def("ExportCsv", &VarianceMap::ExportCsv, (arg("file_name")))
    .def("ExportJson", &VarianceMap::ExportJson, (arg("file_name")))
    .def("GetJsonString", &VarianceMap::GetJsonString)
    .def("GetData", &VarMapGetData)
  ;

  class_<Dist2Mean, Dist2MeanPtr,
         boost::noncopyable>("Dist2Mean", no_init)
    .def("Get", &Dist2Mean::Get, (arg("i_res"), arg("i_str")))
    .def("GetNumResidues", &Dist2Mean::GetNumResidues)
    .def("GetNumStructures", &Dist2Mean::GetNumStructures)
    .def("ExportDat", &Dist2Mean::ExportDat, (arg("file_name")))
    .def("ExportCsv", &Dist2Mean::ExportCsv, (arg("file_name")))
    .def("ExportJson", &Dist2Mean::ExportJson, (arg("file_name")))
    .def("GetJsonString", &Dist2Mean::GetJsonString)
    .def("GetData", &DistToMeanGetData)
  ;
}

////////////////////////////////////////////////////////////////////
// algorithms involving hmms
void export_hmm_algorithms() {

  class_<ContextProfile>("ContextProfile", init<int>())
    .def("SetWeight",&ContextProfile::SetWeight, (arg("pos"), arg("olc"), arg("weight")))
    .def("SetPseudoCount",&ContextProfile::SetPseudoCount, (arg("olc"), arg("count")))
    .def("SetBias",&ContextProfile::SetBias, (arg("bias")))
    .def("GetWeight", &ContextProfile::GetWeight, (arg("pos"), arg("olc")))
    .def("GetPseudoCount", &ContextProfile::GetPseudoCount,(arg("olc")))
    .def("GetBias", &ContextProfile::GetBias)
    .def("GetLength", &ContextProfile::GetLength)
  ;

  class_<ContextProfileDB, ContextProfileDBPtr>("ContextProfileDB", init<>())
    .def("__len__",&ContextProfileDB::size)
    .def("__getitem__",&ContextProfileDB::at,return_value_policy<reference_existing_object>(), (arg("idx")))
    .def("Save", &ContextProfileDB::Save, (arg("filename")))
    .def("Load", &ContextProfileDB::Load, (arg("filename"))).staticmethod("Load")
    .def("FromCRF", &ContextProfileDB::FromCRF, (arg("filename"))).staticmethod("FromCRF")
    .def("AddProfile", &ContextProfileDB::AddProfile, (arg("profile")))
  ;

  def("AddAAPseudoCounts", &AAPseudoCountsSimple, (arg("profile"), 
                                                   arg("a")=1.0,
                                                   arg("b")=1.5,
                                                   arg("c")=1.0));
  def("AddAAPseudoCounts", &AAPseudoCountsAngermueller, (arg("profile"), 
                                                         arg("context_profile_db"),
                                                         arg("a")=0.9,
                                                         arg("b")=4.0,
                                                         arg("c")=1.0));
  def("AddTransitionPseudoCounts", &AddTransitionPseudoCounts, (arg("profile"),
                                                                arg("gapb")=1.0,
                                                                arg("gapd")=0.15,
                                                                arg("gape")=1.0));
  def("AddNullPseudoCounts", &AddNullPseudoCounts, (arg("profile")));
  def("HMMScore", &HMMScore, (arg("profile_0"), arg("profile_1"), arg("alignment"),
                              arg("s_0_idx"), arg("s_1_idx"), 
                              arg("match_score_offset")=-0.03,
                              arg("correl_score_weight")=0.1,
                              arg("del_start_penalty_factor")=0.6,
                              arg("del_extend_penalty_factor")=0.6,
                              arg("ins_start_penalty_factor")=0.6,
                              arg("ins_extend_penalty_factor")=0.6));
}

BOOST_PYTHON_MODULE(_ost_seq_alg)
{
  export_aln_alg();
  export_contact_prediction();
  export_distance_analysis();
  export_hmm_algorithms();
}
