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
#ifndef OST_SEQ_ALG_CONTACT_PREDICTION_SCORE_HH
#define OST_SEQ_ALG_CONTACT_PREDICTION_SCORE_HH
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/module_config.hh>
#include <ost/seq/alg/pair_subst_weight_matrix.hh>
#include <ost/seq/alg/contact_weight_matrix.hh>

/*
  Authors: Niklaus Johner
  Functions and objects used to calculate different scores for predicting contacts
  from a multiple sequence alignment (MSA). Several scores are available:
  1. Mutual information. The average product correction (apc_correction) 
  adds a correction term to account for common phylogenetic origin and inherent
  entropy of each column. There is also a correction for small number of sequences 
  in the alignment and the possibility to use weighted mutual information.
  2. ContactScore. For each residue pair, this score calculates the average interaction
  energy for the amino-acid pairs in the corresponding columns of the MSA
  3. ContactSubstitutionScore. For each residue pair, this score calculates 
  the average energy of the substitutions of amino-acid pairs in the corresponding
  columns of the MSA.
 */
namespace ost { namespace seq { namespace alg {

struct DLLEXPORT_OST_SEQ_ALG ContactPredictionScoreResult {
  std::vector <std::vector <Real> > matrix;
  std::vector <std::pair <int,int> > sorted_indices;
  void SetScore(int i, int j, Real score);
  Real GetScore(int i, int j);
  void RefreshSortedIndices();
  ContactPredictionScoreResult(std::vector <std::vector <Real> >);
};

ContactPredictionScoreResult DLLEXPORT_OST_SEQ_ALG CalculateContactScore(const AlignmentHandle& aln, ContactWeightMatrix w=LoadDefaultContactWeightMatrix());
ContactPredictionScoreResult DLLEXPORT_OST_SEQ_ALG CalculateContactSubstitutionScore(const AlignmentHandle& aln, int ref_seq_index=0, PairSubstWeightMatrix w=LoadDefaultPairSubstWeightMatrix());
ContactPredictionScoreResult DLLEXPORT_OST_SEQ_ALG CalculateMutualInformation(const AlignmentHandle& aln, ContactWeightMatrix w=LoadConstantContactWeightMatrix(),bool apc_correction=true, 
                            bool zpx_transformation=true,float small_number_correction=0.05);

}}}
#endif
