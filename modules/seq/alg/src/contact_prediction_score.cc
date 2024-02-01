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

#include <ost/seq/aligned_column.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/alg/contact_prediction_score.hh>
#include <math.h>

namespace {
  typedef std::pair<float, std::pair<int,int> > val_indices_pair;

  bool pair_comparator ( const val_indices_pair& l, const val_indices_pair& r)
    { return l.first > r.first; }
}

namespace ost { namespace seq { namespace alg {

void ContactPredictionScoreResult::RefreshSortedIndices(){
  int n_col=this->matrix.size();
  std::vector <val_indices_pair> sorting_vec;
  for (int col1=0; col1<n_col; ++col1) {
      for (int col2=col1+1; col2<n_col; ++col2) {
        std::pair <int,int> indices (col1,col2);
        val_indices_pair vi (this->matrix[col1][col2],indices);
        sorting_vec.push_back(vi);
      }
    }
  std::sort(sorting_vec.begin(),sorting_vec.end(),pair_comparator);
  int n_indices=sorting_vec.size();
  std::vector <std::pair <int,int> > sorted_indices;
  for (int i=0; i<n_indices; ++i) {
    sorted_indices.push_back(sorting_vec[i].second);
  }
  this->sorted_indices=sorted_indices;
}

ContactPredictionScoreResult::ContactPredictionScoreResult(std::vector <std::vector <Real> > scores){
  this->matrix=scores;
  this->RefreshSortedIndices();
}

Real ContactPredictionScoreResult::GetScore(int i, int j){
    return this->matrix[i][j];
  }

void ContactPredictionScoreResult::SetScore(int i, int j, Real score){
    this->matrix[i][j]=score;
  }


ContactPredictionScoreResult CalculateContactSubstitutionScore(const AlignmentHandle& aln, int ref_seq_index, PairSubstWeightMatrix w){
  int n_col=aln.GetLength();
  int n_row=aln.GetCount();
  //We transform the alignment in a list of ints
  std::vector< std::vector <int> > ali(n_col, std::vector <int> (n_row) );
  int naa=w.aa_list.size();
  for (int i=0; i<n_col; ++i) {
    for (int j=0; j<n_row; ++j) {
      std::map<char,int>::iterator iter=w.aa_dict.find(aln[i][j]);
      if (iter!=w.aa_dict.end()){
        ali[i][j]=w.aa_dict[aln[i][j]];
      }
      else ali[i][j]=naa;
    }
  }
  std::vector< std::vector <Real> > contact_subst_score(n_col, std::vector <Real> (n_col,0) );
  for (int col1=0; col1<n_col; ++col1) {
    std::vector<int> c1=ali[col1];
    int i1=c1[ref_seq_index];
    if (i1==naa){
      for (int col2=col1+1; col2<n_col; ++col2) {
        contact_subst_score[col1][col2]=0.0;
        contact_subst_score[col2][col1]=0.0;
      }
      continue;
    }
    for (int col2=col1+1; col2<n_col; ++col2) {
      std::vector<int> c2=ali[col2];
      int i2=c2[ref_seq_index];
      if (i2==naa){
        contact_subst_score[col1][col2]=0.0;
        contact_subst_score[col2][col1]=0.0;
        continue;
        }
      Real score=0.0;
      for (int i=0 ; i<n_row; ++i){
        if (i==ref_seq_index) continue;
        int i3=c1[i], i4=c2[i];
        if (i3==naa || i4==naa){
          score+=0.0;
          continue;
        }
        score+=w.weights[i1][i2][i3][i4];
      }
    contact_subst_score[col1][col2]=score/(n_row-1.);
    contact_subst_score[col2][col1]=score/(n_row-1.);
    }
  }
  return ContactPredictionScoreResult(contact_subst_score);
}

ContactPredictionScoreResult CalculateContactScore(const AlignmentHandle& aln, ContactWeightMatrix weights){
  int n_col=aln.GetLength();
  int n_row=aln.GetCount();
  //We transform the alignment in a list of ints
  std::vector< std::vector <int> > ali(n_col, std::vector <int> (n_row) );
  int naa=weights.aa_list.size();
  for (int i=0; i<n_col; ++i) {
    for (int j=0; j<n_row; ++j) {
      std::map<char,int>::iterator iter=weights.aa_dict.find(aln[i][j]);
      if (iter!=weights.aa_dict.end()){
        ali[i][j]=weights.aa_dict[aln[i][j]];
      }
      else ali[i][j]=naa;
    }
  }

  std::vector< std::vector <Real> > mi(n_col, std::vector <Real> (n_col,0) );

  for (int col1=0; col1<n_col; ++col1) {
    std::vector<int> c1=ali[col1];
    for (int col2=col1+1; col2<n_col; ++col2) {
      float score=0.0;
      for (int i=0 ; i<n_row; ++i){
        int aa1=c1[i];
        int aa2=ali[col2][i];
        if (aa1==naa || aa2==naa) continue;
        score+=weights.weights[aa1][aa2];
      }
      mi[col1][col2]=score/n_row;
      mi[col2][col1]=score/n_row;
    }
  }
  return ContactPredictionScoreResult(mi);
}


ContactPredictionScoreResult CalculateMutualInformation(const AlignmentHandle& aln, ContactWeightMatrix weights,
                                  bool apc_correction,bool zpx_transformation, float small_number_correction){
  int n_col=aln.GetLength();
  int n_row=aln.GetCount();
  //We transform the alignment in a list of ints
  std::vector< std::vector <int> > ali(n_col, std::vector <int> (n_row) );
  int naa=weights.aa_list.size();
  for (int i=0; i<n_col; ++i) {
    for (int j=0; j<n_row; ++j) {
      std::map<char,int>::iterator iter=weights.aa_dict.find(aln[i][j]);
      if (iter!=weights.aa_dict.end()){
        ali[i][j]=weights.aa_dict[aln[i][j]];
      }
      else ali[i][j]=naa;
    }
  }

  std::vector< std::vector <Real> > mi(n_col, std::vector <Real> (n_col,0) );

  for (int col1=0; col1<n_col; ++col1) {
    std::vector<int> c1=ali[col1];
    for (int col2=col1+1; col2<n_col; ++col2) {
      std::vector<float> p1(naa,naa*small_number_correction),p2(naa,naa*small_number_correction);
      std::vector<std::vector<float> > p12(naa, std::vector<float> (naa,small_number_correction));
      for (int i=0 ; i<n_row; ++i){
        int aa1=c1[i];
        int aa2=ali[col2][i];
        if (aa1==naa || aa2==naa) continue;
        p1[aa1]++;
        p2[aa2]++;
        p12[aa1][aa2]++;
      }
      float n=naa*small_number_correction;
      for (int aa1=0; aa1!=naa; ++aa1){
        n+=p1[aa1];
      }
      float score=0.0;
      for (int aa1=0; aa1!=naa; ++aa1){
        if(p1[aa1]==0.0) continue;
        for (int aa2=0; aa2!=naa; ++aa2){
          if(p12[aa1][aa2]==0.0 || p2[aa2]==0.0) continue;
          score+=weights.weights[aa1][aa2]*p12[aa1][aa2]/n*std::log(n*p12[aa1][aa2]/(p1[aa1]*p2[aa2]));
        }
      }
      mi[col1][col2]=score;
      mi[col2][col1]=score;
    }
  }
  if (apc_correction==true){
    float average_mi=0.;
    std::vector <float> mi_col(n_col,0.);
    for (int col1=0; col1<n_col; ++col1) {
      for (int col2=0; col2<n_col; ++col2) {
        average_mi+=mi[col1][col2];
        mi_col[col1]+=mi[col1][col2];
      }
    }
    for (int col1=0; col1<n_col; ++col1) {
      for (int col2=col1+1; col2<n_col; ++col2) {
      mi[col1][col2]-=mi_col[col1]*mi_col[col2]/average_mi;
      mi[col2][col1]-=mi_col[col2]*mi_col[col1]/average_mi;
      }
    }
  }
  if (zpx_transformation==true){
    std::vector <float> mi_col(n_col,0.);
    for (int col1=0; col1<n_col; ++col1) {
      for (int col2=0; col2<n_col; ++col2) {
        mi_col[col1]+=mi[col1][col2];
      }
    }
    for (int col1=0; col1<n_col; ++col1) {
        mi_col[col1]/=n_col;
    }
    std::vector <float> mi_std_col(n_col,0.);
    for (int col1=0; col1<n_col; ++col1) {
      for (int col2=0; col2<n_col; ++col2) {
        mi_std_col[col1]+=(mi[col1][col2]-mi_col[col1])*(mi[col1][col2]-mi_col[col1]);
      }
    }
    for (int col1=0; col1<n_col; ++col1) {
        mi_std_col[col1]=sqrt(mi_std_col[col1]/n_col);
    }
    for (int col1=0; col1<n_col; ++col1) {
      for (int col2=col1+1; col2<n_col; ++col2) {
      mi[col1][col2]=(mi[col1][col2]-mi_col[col1])/mi_std_col[col1]+(mi[col1][col2]-mi_col[col2])/mi_std_col[col2];
      mi[col2][col1]=mi[col1][col2];
      }
    }
  }  
  return ContactPredictionScoreResult(mi);
}

}}}

