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
#include <ost/seq/alg/contact_weight_matrix.hh>
#include <ost/seq/alg/default_contact_weight_matrix.hh>

namespace {
std::map<char,int> init_dict(std::vector <char> aal){
  std::map<char,int> aa_dict;
  const int aal_size = static_cast<int>(aal.size());
  for (int i=0;i!=aal_size;++i){
    aa_dict[aal[i]]=i;
  }
  return aa_dict;
}
}

namespace ost { namespace seq { namespace alg {

ContactWeightMatrix LoadConstantContactWeightMatrix(){
  std::vector <std::vector<Real> > w(21, std::vector<Real>(21,1.0));
  std::vector <char> aal(21);
  for (int i=0;i!=21;i++){
    aal[i]=RAW_CONTACT_WEIGHT_MATRIX_RES_LIST[i];
  }
  for (int i1=0;i1!=21;i1++){
    w[20][i1]=0.0;
    w[i1][20]=0.0;
  }
  return ContactWeightMatrix(w,aal);
}

ContactWeightMatrix LoadDefaultContactWeightMatrix(){
  std::vector <std::vector<Real> > w(21, std::vector<Real>(21));
  std::vector <char> aal(21);
  for (int i=0;i!=21;i++){
    aal[i]=RAW_CONTACT_WEIGHT_MATRIX_RES_LIST[i];
  }
  for (int i1=0;i1!=21;i1++){
    for (int i2=0;i2!=21;i2++){
      w[i1][i2]=RAW_CONTACT_WEIGHT_MATRIX[i1][i2];
    }
  }
  return ContactWeightMatrix(w,aal);
}

ContactWeightMatrix::ContactWeightMatrix(){
  std::vector <std::vector <Real> > w;
  std::vector <char> aal;
  this->aa_dict=init_dict(this->aa_list);
  this->weights=w;
}

ContactWeightMatrix::ContactWeightMatrix(std::vector <std::vector <Real> > w,std::vector <char> aal){
  this->aa_list=aal;
  this->naa_=this->aa_list.size();
  this->aa_dict=init_dict(this->aa_list);
  this->weights=w;
}


}}}

