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
#ifndef OST_SEQ_CONTACT_WEIGHT_MATRIX_HH
#define OST_SEQ_CONTACT_WEIGHT_MATRIX_HH

#include <map>
#include <vector>

#include <ost/base.hh>
#include <ost/seq/alg/module_config.hh>

/*
  Author: Niklaus Johner
 */
namespace ost { namespace seq { namespace alg {

  
/// \brief matrix to weight pairs of amino-acids.
struct DLLEXPORT_OST_SEQ_ALG ContactWeightMatrix {
  std::vector <std::vector <Real> > weights;
  std::vector <char> aa_list;
  std::map<char,int> aa_dict;
  int naa_;

  ContactWeightMatrix();
  ContactWeightMatrix(std::vector <std::vector <Real> >,std::vector <char>);
};

/// \brief statistical potential matrix containing interaction pseudo energies
ContactWeightMatrix DLLEXPORT_OST_SEQ_ALG LoadDefaultContactWeightMatrix();
/// \brief weight of 1 for all amino-acid pairs and 0 for gaps.
ContactWeightMatrix DLLEXPORT_OST_SEQ_ALG LoadConstantContactWeightMatrix();

}}}

#endif