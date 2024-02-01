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
#ifndef OST_SEC_STRUCT_HH
#define OST_SEC_STRUCT_HH

#include <vector>

#include <ost/base.hh>
#include <ost/geom/geom.hh>
#include <ost/message.hh>
#include <ost/mol/mol.hh>

namespace ost { namespace mol{ namespace alg{

inline Real DSSPHBondEnergy(const geom::Vec3& h_pos, const geom::Vec3& n_pos,
                            const geom::Vec3& c_pos, const geom::Vec3& o_pos) {
  Real on = 1.0/geom::Distance(o_pos,n_pos);
  Real ch = 1.0/geom::Distance(c_pos,h_pos);
  Real oh = 1.0/geom::Distance(o_pos,h_pos);
  Real cn = 1.0/geom::Distance(c_pos,n_pos);
  return 27.888 * (on+ch-oh-cn);
}

// Raw estimation of secondary structure
//
// This function is not intended for Python export, since the input has to be 
// prepared carefully. It basically estimates the secondary structure of a 
// stretch of amino acids based on the hydrogen bond pattern as described for 
// the dssp tool. 
//
// To define the hydrogen bonds you need to provide the vectors donor_for_one
// and donor_for two. The index of an acceptor residue appears in one of the
// two vectors if the corresponding hbond energy is < -0.5 
// (be aware of prolines that can't be donors!).
// There are two of those vectors, because for every residue we store the two
// lowest energy acceptors. If there is no acceptor available, the value
// at this position must be -1.
// The connected_to_next contains zeros and ones that defines, whether there
// is a peptide bond towards the next residue and the ca_positions vector
// is self explaining.
// As an additional feature you can also provide the according data for the
// full structure but only estimate the secondary structure for a small
// stretch, that gets defined by start_idx and size.
// For an example usage have a look at the AssignSecStruct functions.

String RawEstimateSS(const std::vector<geom::Vec3>& ca_positions, 
                     int start_idx, int size,
                     const std::vector<int>& donor_for_one,
                     const std::vector<int>& donor_for_two, 
                     const std::vector<int>& connected_to_next);

void PrepareSSData(const ost::mol::ResidueViewList& res_list,
                   std::vector<int>& res_indices,
                   std::vector<geom::Vec3>& ca_positions,
                   std::vector<int>& donor_for_one,
                   std::vector<int>& donor_for_two,
                   std::vector<int>& connected_to_next);

void AssignSecStruct(ost::mol::EntityView& ent);

void AssignSecStruct(ost::mol::EntityHandle& ent);

}}} //ns

#endif
