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
  Author: Marco Biasini, Juergen Haas, Gerardo Tauriello
 */

#include <ost/mol/mol.hh>
#include <ost/seq/alg/distance_map.hh>

namespace ost { namespace seq { namespace alg {

DistanceMapPtr CreateDistanceMap(const seq::AlignmentHandle& alignment) {

  // setup and check
  if (alignment.GetCount() <= 1) {
    throw IntegrityError("At least two sequences required!");
  }
  seq::ConstSequenceHandle ref_seq = alignment.GetSequence(0);
  String ref_seq_str = ref_seq.GetString();
  uint len = ref_seq.GetGaplessString().size();
  DistanceMapPtr dist_map(new DistanceMap(len, alignment.GetCount() - 1));

  // parse all non-ref sequences
  for (int j = 1; j < alignment.GetCount(); ++j) {
    seq::ConstSequenceHandle s = alignment.GetSequence(j);
    if (!s.HasAttachedView()) {
      throw IntegrityError("All non-ref sequences must have an attached view!");
    }
    // parse all residues in ref seq.
    for (uint x = 0; x < ref_seq_str.size(); ++x) {
      // skip gaps
      if (ref_seq_str[x] == '-') continue;
      // find CA and skip gaps and weird residues
      mol::ResidueView res_a = s.GetResidue(x);
      if (!res_a) continue;
      mol::AtomView ca_a = res_a.FindAtom("CA");
      if (!ca_a) continue;
      // compare with all others (considering symmetry)
      for (uint y = x+1; y < ref_seq_str.size(); ++y) {
        // skip gaps
        if (ref_seq_str[y]=='-') continue;
        // find CA and skip gaps and weird residues
        mol::ResidueView res_b = s.GetResidue(y);
        if (!res_b) continue;
        mol::AtomView ca_b = res_b.FindAtom("CA");
        if (!ca_b) continue;
        // add it to the map
        dist_map->AddDistance(ref_seq.GetResidueIndex(x), 
                              ref_seq.GetResidueIndex(y), 
                              geom::Distance(ca_a.GetPos(), ca_b.GetPos()),
                              j);
      }
    }
  }
  return dist_map;
} 
 
}}}
