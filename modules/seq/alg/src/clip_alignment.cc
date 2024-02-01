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
  Author: Gerardo Tauriello, Juergen Haas
 */

#include <vector>
#include <ost/seq/alg/clip_alignment.hh>

namespace ost { namespace seq { namespace alg {

int ClipAlignment(AlignmentHandle& aln, uint n_seq_thresh,
                  bool set_offset, bool remove_empty) {
  // count entries per column
  const int num_residues = aln.GetLength();
  const int num_seqs = aln.GetCount();
  std::vector<uint> counts(num_residues, 0);
  // parse all sequences
  for (int i_s = 0; i_s < num_seqs; ++i_s) {
    String s = aln.GetSequence(i_s).GetString();
    // only consider non-gaps
    const size_t from = s.find_first_not_of("-");
    if (from == String::npos) continue;
    const size_t last = s.find_last_not_of("-");
    for (size_t i = from; i <= last; ++i) {
      if (s[i] != '-') ++counts[i];
    }
  }
  
  // find clip_start and clip_end
  int clip_start = 0;
  while (clip_start < num_residues && counts[clip_start] < n_seq_thresh) {
    ++clip_start;
  }
  // check if we clipped away everything -> ABORT
  if (clip_start == num_residues) {
    if (remove_empty) aln = seq::CreateAlignment();
    return -1;
  }
  int clip_end = num_residues - 1;
  while (clip_end > clip_start && counts[clip_end] < n_seq_thresh) {
    --clip_end;
  }

  // update offsets
  if (set_offset) {
    for (int i_s = 0; i_s < num_seqs; ++i_s) {
      // only relevant if a view attached...
      if (!aln.GetSequence(i_s).HasAttachedView()) continue;
      // count non-gapped part to be cut
      String s = aln.GetSequence(i_s).GetString();
      const size_t from = s.find_first_not_of("-");
      if (from == String::npos) continue;
      uint cut_res = 0;
      for (int i = from; i < clip_start; ++i) {
        if (s[i] != '-') ++cut_res;
      }
      // update offset
      aln.SetSequenceOffset(i_s, aln.GetSequenceOffset(i_s) + cut_res);
    }
  }

  // cut alignment
  if (clip_end < num_residues - 1) aln.Cut(clip_end + 1, num_residues);
  if (clip_start > 0) aln.Cut(0, clip_start);

  // remove empty sequences
  if (remove_empty) {
    int i_s = 0;
    while (i_s < aln.GetCount()) {
      if (aln.GetSequence(i_s).GetFirstNonGap() == -1) {
        aln.RemoveSequence(i_s);
      } else {
        ++i_s;
      }
    }
  }

  return clip_start;
}

}}}
