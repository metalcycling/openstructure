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
#include <limits> 
#include <cstring>
#include "entropy.hh"

#include <ost/seq/aligned_column.hh>
#include <ost/seq/aligned_column_iterator.hh>

namespace ost { namespace seq { namespace alg {

int OneLetterCodeToIndex(char olc, bool ignore_gaps) 
{
  if (olc >= 'A' && olc <= 'Z') return (olc-'A');
  if (olc >= 'a' && olc <= 'z') return (olc-'a');
  if (olc == '-' && !ignore_gaps) return 26;
  return -1;
}

std::vector<Real> ShannonEntropy(const AlignmentHandle& aln, bool ignore_gaps) {
  std::vector<Real> entropies(aln.GetLength(), 0);
  int aa_counts[27];
  int k = 0;
  for (AlignmentHandle::iterator 
       i = aln.begin(), e = aln.end(); i != e; ++i, ++k) {
    const AlignedColumn& col = *i;
    memset(aa_counts, 0, sizeof(int)*27);
    int counts = 0;
    for (int j=0; j<aln.GetCount(); ++j) {
      int index = OneLetterCodeToIndex(col[j], ignore_gaps);
      if (index==-1) {
        continue;
      }
      counts += 1;
      aa_counts[index] += 1;
    }
    Real entropy = 0.0;
    for (int j = 0; j < 27; ++j) {
      if (aa_counts[j]) {
        Real freq = static_cast<Real>(aa_counts[j])/counts;
        entropy-= log(freq)*freq;
      }
    }
    if (counts) 
      entropies[k] = entropy;
    else
      entropies[k] = std::numeric_limits<Real>::quiet_NaN();
  }
  return entropies;
}

}}}
