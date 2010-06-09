//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
#include <ost/seq/alg/conservation.hh>

namespace ost { namespace seq { namespace alg {

// taken from Armon et al., J. Mol. Biol. (2001) 307, 447-463
static float CHEM_DISSIM[][24]={
  {0.00, 1.33, 1.39, 2.22, 1.84, 1.45, 2.48, 3.26, 2.83, 3.48, 2.56, 3.27, 3.06, 
   0.86, 1.65, 1.63, 1.46, 2.24, 2.38, 3.34, 6.00, 6.00, 2.87, 3.15},
  {0.00, 0.06, 0.97, 0.56, 0.87, 1.92, 2.48, 1.80, 2.40, 2.15, 2.94, 2.90, 1.79, 
   2.70, 2.62, 2.36, 3.17, 3.12, 4.17, 6.00, 6.00, 2.20, 2.10,-1.00},
  {0.00, 0.91, 0.51, 0.90, 1.92, 2.46, 1.78, 2.37, 2.78, 3.54, 3.58, 2.76, 3.67, 
   3.60, 3.34, 4.14, 4.08, 5.13, 6.00, 6.00, 2.19, 2.08,-1.00,-1.00},
  {0.00, 0.85, 1.70, 2.48, 2.78, 1.96, 2.37, 2.78, 3.54, 3.58, 2.76, 3.67, 3.60, 
   3.34, 4.14, 4.08, 5.13, 6.00, 6.00, 2.63, 2.17,-1.00,-1.00,-1.00},
  {0.00, 0.89, 1.65, 2.06, 1.31, 1.87, 1.94, 2.71, 2.74, 2.15, 3.04, 2.95, 2.67, 
   3.45, 3.33, 4.38, 6.00, 6.00, 1.85, 1.59,-1.00,-1.00,-1.00,-1.00},
  {0.00, 1.12, 1.83, 1.40, 2.05, 1.32, 2.10, 3.03, 1.42, 2.25, 2.14, 1.86, 2.60, 
   2.45, 3.50, 6.00, 6.00, 1.47, 1.73,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.84, 0.99, 1.47, 0.32, 1.06, 1.30, 2.13, 2.70, 2.57, 2.30, 2.81, 2.48, 
   3.42, 6.00, 6.00, 0.42, 1.23,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.85, 0.90, 0.96, 1.14, 1.45, 2.97, 3.53, 3.39, 3.13, 3.59, 3.22, 4.08, 
   6.00, 6.00, 0.42, 0.88,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.65, 1.29, 1.84, 2.04, 2.76, 3.49, 3.37, 3.08, 3.70, 3.42, 4.39, 6.00, 
   6.00, 0.92, 0.33,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 1.72, 2.05, 2.34, 3.40, 4.10, 3.98, 3.69, 4.27, 3.95, 4.88, 6.00, 6.00, 
   1.18, 0.33,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.79, 0.82, 2.11, 2.59, 2.45, 2.19, 2.63, 2.27, 3.16, 6.00, 6.00, 0.64, 
   1.51,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.40, 2.70, 2.98, 2.84, 2.63, 2.85, 2.42, 3.11, 6.00, 6.00, 1.10, 1.95,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 2.43, 2.62, 2.49, 2.29, 2.47, 2.02, 2.72, 6.00, 6.00, 1.29, 2.19,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.91, 0.85, 0.62, 1.43, 1.52, 2.51, 6.00, 6.00, 2.55, 3.08,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.14, 0.41, 0.63, 0.94, 1.73, 6.00, 6.00, 3.11, 3.80,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.29, 0.61, 0.86, 1.72, 6.00, 6.00, 2.98, 3.68,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.82, 0.93, 1.89, 6.00, 6.00, 2.71, 3.39,-1.00,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 0.48, 0.11, 6.00, 6.00, 3.20, 3.99,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 1.06, 6.00, 6.00, 2.85, 3.67,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 6.00, 6.00, 3.75, 4.64,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.50, 6.00, 6.00, 6.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 6.00, 6.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00, 1.05,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00},
  {0.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,
  -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00}
};

float PhysicoChemicalDissim(char c1, char c2)
{
  static int indices[]={2, 23, 0, 9, 7, 17, 3, 10, 15, -1, 
             11, 14, 16, 8, -1, 1, 6, 12, 4, 5, 
             -1, 13, 19, 21, 18, 22};
  int idx_a=(c1=='-' || c1<'A' || c1>'Z') ? 20 : indices[c1-'A'];
  int idx_b=(c2=='-' || c2<'A' || c2>'Z') ? 20 : indices[c2-'A'];  
  assert(idx_a>=0);
  assert(idx_b>=0);
  if (idx_a<0) {
    idx_a=20;
  }
  if (idx_b<0) {
    idx_b=20;
  }
  float s=0.0;
  if (idx_a>idx_b)
    s=CHEM_DISSIM[idx_b][idx_a-idx_b];
  else
    s=CHEM_DISSIM[idx_a][idx_b-idx_a];
  assert(s>=0.0);
  return s;
}

std::vector<Real> Conservation(const AlignmentHandle& aln, bool assign, 
                               const String& prop)
{
  std::vector<Real> cons(aln.GetLength(), 0.0);
  int comb=aln.GetCount()*(aln.GetCount()-1)/2;
  for (int col=0; col<aln.GetLength(); ++col) {
    float score=0.0;
    AlignedColumn c=aln[col];
    for (int i=0; i<aln.GetCount(); ++i) {
      for (int j=i+1; j<aln.GetCount(); ++j) {
        score+=PhysicoChemicalDissim(c[i], c[j]);
      }
    }
    score=1.0-score/(6.0*comb);
    cons[col]=score;
    if (assign) {
      for (int i=0; i<aln.GetCount(); ++i) {
        if (c[i]!='-') {
          mol::ResidueView r=c.GetResidue(i);
          if (r.IsValid()) {
            r.SetFloatProp(prop, score);
          }
        }
      }
    }
  }
  return cons;
}

}}}

