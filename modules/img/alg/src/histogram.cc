//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  Author: Ansgar Philippsen
*/

#include <ost/message.hh>

#include "stat.hh"
#include "histogram.hh"

namespace ost { namespace img { namespace alg {

HistogramBase::HistogramBase():
  bin_count_(10),
  min_(0.0),
  max_(1.0),
  cfac_(1.0),
  bins_(10,0)
{
}

HistogramBase::HistogramBase(int bin_count, Real mn, Real mx):
  bin_count_(bin_count),
  min_(mn),
  max_(mx),
  cfac_(1.0),
  bins_(bin_count,0)
{
  if(bin_count_<2) throw HistogramError("bin_count must be larger than 1");
  if(max_==min_) throw HistogramError("maximummust not be equal to minimum");

  cfac_=Real(bin_count-1)/(max_-min_);
}


const HistogramBase::Bins& HistogramBase::GetBins() const {return bins_;}

template <typename T, class D>
void HistogramBase::VisitState(const ImageStateImpl<T,D>& isi)
{
  bins_=Bins(bin_count_,0);
  for(const T* ptr = isi.Data().GetData(); ptr<isi.Data().GetEnd(); ++ptr) {
    Real val=Val2Val<T,Real>(*ptr);
    val = std::max(min_,val);
    val=std::min(max_,val);
    int bin = int(floor((val-min_)*cfac_));
    bins_[bin]++;
  }
}

}

namespace image_state {

template class TEMPLATE_DEF_EXPORT ImageStateNonModAlgorithm<alg::HistogramBase>;

}}} // namespaces
