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
  Authors: Andreas Schenk,  Ansgar Philippsen
*/

#include <utility>
#include <cmath>
#include <iostream>

#include <ost/message.hh>
#include <ost/img/value_util.hh>

#include "stat_min_max.hh"

namespace ost { namespace img { namespace alg {

typedef std::pair<Real,Point> ValIndex;

template <typename T, class D>
void StatMinMaxBase::VisitState(const ImageStateImpl<T,D>& isi)
{
  Real n = isi.GetSize().GetVol();

  if(n==0.0)  return;

  ValIndex minindex(std::numeric_limits<Real>::max(),Point(0,0,0));
  ValIndex maxindex(-std::numeric_limits<Real>::max(),Point(0,0,0));
  min_ = std::numeric_limits<Real>::max();
  max_ = -std::numeric_limits<Real>::max();

  int wi=isi.GetSize()[0];
  int he=isi.GetSize()[1];
  int de=isi.GetSize()[2];
  for(int i=0;i<wi;++i) {
    for(int j=0;j<he;++j) {
      for(int k=0;k<de;++k) {
	      Real val=Val2Val<T,Real>(isi.Value(Index(i,j,k)));
        ValIndex vi(val,Point(i,j,k));
        minindex=std::min<ValIndex>(vi,minindex);
        maxindex=std::max<ValIndex>(vi,maxindex);
      }
    }
  }
  min_=minindex.first;
  max_=maxindex.first;
  minpos_=minindex.second+isi.GetExtent().GetStart();
  maxpos_=maxindex.second+isi.GetExtent().GetStart();  
}


std::ostream& operator<<(std::ostream& o, const StatMinMax& s)
{
  o << "min=" << s.GetMinimum() << " " << s.GetMinimumPosition() << " ";
  o << "max=" << s.GetMaximum() << " "<< s.GetMaximumPosition() << " ";
  return o;
}


}

namespace image_state {
  
template class TEMPLATE_DEF_EXPORT ImageStateNonModAlgorithm<alg::StatMinMaxBase>;  

}}} // ns
