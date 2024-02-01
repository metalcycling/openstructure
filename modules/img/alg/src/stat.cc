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

#include "stat.hh"
#include "stat_accumulator.hh"

namespace ost { namespace img { namespace alg {

namespace {
  typedef std::pair<Real,Point> ValIndex;
}

template <typename T, class D>
void StatBase::VisitState(const ImageStateImpl<T,D>& isi)
{
  sum_=0.0;
  mean_=0.0;

  Real n = static_cast<Real>(isi.GetSize().GetVol());
  if(n==0.0){
    var_=0.0;
    std_dev_=0.0;
    min_=0.0;
    max_=0.0;
    maxpos_=Point(0,0,0),
    minpos_=Point(0,0,0),
    rms_=0.0;
    skewness_=0.0;
    kurtosis_=0.0;
    center_of_mass_=Vec3(0.0,0.0,0.0);
    return;
  }

  Vec3 sumcenter(0.0,0.0,0.0);
  ValIndex minindex(std::numeric_limits<Real>::max(),Point(0,0,0));
  ValIndex maxindex(-std::numeric_limits<Real>::max(),Point(0,0,0));
  min_ = std::numeric_limits<Real>::max();
  max_ = -std::numeric_limits<Real>::max();
  StatAccumulator<> acc;

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
        sumcenter+=Vec3(i,j,k)*val;
        acc(val);
      }
    }
  }
  min_=minindex.first;
  max_=maxindex.first;
  minpos_=minindex.second+isi.GetExtent().GetStart();
  maxpos_=maxindex.second+isi.GetExtent().GetStart();  
  var_=acc.GetVariance();
  std_dev_=acc.GetStandardDeviation();
  rms_=acc.GetRootMeanSquare();
  skewness_=acc.GetSkewness();
  kurtosis_= acc.GetKurtosis();
  sum_=acc.GetSum();
  mean_=acc.GetMean();
  if(sum_!=0.0){
    center_of_mass_=sumcenter/sum_+isi.GetExtent().GetStart().ToVec3();
  }else{
    center_of_mass_=Vec3(0.0,0.0,0.0);
  }
}

std::ostream& operator<<(std::ostream& o, const Stat& s)
{
  o << "mean=" << s.GetMean() << " ";
  o << "var=" << s.GetVariance() << " ";
  o << "std-dev=" << s.GetStandardDeviation() << " ";
  o << "min=" << s.GetMinimum() << " " << s.GetMinimumPosition() << " ";
  o << "max=" << s.GetMaximum() << " "<< s.GetMaximumPosition() << " ";
  o << "rms=" << s.GetRootMeanSquare() << " ";
  o << "skewness=" << s.GetSkewness() << " ";
  o << "kurtosis=" << s.GetKurtosis() << " ";
  o << "center of mass=" << s.GetCenterOfMass() << " ";
  return o;
}


}
namespace image_state {
  
template class TEMPLATE_DEF_EXPORT ImageStateNonModAlgorithm<alg::StatBase>;

}}} // ns
