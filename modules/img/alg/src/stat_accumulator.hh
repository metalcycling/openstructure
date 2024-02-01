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

#ifndef OST_STAT_ACCUMULATOR_HH
#define OST_STAT_ACCUMULATOR_HH

#include <limits>
#include <boost/math/special_functions/binomial.hpp>
#include <ost/base.hh>
#include <ost/message.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

namespace {

template<class D>
struct MinMax{
  static bool less_cmp_(const D& v1, const D& v2){return v1<v2;}
};
template<>
struct MinMax<Complex>{
  static bool less_cmp_(const Complex& v1, const Complex& v2){return abs(v1)<abs(v2);}
};

} //anon ns

template<unsigned int MAX_MOMENT=4,typename DATATYPE=Real>
class  StatAccumulator
{
public:
  StatAccumulator():
    sum_(0.0),
    sum2_(0.0),
    max_(-std::numeric_limits<DATATYPE>::max()),
    min_(std::numeric_limits<DATATYPE>::max()),
    m_(),
    w_(0.0),
    n_(0)
  {
    for(unsigned int i=0;i<MAX_MOMENT;++i){
      m_[i]=0.0;
    }

  }

  StatAccumulator(DATATYPE val, Real w=1.0):
    sum_(val),
    sum2_(val*val),
    max_(val),
    min_(val),
    m_(),
    w_(w),
    n_(1)
  {
    if(MAX_MOMENT>0){
      m_[0]=val;
      for(unsigned int i=1;i<MAX_MOMENT;++i){
        m_[i]=0.0;
      }
    }
  }

  void operator()(DATATYPE val, Real w=1.0)
  {
    *this+=StatAccumulator(val,w);
  }

  StatAccumulator<MAX_MOMENT,DATATYPE> operator+(const StatAccumulator<MAX_MOMENT,DATATYPE>& acc2) const
  {
    StatAccumulator<MAX_MOMENT,DATATYPE> acc(acc2);
    acc+=*this;
    return acc;
  }

  StatAccumulator<MAX_MOMENT,DATATYPE>& operator+=(const StatAccumulator<MAX_MOMENT,DATATYPE>& acc)
  {
    if(0.0>=w_){
      *this=acc;
      return *this;
    }
    if(0.0>=acc.w_){
      return *this;
    }
    n_+=acc.n_;
    Real wn=w_+acc.w_;

    sum_+=acc.sum_;
    sum2_+=acc.sum2_;
    max_=std::max<DATATYPE>(max_,acc.max_,MinMax<DATATYPE>::less_cmp_);
    min_=std::min<DATATYPE>(min_,acc.min_,MinMax<DATATYPE>::less_cmp_);
    if(MAX_MOMENT>0){
      DATATYPE delta=acc.m_[0]-m_[0];
      DATATYPE delta_w=delta/wn;
      if(MAX_MOMENT>2){
        DATATYPE w1w2_delta_w=w_*acc.w_*delta_w;
        DATATYPE w1w2_delta_wp=w1w2_delta_w*w1w2_delta_w;
        Real iw2=1.0/acc.w_;
        Real iw2pm1=iw2;
        Real miw1=-1.0/w_;
        Real miw1pm1=miw1;
        DATATYPE mn[MAX_MOMENT]; // only MAX_MOMENT-2 values needed but overdimensioned to allow compilation for MAX_MOMENT==0 (even though it gets kicked out in the end by the dead code elimination)
        for(unsigned int p=3;p<=MAX_MOMENT;++p){
          w1w2_delta_wp*=w1w2_delta_w;
          iw2pm1*=iw2;
          miw1pm1*=miw1;
          DATATYPE delta_wk=1.0;
          DATATYPE s=0.0;
          Real mw2k=1.0;
          Real w1k=1.0;
          for(unsigned int k=1;k<=p-2;++k){
            w1k*=w_;
            mw2k*=-acc.w_;
            delta_wk*=delta_w;
            s+=boost::math::binomial_coefficient<Real>(p, k)*(mw2k*m_[p-k-1]+w1k*acc.m_[p-k-1])*delta_wk;
          }
          mn[p-3]=acc.m_[p-1]+s+w1w2_delta_wp*(iw2pm1-miw1pm1);
        }
        for(unsigned int p=3;p<=MAX_MOMENT;++p){
          m_[p-1]+=mn[p-3];
        }
      }
      if(MAX_MOMENT>1){
        m_[1]+=acc.m_[1]+delta_w*delta*acc.w_*w_;
      }
      m_[0]+=delta_w*acc.w_;
      w_=wn;
    }
    return *this;
  }

  Real GetNumSamples() const
  {
    return n_;
  }

  DATATYPE GetMaximum() const
  {
    return max_;
  }

  DATATYPE GetMinimum() const
  {
    return min_;
  }

  Real GetWeight() const
  {
    return w_;
  }

  DATATYPE GetMoment(unsigned int i) const
  {
    if(1>i){
      throw Error("Invalid moment.");
    }
    if(MAX_MOMENT<i){
      throw Error("Moment was not calculated.");
    }
    return m_[i-1];
  }

  DATATYPE GetMean() const
  {
    if(MAX_MOMENT<1){
      throw Error("Mean was not calculated.");
    }
    return m_[0];
  }

  DATATYPE GetSum() const
  {
    return sum_;
  }

  DATATYPE GetVariance() const
  {
    if(MAX_MOMENT<2){
      throw Error("Variance was not calculated.");
    }
    if(n_==0.0){
      return 0.0;
    }
    return m_[1]/(w_);
  }

  DATATYPE GetStandardDeviation() const
  {
    return sqrt(GetVariance());
  }

  DATATYPE GetRootMeanSquare() const
  {
    if(n_==0.0){
      return 0.0;
    }
    return sqrt(sum2_/(w_));
  }

  DATATYPE GetSkewness() const
  {
    if(MAX_MOMENT<3){
      throw Error("Skewness was not calculated.");
    }
    if(m_[1]<1e-20){
      return 0.0;
    }else{
      return m_[2]/sqrt(m_[1]*m_[1]*m_[1]);
    }
  }

  DATATYPE GetKurtosis() const
  {
    if(MAX_MOMENT<4){
      throw Error("Kurtosis was not calculated.");
    }
    if(m_[1]<1e-20){
      return 0.0;
    }else{
      return w_*m_[3] / (m_[1]*m_[1]);
    }
  }

private:
  DATATYPE sum_;
  DATATYPE sum2_;
  DATATYPE max_;
  DATATYPE min_;
  DATATYPE m_[MAX_MOMENT];
  Real w_;
  unsigned int n_;
};

}}} //ns
#endif // OST_STAT_ACCUMULATOR_HH
