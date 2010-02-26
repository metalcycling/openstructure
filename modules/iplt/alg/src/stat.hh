//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
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
  Authors: Andreas Schenk, Ansgar Philippsen
*/

#ifndef IPLT_ALG_STAT_H
#define IPLT_ALG_STAT_H

#include <iosfwd>

#include <ost/iplt/algorithm.hh>
#include <ost/iplt/image_state.hh>
#include <ost/iplt/alg/module_config.hh>

namespace ost { namespace iplt { namespace alg { 

/*!
  collect statistical information about data

  Since this algorithm is implemented as a combined image stage visitor
  and algorithm, the main workhorse is this class StatBase, which will
  act as the parent class of the actual algorithm class, Stat
*/

class DLLEXPORT_IPLT_ALG StatBase
{
 public:
  StatBase():
    mean_(0.0),
    var_(0.0),
    std_dev_(0.0),
    sum_(0.0),
    min_(0.0),
    max_(0.0),
    maxpos_(),
    minpos_(),
    rms_(0.0),
    skewness_(0.0),
    kurtosis_(0.0),
    center_of_mass_(0.0,0.0,0.0)
  {}

  // image state algorithm interface
  template <typename T, class D>
  void VisitState(const ImageStateImpl<T,D>& isi);

  void VisitFunction(const Function& f);

  static String GetAlgorithmName() {return "Stat";}

  // other interface
  
  Real GetMean() const {return mean_;}
  void SetMean(Real m) {mean_=m;}
  Real GetMinimum() const {return min_;}
  Point GetMinimumPosition() const {return minpos_;}
  void SetMinimum(Real m) {min_=m;}
  Real GetMaximum() const {return max_;}
  Point GetMaximumPosition() const {return maxpos_;}
  void SetMaximum(Real m) {max_=m;}
  Real GetSum() const {return sum_;}
  void SetSum(Real s) {sum_=s;}
  Real GetVariance() const {return var_;}
  void SetVariance(Real v) {var_=v;}
  Real GetStandardDeviation() const {return std_dev_;}
  void SetStandardDeviation(Real s) {std_dev_=s;} 
  Real GetRootMeanSquare() const {return rms_;}
  Real GetSkewness() const {return skewness_;}
  Real GetKurtosis() const {return kurtosis_;}
  Vec3 GetCenterOfMass() const {return center_of_mass_;}
protected:
  Real mean_, var_, std_dev_;
  Real sum_, min_, max_;
  Point maxpos_,minpos_;
  Real rms_,skewness_,kurtosis_;
  Vec3 center_of_mass_;
};

typedef ImageStateNonModAlgorithm<StatBase> Stat;

DLLEXPORT_IPLT_ALG std::ostream& operator<<(std::ostream& o, const Stat& s);
}

OST_IPLT_ALG_EXPLICIT_INST_DECL(class,ImageStateNonModAlgorithm<alg::StatBase>)

}} // ns

#endif

