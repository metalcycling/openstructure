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
  Authors: Andreas Schenk, Ansgar Philippsen
*/

#ifndef IMG_ALG_STAT_MIN_MAX_H
#define IMG_ALG_STAT_MIN_MAX_H

#include <iosfwd>

#include <ost/img/algorithm.hh>
#include <ost/img/image_state.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg { 

/*!
  collect statistical information about data
  this minimal version only calculated min and max

  Since this algorithm is implemented as a combined image stage visitor
  and algorithm, the main workhorse is this class StatMinMaxBase, which will
  act as the parent class of the actual algorithm class, Stat
*/

class DLLEXPORT_IMG_ALG StatMinMaxBase
{
 public:
  StatMinMaxBase():
    min_(0.0),
    max_(0.0),
    maxpos_(),
    minpos_()
  {}

  // image state algorithm interface
  template <typename T, class D>
  void VisitState(const ImageStateImpl<T,D>& isi);

  static String GetAlgorithmName() {return "StatMinMax";}

  // other interface
  
  Real GetMinimum() const {return min_;}
  Point GetMinimumPosition() const {return minpos_;}
  void SetMinimum(Real m) {min_=m;}
  Real GetMaximum() const {return max_;}
  Point GetMaximumPosition() const {return maxpos_;}
  void SetMaximum(Real m) {max_=m;}
protected:
  Real min_, max_;
  Point maxpos_,minpos_;
};

typedef ImageStateNonModAlgorithm<StatMinMaxBase> StatMinMax;

DLLEXPORT_IMG_ALG std::ostream& operator<<(std::ostream& o, const StatMinMax& s);

}

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateNonModAlgorithm<alg::StatMinMaxBase>)


  }} // ns

#endif

