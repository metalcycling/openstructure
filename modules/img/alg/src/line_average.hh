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
#ifndef IMG_ALG_LINE_AVERAGE_HH
#define IMG_ALG_LINE_AVERAGE_HH

#include <ost/img/extent.hh>
#include <ost/img/image_state.hh>
#include <ost/img/image_state/image_state_algorithm.hh>
#include <ost/dllexport.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

//! computes mean value of columns or lines
class DLLEXPORT_IMG_ALG LineAverageBase
{
 public:

  LineAverageBase(): geometry_(0) {};
  LineAverageBase(unsigned int geometry): geometry_(geometry) {};

  // image state algorithm interface
  template <typename T, class D>
  img::ImageStateBasePtr VisitState(const ImageStateImpl<T,D>& s) const;

  static String GetAlgorithmName() {return "LineAverage";}

  void SetGeometry (unsigned int geometry) { geometry_ = geometry; }
  unsigned int GetGeometry () { return geometry_; }

private:

   int geometry_;

};

typedef image_state::ImageStateConstModOPAlgorithm<alg::LineAverageBase> LineAverage;

}
OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateConstModOPAlgorithm<alg::LineAverageBase>)
}} // ns

#endif // IPLT_ALG_LINE_AVERAGE_HH
