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
#ifndef IMG_ALG_COMMON_LOCAL_SIGMA_THRES
#define IMG_ALG_COMMON_LOCAL_SIGMA_THRES

#include <ost/img/image_state.hh>
#include <ost/img/image_state/image_state_algorithm.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

/*
  Binary threshold based on local standard deviation; for each point,
  the standard deviation in a window of size*2+1 is checked against the
  given level, and if it surpasses that, the value in the resulting image
  is set to one, otherwise zero

  The output image is slightly smaller due to the edge being removed
*/

class DLLEXPORT_IMG_ALG LocalSigmaThresholdBase {
public:
  LocalSigmaThresholdBase();
  LocalSigmaThresholdBase(int size,Real level);

  template <typename T, class D>
  ImageStateBasePtr VisitState(const ImageStateImpl<T,D>& s) const;

  static String GetAlgorithmName() {return "LocalSigmaThreshold";}
private:
  int size_;
  Real level_;
};

typedef ImageStateConstModOPAlgorithm<LocalSigmaThresholdBase> LocalSigmaThreshold;


}

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateConstModOPAlgorithm<alg::LocalSigmaThresholdBase>)

  }} // ns


#endif
