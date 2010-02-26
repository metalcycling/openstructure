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
  Author: Valerio Mariani
*/

#ifndef SMOOTH_MASK_IMAGE_HH_
#define SMOOTH_MASK_IMAGE_HH_

#include <ost/base.hh>
#include <ost/iplt/algorithm.hh>
#include <ost/iplt/alg/mask_image.hh>
#include <ost/iplt/alg/module_config.hh>

namespace ost { namespace iplt { namespace alg {


class DLLEXPORT_IPLT_ALG SmoothMaskImageBase
{
public:
  SmoothMaskImageBase();  // To make the compiler happy

  SmoothMaskImageBase(const MaskPtr& m, Real sigma, bool float_bg = false);

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is);

  static String GetAlgorithmName() {return "SmoothMask";}

private:
  MaskPtr mask_;
  Real sigma_;
  bool float_bg_;
};

typedef ImageStateModIPAlgorithm<SmoothMaskImageBase> SmoothMaskImage;

}

OST_IPLT_ALG_EXPLICIT_INST_DECL(class,ImageStateModIPAlgorithm<alg::SmoothMaskImageBase>)

  }} // namespaces

#endif /*SMOOTH_MASK_IMAGE_HH_*/
