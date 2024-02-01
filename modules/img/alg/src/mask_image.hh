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
  Author: Andreas Schenk
*/

#ifndef MASK_IMAGE_HH_
#define MASK_IMAGE_HH_

#include <ost/base.hh>
#include <ost/img/algorithm.hh>
#include <ost/img/mask.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

class DLLEXPORT_IMG_ALG MaskImageBase
{
public:
  MaskImageBase(); // To make the compiler happy

  MaskImageBase(const MaskPtr& m, bool float_bg = false);
 
  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is) const;

  static String GetAlgorithmName() {return "MaskImage";}

private:
  MaskPtr mask_;
  bool float_bg_;
};

typedef ImageStateConstModIPAlgorithm<MaskImageBase> MaskImage;


}


}} // namespaces

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateConstModIPAlgorithm<ost::img::alg::MaskImageBase>)

#endif /*MASK_IMAGE_HH_*/

