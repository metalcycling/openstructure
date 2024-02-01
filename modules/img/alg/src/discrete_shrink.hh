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

#ifndef IMG_ALG_DISCRETE_SHRINK_HH
#define IMG_ALG_DISCRETE_SHRINK_HH

#include <ost/img/image.hh>
#include <ost/img/image_state.hh>

#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

/// \brief shrink image by integral amount
/// 
/// The algorithm performs a scaling of the orignal image by merging adjacent 
/// blocks of pixels. The block size can be passed in the constructor or
/// changed with SetBlocksize(). The spatial sampling of the image is updated, 
/// such that the original size is virtually unchanged (assuming, that the 
/// block size is a divisor of the image size)
struct DLLEXPORT_IMG_ALG DiscreteShrinkFnc {
  // initialize with blocksize
  DiscreteShrinkFnc(const Size& blocksize=Size(2,2,2)): bs_(blocksize) {}
  
  // set blocksize
  void SetBlocksize(const Size& bs) {bs_=bs;}
  Size GetBlocksize() const {return bs_;}
    
  template <typename T, class D>
  ImageStateBasePtr VisitState(const ImageStateImpl<T,D>& isi) const;

  static String GetAlgorithmName() {return "DiscreteShrink";}
protected:
  Size bs_;
};

typedef ImageStateConstModOPAlgorithm<DiscreteShrinkFnc> DiscreteShrink;

}

}} // ns

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateConstModOPAlgorithm<ost::img::alg::DiscreteShrinkFnc>)

#endif
