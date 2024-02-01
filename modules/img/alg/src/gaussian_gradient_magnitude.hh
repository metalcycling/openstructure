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
  Authors: Giani Signorell, Ansgar Philippsen, Valerio Mariani
*/

#ifndef IMG_ALG_FILTER_GAUSSIAN_GRADIENT_MAGNITUDE_H
#define IMG_ALG_FILTER_GAUSSIAN_GRADIENT_MAGNITUDE_H

#include <ost/img/algorithm.hh>
#include <ost/img/image_state.hh>
#include <ost/img/image_state/image_state_visitor.hh>
#include <ost/message.hh>

#include "line_iterator.hh"

namespace ost { namespace img { namespace alg {

//! First Derivative Gaussian Filter
/*!
  Applies a gaussian filter and at the same time a 
  first derivative filter to the supplied image. 
  Implemented after 
  I.T.Young, L.J. van Vliet,"Recursive implementation of the 
  Gaussian filter",Signal Processing, 44(1995), 139-151
*/
class DLLEXPORT_IMG_ALG GaussianGradientMagnitudeFilterBase
{
public:
  //! Initialization takes smoothing radius
  GaussianGradientMagnitudeFilterBase(Real sigma=1.0);

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& s);

  static String GetAlgorithmName() {return "GaussianGradientMagnitudeFilter";}

  //! Set smoothing radius  
  void SetSigma(Real sigma) {
    calcBVals( calcQ(sigma) );
  }
  //! Set q value directly
  void SetQ(Real q) {
    calcBVals( q );
  }
		
protected:
  Real b0_, b1_, b2_, b3_, bb_;
		
private:
  Real calcQ( Real sigma );
  void calcBVals( Real q );

  template <typename T, class D>
  void lineTransformFW(ImageStateImpl<T,D>& s, ExtentIterator& it);
  
  template <typename T, class D>
  void lineTransformBW(ImageStateImpl<T,D>& s, ExtentIterator& it);
};

typedef ImageStateModIPAlgorithm<GaussianGradientMagnitudeFilterBase> GaussianGradientMagnitudeFilter;

}

OST_IMG_ALG_EXPLICIT_INST_DECL(class,ImageStateModIPAlgorithm<alg::GaussianGradientMagnitudeFilterBase>)

}} // namespaces

#endif
