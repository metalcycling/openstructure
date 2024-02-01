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
  generic fft functions for image state

  Author: Ansgar Philippsen
*/

#include <iostream>

#include "image_state_fft.hh"
#include "image_state_def.hh"

#include "fft_impl.hh"
using namespace ::img::fft_impl;

namespace ost { namespace img { namespace image_state { namespace fft {

/*
  BAD: too much implicit format specification in these routines
  and the fft implementation, especially regarding the half-complex
  format.
*/

namespace {

ImageStateBase* create_half_state(Extent ex, const PixelSampling& ps, Complex* ptr) 
{
  if(ex.GetDim()==1) {
    return new ComplexFrequencyHalf1DImageState(ex,ps,ptr);
  } else if(ex.GetDim()==2) {
    return new ComplexFrequencyHalf2DImageState(ex,ps,ptr);
  } else {
    return new ComplexFrequencyHalf3DImageState(ex,ps,ptr);
  }
}

} // anon ns

// out-of-place real to half-complex transform
template<>
ImageStateBase* FFT<Real, SpatialDomainPolicy>(ImageStateImpl<Real, SpatialDomainPolicy>* in_state, bool use_origin)
{
  Extent ex(in_state->GetInternalExtent());

  Complex* out_ptr = FFT_R2H(in_state->ValueHolder().GetData(), 
			     ex.GetSize());

  PixelSampling ps = in_state->GetSampling();
  ps.SetDomain(FREQUENCY);

  ImageStateBase* bp = create_half_state(ex,ps,out_ptr);

  if(use_origin) bp->AdjustPhaseOrigin(ex.GetStart());

  return bp;
}

// out-of-place half-complex to real transform
template<>
ImageStateBase* FFT<Complex, HalfFrequencyDomainPolicy>(ImageStateImpl<Complex, HalfFrequencyDomainPolicy>* in_state, bool use_origin)
{
  Extent ex(in_state->GetInternalExtent());

  if(use_origin) in_state->AdjustPhaseOrigin(-ex.GetStart());

  Real* out_ptr = FFT_H2R(in_state->ValueHolder().GetData(), 
			  ex.GetSize());

  if(use_origin) in_state->AdjustPhaseOrigin(ex.GetStart());

  PixelSampling ps = in_state->GetSampling();
  ps.SetDomain(SPATIAL);
  
  return new RealSpatialImageState(ex, ps, out_ptr);
}

// out-of-place complex spatial to complex frequency transform
template<>
ImageStateBase* FFT<Complex, SpatialDomainPolicy>(ImageStateImpl<Complex, SpatialDomainPolicy>* in_state, bool use_origin)
{
  Extent ex(in_state->GetInternalExtent());

  Complex* out_ptr = FFT_C2C(in_state->ValueHolder().GetData(),
			     ex.GetSize(), 
			     -1);

  PixelSampling ps = in_state->GetSampling();
  ps.SetDomain(FREQUENCY);

  ImageStateBase* bp = new ComplexFrequencyImageState(ex, ps, out_ptr);

  bp->AdjustPhaseOrigin(ex.GetStart());
  return bp;
}

// out-of-place complex frequency to complex spatial transform
template<>
ImageStateBase* FFT<Complex, FrequencyDomainPolicy>(ImageStateImpl<Complex, FrequencyDomainPolicy>* in_state, bool use_origin)
{
  Extent ex(in_state->GetInternalExtent());

  Complex* out_ptr = FFT_C2C(in_state->ValueHolder().GetData(),
			     ex.GetSize(), 
			     +1);

  PixelSampling ps = in_state->GetSampling();
  ps.SetDomain(SPATIAL);

  return new ComplexSpatialImageState(ex, ps, out_ptr);
}

#ifdef OSX
// workaround for bug on OSX that needs this symbol explicitely
template<>
ImageStateBase* FFT<Real, FrequencyDomainPolicy>(ImageStateImpl<Real, FrequencyDomainPolicy>* in_state, bool use_origin)
{
  throw InvalidFFT(String("FFT using state ") + typeid(ImageStateImpl<Real,FrequencyDomainPolicy>).name() + String(" is not possible"));
}
#endif

// default out-if-place variant
template<typename V, class DP>
ImageStateBase* FFT(ImageStateImpl<V,DP>* isi, bool use_origin) 
{
  throw InvalidFFT(String("FFT using state ") + typeid(ImageStateImpl<V,DP>).name() + String(" is not possible"));
}

// handle all in-place the same way
template<typename V, class DP>
ImageStateBase* FFTIP(ImageStateImpl<V,DP>* isi, bool use_origin)
{
  ImageStateBase* nb = FFT(isi, use_origin); // might throw InvalidFFT
  delete isi;
  return nb;
}

}}}} // namespaces
