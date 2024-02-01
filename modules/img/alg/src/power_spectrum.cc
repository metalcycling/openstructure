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

#include <ost/message.hh>
#include <ost/dyn_cast.hh>
#include "fft.hh"
#include "power_spectrum.hh"

namespace ost { namespace img { namespace alg {

PowerSpectrumFnc::PowerSpectrumFnc()
{}

template <>
ImageStateBasePtr PowerSpectrumFnc::VisitState<Real,SpatialDomain>(const RealSpatialImageState& in_state) const
{
  ImageStateBasePtr fft_state_bp=FFT().VisitState(in_state);
  boost::shared_ptr<ComplexHalfFrequencyImageState> fft_state=dyn_cast<ComplexHalfFrequencyImageState>(fft_state_bp);
  if(!fft_state) {
    throw Error("unexpected dynamic_cast failure");
  }

  boost::shared_ptr<RealSpatialImageState> out_state(new RealSpatialImageState(fft_state->GetLogicalExtent(),fft_state->GetSampling()));

  /*
    this could be faster, but involved complicated half_complex conversion
    that is best dealt with using the Complex() method
  */
  for(ExtentIterator it(out_state->GetExtent());!it.AtEnd();++it) {
    out_state->Value(it)=std::abs(fft_state->GetComplex(it));
  }

  return out_state;
}

template <>
ImageStateBasePtr PowerSpectrumFnc::VisitState<Complex,SpatialDomain>(const ComplexSpatialImageState& in_state) const
{
  ImageStateBasePtr fft_state_bp=FFT().VisitState(in_state);
  boost::shared_ptr<ComplexFrequencyImageState> fft_state=dyn_cast<ComplexFrequencyImageState>(fft_state_bp);
  if(!fft_state) {
    throw Error("unexpected dynamic_cast failure");
  }

  boost::shared_ptr<RealSpatialImageState> out_state(new RealSpatialImageState(fft_state->GetExtent(),fft_state->GetSampling()));

  for(ExtentIterator it(out_state->GetExtent());!it.AtEnd();++it) {
    out_state->Value(it)=std::abs(fft_state->Value(it));
  }
  
  return out_state;
}

template <>
ImageStateBasePtr PowerSpectrumFnc::VisitState<Complex,FrequencyDomain>(const ComplexFrequencyImageState& in_state) const
{
  // simply convert to real spatial image state
  boost::shared_ptr<RealSpatialImageState> out_state(new RealSpatialImageState(in_state.GetExtent(),in_state.GetSampling()));
  
  for(ExtentIterator it(out_state->GetExtent());!it.AtEnd();++it) {
    out_state->Value(it)=std::abs(in_state.Value(it));
  }
  
  return out_state;
}

template <>
ImageStateBasePtr PowerSpectrumFnc::VisitState<Complex,HalfFrequencyDomain>(const ComplexHalfFrequencyImageState& in_state) const
{
  boost::shared_ptr<RealSpatialImageState> out_state(new RealSpatialImageState(in_state.GetLogicalExtent(),in_state.GetSampling()));

  /*
    see note above about speed of this
  */
  for(ExtentIterator it(out_state->GetExtent());!it.AtEnd();++it) {
    out_state->Value(it)=std::abs(in_state.GetComplex(it));
  }

  return out_state;
}

template <typename V, class D>
ImageStateBasePtr PowerSpectrumFnc::VisitState(const ImageStateImpl<V,D>& isi) const
{
  throw Error("PowerSpectrum not supported for this image state");
  return ImageStateBasePtr();
}

}

namespace image_state {

template class TEMPLATE_DEF_EXPORT ImageStateConstModOPAlgorithm<alg::PowerSpectrumFnc>;

}}} //

