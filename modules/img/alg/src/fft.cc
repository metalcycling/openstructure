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


#include <sstream>
#include <cassert>

#include <boost/shared_ptr.hpp>

#include <fftw3.h>
#if OST_INFO_ENABLED
#include <QThread>
#define IDEAL_NUMBER_OF_THREADS() QThread::idealThreadCount()
#else
#define IDEAL_NUMBER_OF_THREADS() 1
#endif
#include <ost/message.hh>
#include <ost/img/value_util.hh>
#include <ost/img/image_state/image_state_def.hh>
#include "fft.hh"
#include <ost/img/alg/fftw_helper.hh>


namespace ost { namespace img { namespace alg {
namespace {

template <class D>
String Domain2String();

template<>
String Domain2String<image_state::SpatialDomain>() {return "SPATIAL";}

template<>
String Domain2String<image_state::FrequencyDomain>() {return "FREQUENCY";}

template <class D>
String Domain2String() {return "HALF_FREQUENCY";}

int half_plus_one(int n) {
  return (n&0x1) ? (n-1)/2+1 : n/2+1;
}

} // anon ns

FFTFnc::FFTFnc(): ori_flag_(false)
{
  OST_FFTW_fftw_init_threads();
}
FFTFnc::FFTFnc(bool f): ori_flag_(f)
{
  OST_FFTW_fftw_init_threads();
}

FFTFnc::~FFTFnc()
{
  OST_FFTW_fftw_cleanup();
}

// real spatial -> complex half-frequency
template <>
ImageStateBasePtr FFTFnc::VisitState<Real,SpatialDomain>(const RealSpatialImageState& in_state) const
{
  int rank,n[3]; // for actual fftw call

  Size in_size = in_state.GetExtent().GetSize();
  n[0] = in_size.GetWidth();
  n[1] = in_size.GetHeight();
  n[2] = in_size.GetDepth();
  rank=in_size.GetDim();

  PixelSampling ps=in_state.GetSampling();
  ps.SetDomain(FREQUENCY);

  boost::shared_ptr<ComplexHalfFrequencyImageState> out_state(new ComplexHalfFrequencyImageState(in_size,ps));

  out_state->SetSpatialOrigin(in_state.GetSpatialOrigin());
  out_state->SetAbsoluteOrigin(in_state.GetAbsoluteOrigin());

  assert(sizeof(OST_FFTW_fftw_complex)==sizeof(Complex));
  OST_FFTW_fftw_complex* fftw_out = reinterpret_cast<OST_FFTW_fftw_complex*>(out_state->Data().GetData());
  Real* fftw_in = const_cast<Real*>(in_state.Data().GetData());
  OST_FFTW_fftw_plan_with_nthreads(std::max<int>(1,IDEAL_NUMBER_OF_THREADS()));
  OST_FFTW_fftw_plan plan = OST_FFTW_fftw_plan_dft_r2c(rank,n,fftw_in,fftw_out,FFTW_ESTIMATE|FFTW_PRESERVE_INPUT);
  OST_FFTW_fftw_execute(plan);
  OST_FFTW_fftw_destroy_plan(plan);

  if(ori_flag_) out_state->AdjustPhaseOrigin(in_state.GetSpatialOrigin());

  return out_state;
}

// complex half-frequency -> real spatial
template <>
ImageStateBasePtr FFTFnc::VisitState<Complex,HalfFrequencyDomain>(const ComplexHalfFrequencyImageState& in_state) const
{
  int rank,n[3]; // for actual fftw call

  Size in_size=in_state.GetLogicalExtent().GetSize();
  n[0]=in_size.GetWidth();
  n[1]=in_size.GetHeight();
  n[2]=in_size.GetDepth();
  // copy
  ComplexHalfFrequencyImageState tmp_state(in_state);
  // correct phase origin
  if (ori_flag_) tmp_state.AdjustPhaseOrigin(-tmp_state.GetSpatialOrigin());

  switch (in_size.GetDim()) {
  case 1: 
    rank=1; 
    break;
  case 2: 
    rank=2; 
    // complete zero line
    for(int i=1;i<half_plus_one(n[0]);++i){
        tmp_state.Data().Value(Index(n[0]-i,0,0))=conj(tmp_state.Data().Value(Index(i,0,0)));
    }
    break;
  case 3:
    rank=3;
    // complete zero line
    for(int i=1;i<half_plus_one(n[0]);++i){
        tmp_state.Data().Value(Index(n[0]-i,0,0))=conj(tmp_state.Data().Value(Index(i,0,0)));
    }
    // complete rest of zero plane
    for(int i=0;i<n[0];++i){
      for(int j=1;j<half_plus_one(n[1]);++j){
        tmp_state.Data().Value(Index((n[0]-i)%n[0],n[1]-j,0))=conj(tmp_state.Data().Value(Index(i,j,0)));
      }
    }
    break;
  default:
    throw(FFTException("unexpected dimension in FFT C2R"));
  }


  PixelSampling ps=in_state.GetSampling();
  ps.SetDomain(SPATIAL);

  Size out_size = in_state.GetLogicalExtent().GetSize();
  boost::shared_ptr<RealSpatialImageState> out_state(new RealSpatialImageState(out_size,ps ));
  out_state->SetSpatialOrigin(in_state.GetSpatialOrigin());
  out_state->SetAbsoluteOrigin(in_state.GetAbsoluteOrigin());


  Real* out_ptr = out_state->Data().GetData();
  Real* fftw_out = reinterpret_cast<Real*>(out_ptr);
  
  assert(sizeof(OST_FFTW_fftw_complex)==sizeof(Complex));
  OST_FFTW_fftw_plan_with_nthreads(std::max<int>(1,IDEAL_NUMBER_OF_THREADS()));
  Complex* in_ptr = tmp_state.Data().GetData();
  OST_FFTW_fftw_complex* fftw_in =reinterpret_cast<OST_FFTW_fftw_complex*>(in_ptr);

  OST_FFTW_fftw_plan plan = OST_FFTW_fftw_plan_dft_c2r(rank,n,fftw_in,fftw_out,FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

  OST_FFTW_fftw_execute(plan);
  OST_FFTW_fftw_destroy_plan(plan);


  Real fac = 1.0/static_cast<Real>(out_state->GetSize().GetVolume());
  for(Real* ptr = out_state->Data().GetData(); ptr<out_state->Data().GetEnd(); ++ptr) {
    (*ptr)*=fac;
  }

  return out_state;
}


// complex spatial -> complex frequency
template <>
ImageStateBasePtr FFTFnc::VisitState<Complex,SpatialDomain>(const ComplexSpatialImageState& in_state) const
{
  assert(sizeof(OST_FFTW_fftw_complex)==sizeof(Complex));

  Size size=in_state.GetExtent().GetSize();
  PixelSampling ps=in_state.GetSampling();
  ps.SetDomain(FREQUENCY);
  int rank = size.GetDim();
  int n[3] = {static_cast<int>(size[0]),static_cast<int>(size[1]),static_cast<int>(size[2])};
  int dir = FFTW_FORWARD;
  OST_FFTW_fftw_complex* fftw_in = const_cast<OST_FFTW_fftw_complex*>(reinterpret_cast<const OST_FFTW_fftw_complex*>(in_state.Data().GetData()));

  boost::shared_ptr<ComplexFrequencyImageState> out_state(new ComplexFrequencyImageState(size,ps));
  out_state->SetSpatialOrigin(in_state.GetSpatialOrigin());
  out_state->SetAbsoluteOrigin(in_state.GetAbsoluteOrigin());
  OST_FFTW_fftw_complex* fftw_out = reinterpret_cast<OST_FFTW_fftw_complex*>(out_state->Data().GetData());

  // out of place transform
  OST_FFTW_fftw_plan_with_nthreads(std::max<int>(1,IDEAL_NUMBER_OF_THREADS()));
  OST_FFTW_fftw_plan plan = OST_FFTW_fftw_plan_dft(rank,n,fftw_in, fftw_out,dir,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
  OST_FFTW_fftw_execute(plan);
  OST_FFTW_fftw_destroy_plan(plan);

  if(ori_flag_) out_state->AdjustPhaseOrigin(in_state.GetSpatialOrigin());

  return out_state;
}

// complex frequency -> complex spatial
template <>
ImageStateBasePtr FFTFnc::VisitState<Complex,FrequencyDomain>(const ComplexFrequencyImageState& in_state) const
{
  assert(sizeof(OST_FFTW_fftw_complex)==sizeof(Complex));

  Size size=in_state.GetExtent().GetSize();
  PixelSampling ps=in_state.GetSampling();
  ps.SetDomain(SPATIAL);
  int rank = size.GetDim();
  int n[3] = {static_cast<int>(size[0]),static_cast<int>(size[1]),static_cast<int>(size[2])};
  int dir = FFTW_BACKWARD;

  boost::shared_ptr<ComplexSpatialImageState> out_state(new ComplexSpatialImageState(size,ps));
  out_state->SetSpatialOrigin(in_state.GetSpatialOrigin());
  out_state->SetAbsoluteOrigin(in_state.GetAbsoluteOrigin());
  OST_FFTW_fftw_complex* fftw_out = reinterpret_cast<OST_FFTW_fftw_complex*>(out_state->Data().GetData());

  if(ori_flag_){
    // copy in state
    ComplexFrequencyImageState tmp_state(in_state);
    tmp_state.AdjustPhaseOrigin(-in_state.GetSpatialOrigin());
    OST_FFTW_fftw_complex* fftw_in = reinterpret_cast<OST_FFTW_fftw_complex*>(tmp_state.Data().GetData());
    // out of place transform
    OST_FFTW_fftw_plan_with_nthreads(std::max<int>(1,IDEAL_NUMBER_OF_THREADS()));
    OST_FFTW_fftw_plan plan = OST_FFTW_fftw_plan_dft(rank,n, fftw_in, fftw_out,dir,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    OST_FFTW_fftw_execute(plan);
    OST_FFTW_fftw_destroy_plan(plan);
  }else{
    // use in state
    OST_FFTW_fftw_complex* fftw_in = const_cast<OST_FFTW_fftw_complex*>(reinterpret_cast<const OST_FFTW_fftw_complex*>(in_state.Data().GetData()));
    // out of place transform
    OST_FFTW_fftw_plan_with_nthreads(std::max<int>(1,IDEAL_NUMBER_OF_THREADS()));
    OST_FFTW_fftw_plan plan = OST_FFTW_fftw_plan_dft(rank,n, fftw_in, fftw_out,dir,FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    OST_FFTW_fftw_execute(plan);
    OST_FFTW_fftw_destroy_plan(plan);
  }

  Real fac = 1.0/static_cast<Real>(size.GetVolume());
  for(Complex* ptr = out_state->Data().GetData(); ptr<out_state->Data().GetEnd(); ++ptr) {
    (*ptr)*=fac;
  }

  return out_state;
}

// default impl
template <typename T, class D>
ImageStateBasePtr FFTFnc::VisitState(const ImageStateImpl<T,D>& s) const
{
  std::ostringstream ostr;
  ostr << "FFT  not supported for <" << Val2String<T>() << "," << Domain2String<D>() << ">";
  throw FFTException(ostr.str());
}

// force explicit instantiation to make intel compiler happy
template ImageStateBasePtr FFTFnc::VisitState<Real,FrequencyDomain>(const ImageStateImpl<Real,FrequencyDomain>& ) const;
template ImageStateBasePtr FFTFnc::VisitState<unsigned short,SpatialDomain>(const ImageStateImpl<unsigned short,SpatialDomain>& ) const;

}  //alg ns

template class TEMPLATE_DEF_EXPORT image_state::ImageStateConstModOPAlgorithm<alg::FFTFnc>;

}} // img ns

