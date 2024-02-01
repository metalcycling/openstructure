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
  authors: Valerio Mariani, Andreas Schenk
*/

#include <cmath>

#include <ost/img/alg/fft.hh>
#include <ost/img/data.hh>
#include "fourier_filters.hh"

namespace ost { namespace img { namespace alg { namespace filter_detail {

Real calcF(const Point& point, const Vec3& sampl)
{
  Real F = sqrt((sampl[0]*static_cast<Real>(point[0])*sampl[0]*static_cast<Real>(point[0]))+
                  (sampl[1]*static_cast<Real>(point[1])*sampl[1]*static_cast<Real>(point[1]))+
                  (sampl[2]*static_cast<Real>(point[2])*sampl[2]*static_cast<Real>(point[2])));
  return F;
}



// LowPass

class LPBase 
{

public:

  LPBase (Real limit=0.0):
    limit_(limit)
  {}

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is) const 
  {
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      Vec3 d=is.GetSampling().GetPixelSampling();
      Real F=calcF(Point(it),d);
      Real value = (F <= limit_ ? 1.0 : 0.0);
      is.Value(it)*=Val2Val<Real,T>(value);
    }
  }

  static String GetAlgorithmName() {return "LP"; }

private:
  Real limit_;
};

typedef ImageStateConstModIPAlgorithm<LPBase> LP;




// HighPass

class HPBase 
{

public:

  HPBase (Real limit=0.0):
    limit_(limit)
  {}

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is) const 
  {
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      Vec3 d=is.GetSampling().GetPixelSampling();
      Real F=calcF(Point(it),d);
      Real value= (F >= limit_ ? 1.0 : 0.0);
      is.Value(it)*=Val2Val<Real,T>(value);
    }
  }

  static String GetAlgorithmName() {return "HP"; }

private:
  Real limit_;
};

typedef ImageStateConstModIPAlgorithm<HPBase> HP;



// Gaussian LowPass

class GaussianLPBase 
{

public:

  GaussianLPBase (Real limit=0.0):
    limit_(limit)
  {}

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is) const 
  {
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      Vec3 d=is.GetSampling().GetPixelSampling();
      Real F=calcF(Point(it),d);
      Real value=exp(-F*F/(2.0*limit_*limit_));
      is.Value(it)*=Val2Val<Real,T>(value);
    }
  }

  static String GetAlgorithmName() {return "GaussianLP"; }

private:
  Real limit_;
};

typedef ImageStateConstModIPAlgorithm<GaussianLPBase> GaussianLP;



// Gaussian HighPass

class GaussianHPBase 
{

public:

  GaussianHPBase (Real limit=0.0):
    limit_(limit)
  {}

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is) const 
  {
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      Vec3 d=is.GetSampling().GetPixelSampling();
      Real F=calcF(Point(it),d);
      Real value=1.0-exp(-F*F/(2.0*limit_*limit_));
      is.Value(it)*=Val2Val<Real,T>(value);
    }
  }

  static String GetAlgorithmName() {return "GaussianHP"; }

private:
  Real limit_;
};

typedef ImageStateConstModIPAlgorithm<GaussianHPBase> GaussianHP;




// Fermi LowPass

class FermiLPBase 
{

public:

  FermiLPBase (Real limit=0.0, Real t=0.0):
    limit_(limit),
    t_(t)
  {}

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is) const 
  {
    Vec3 d=is.GetSampling().GetPixelSampling();
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      Vec3 v=Point(it).ToVec3();
      Real l2=Length2(v);
      if(l2==0.0){
        is.Value(it)*=(t_>=0?1.0:0.0);
      }else{
        Real F=calcF(Point(it),d);
        Real teff=t_*t_/l2*F*(t_>=0?1.0:-1.0);
        Real value=1.0/(1.0+exp((F-limit_)/(teff)));
        is.Value(it)*=Val2Val<Real,T>(value);
      }
    }
  }

  static String GetAlgorithmName() {return "FermiLP"; }

private:
  Real limit_;
  Real t_;
};

typedef ImageStateConstModIPAlgorithm<FermiLPBase> FermiLP;




// Butterworth Low pass

class ButterworthLPBase 
{

public:

  ButterworthLPBase (Real limit=0.0, Real stop=0.0, Real eps=0.882, Real a=10.624):
    limit_(limit),
    stop_(stop),
    eps_(eps),
    a_(a)
  {}

  template <typename T, class D>
  void VisitState(ImageStateImpl<T,D>& is) const 
  {
    for(ExtentIterator it(is.GetExtent());!it.AtEnd();++it) {
      Real order=(2.0*log10(eps_/sqrt(a_*a_-1)))/Real(log10(limit_/stop_));
      Real rad=limit_/pow(eps_,Real(2.0)/order);
      Vec3 d=is.GetSampling().GetPixelSampling();
      Real F=calcF(Point(it),d);
      Real value=(1.0/(sqrt(1+pow(F/rad,Real(2.0)*order))));
      is.Value(it)*=Val2Val<Real,T>(value);
    }
  }

  static String GetAlgorithmName() {return "ButterworthLP"; }

private:
  Real limit_;
  Real stop_;
  Real eps_;
  Real a_;      
};

typedef ImageStateConstModIPAlgorithm<ButterworthLPBase> ButterworthLP;


}}}} // ns

#ifndef NO_EXPL_INST
namespace ost { namespace img { namespace image_state {
template class ImageStateConstModIPAlgorithm<ost::img::alg::filter_detail::LPBase>;
template class ImageStateConstModIPAlgorithm<ost::img::alg::filter_detail::HPBase>;
template class ImageStateConstModIPAlgorithm<ost::img::alg::filter_detail::GaussianLPBase>;
template class ImageStateConstModIPAlgorithm<ost::img::alg::filter_detail::GaussianHPBase>;
template class ImageStateConstModIPAlgorithm<ost::img::alg::filter_detail::FermiLPBase>;
template class ImageStateConstModIPAlgorithm<ost::img::alg::filter_detail::ButterworthLPBase>;
}}} // ns
#endif



namespace ost { namespace img { namespace alg {

void LowPassFilter::Visit(ImageHandle& image) const 
{
  if(image.IsSpatial()) {
     filter_detail::LP lp_filter(limit_); 
     image.ApplyIP(FFT());
     image.ApplyIP(lp_filter); 
     image.ApplyIP(FFT());
   } else {
     filter_detail::LP lp_filter(limit_); 
     image.ApplyIP(lp_filter); 
   }
 }                



void HighPassFilter::Visit(ImageHandle& image) const 
{
  if(image.IsSpatial()) {
    filter_detail::HP hp_filter(limit_); 
    image.ApplyIP(FFT());
    image.ApplyIP(hp_filter); 
    image.ApplyIP(FFT());
  } else {
    filter_detail::HP hp_filter(limit_); 
    image.ApplyIP(hp_filter); 
  }
}                



void GaussianLowPassFilter::Visit(ImageHandle& image) const 
{
  if(image.IsSpatial()) {
    filter_detail::GaussianLP gauss_lp_filter(limit_); 
    image.ApplyIP(FFT());
    image.ApplyIP(gauss_lp_filter); 
    image.ApplyIP(FFT());
  } else {
    filter_detail::GaussianLP gauss_lp_filter(limit_); 
    image.ApplyIP(gauss_lp_filter); 
  }
}                




void GaussianHighPassFilter::Visit(ImageHandle& image) const 
{
  if(image.IsSpatial()) {
    filter_detail::GaussianHP gauss_hp_filter(limit_); 
    image.ApplyIP(FFT());
    image.ApplyIP(gauss_hp_filter); 
    image.ApplyIP(FFT());
  } else {
    filter_detail::GaussianHP gauss_hp_filter(limit_); 
    image.ApplyIP(gauss_hp_filter); 
  }
}                



void FermiLowPassFilter::Visit(ImageHandle& image) const 
{
  if(image.IsSpatial()) {
    filter_detail::FermiLP fermi_lp_filter(limit_,t_); 
    image.ApplyIP(FFT());
    image.ApplyIP(fermi_lp_filter); 
    image.ApplyIP(FFT());
  } else {
    filter_detail::FermiLP fermi_lp_filter(limit_,t_); 
    image.ApplyIP(fermi_lp_filter); 
  }
}                



void FermiHighPassFilter::Visit(ImageHandle& image) const 
{
  if(image.IsSpatial()) {
    filter_detail::FermiLP fermi_lp_filter(limit_,t_); 
    image.ApplyIP(FFT());
    image.ApplyIP(fermi_lp_filter); 
    image.ApplyIP(FFT());
  } else {
    filter_detail::FermiLP fermi_lp_filter(limit_,t_); 
    image.ApplyIP(fermi_lp_filter); 
  }
}                



void ButterworthLowPassFilter::Visit(ImageHandle& image) const 
{
   if(image.IsSpatial()) {
     filter_detail::ButterworthLP butt_lp_filter(limit_,stop_,eps_,a_); 
     image.ApplyIP(FFT());
     image.ApplyIP(butt_lp_filter); 
     image.ApplyIP(FFT());
   } else {
     filter_detail::ButterworthLP butt_lp_filter(limit_,stop_,eps_,a_); 
     image.ApplyIP(butt_lp_filter); 
   }
} 
               
        

}}} // ns

