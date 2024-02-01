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
  generic raster image

  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include <stdlib.h>
#ifdef _WIN32
#include <malloc.h>
#endif
#include <cmath>
#include <iostream>
#include <cassert>

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/preprocessor/repetition.hpp>
#include "raster_image.hh"
#include <ost/base.hh>
#include <ost/img/data.hh>
#include <ost/img/point.hh>
#include <ost/img/normalizer_impl.hh>
#include <ost/img/algorithm.hh>
#include <ost/img/image_state.hh>

//#define USE_RASTER_IMAGE_STATE_ALG

#define MAXBINNING 4

namespace ost { namespace img {

namespace {

/*
  adapted from http://www.cs.rit.edu/~ncs/color/t_convert.html
*/
static void HSVtoRGB(float hsv[3],float rgb[3])
{
  if( hsv[1] == 0 ) {
    // achromatic (grey)
    rgb[0] = rgb[1] = rgb[2] = hsv[2];
    return;
  }

  hsv[0] /= 60.0;      // sector 0 to 5
  int i = (int)floor( hsv[0] );
  float f = hsv[0] - (float)i;  // factorial part of h
  float p = hsv[2] * ( 1.0 - hsv[1] );
  float q = hsv[2] * ( 1.0 - hsv[1] * f );
  float t = hsv[2] * ( 1.0 - hsv[1] * ( 1.0 - f ) );

  switch( i ) {
  case 0:
    rgb[0] = hsv[2];
    rgb[1] = t;
    rgb[2] = p;
    break;
  case 1:
    rgb[0] = q;
    rgb[1] = hsv[2];
    rgb[2] = p;
    break;
  case 2:
    rgb[0] = p;
    rgb[1] = hsv[2];
    rgb[2] = t;
    break;
  case 3:
    rgb[0] = p;
    rgb[1] = q;
    rgb[2] = hsv[2];
    break;
  case 4:
    rgb[0] = t;
    rgb[1] = p;
    rgb[2] = hsv[2];
    break;
  default:    // case 5:
    rgb[0] = hsv[2];
    rgb[1] = p;
    rgb[2] = q;
    break;
  }
}
namespace detail{
using namespace boost::mpl;


template<class T>
struct bin{
};

template<>
struct bin<Complex>
{
  unsigned int x_;
  const Point& p_;
  Complex& v_;
  const Data& d_;
  bin(unsigned int x,const Point& p,Complex& v,const Data& d):x_(x),p_(p),v_(v),d_(d){}
  template< typename U > void operator()(U y)
  {
    v_ += d_.GetComplex(p_+Point(x_,y));
  }
};

template<>
struct bin<Real>
{
  unsigned int x_;
  const Point& p_;
  Real& v_;
  const Data& d_;
  bin(unsigned int x,const Point& p,Real& v,const Data& d):x_(x),p_(p),v_(v),d_(d){}
  template< typename U > void operator()(U y)
  {
    v_ += d_.GetReal(p_+Point(x_,y));
  }
};

template<class T,unsigned int numbin>
struct bin_inner_loop
{
  const Point& p_;
  T& v_;
  const Data& d_;
  bin_inner_loop(const Point& p,T& v,const Data& d):p_(p),v_(v),d_(d){}
  template< typename U > void operator()(U x)
  {
    for_each< range_c<int,0,numbin> >( bin<T>(x,p_,v_,d_) );
  }
};

}//ns

template <bool bin, bool intpol, bool phasecolor, bool signcolor,unsigned int b>
RasterImage::Pixel get_value_t(const Data& d, const Point& p, const geom::Vec3& v,
             const NormalizerPtr& n)
{
  if(phasecolor) {
    float hsv[3];
    float rgb[3];
    
    Complex cv(0);
    if(bin) {
      boost::mpl::for_each< boost::mpl::range_c<int,0,b> >( detail::bin_inner_loop<Complex,b>(p,cv,d) );
      cv=n->Convert(cv/static_cast<Real>(b*b));
    } else {
      cv = intpol? d.GetIntpolComplex(v) : d.GetComplex(p);
    }

    hsv[0] = fmod(std::arg(cv)*180.0/M_PI+180.0,360.0); // phase is hue
    hsv[1] = 1.0; // saturation is always 1
    hsv[2] = n->Convert(std::abs(cv))/255.0; // brightness is amplitude
    if(hsv[2]>1.0) hsv[2]=1.0;
    
    HSVtoRGB(hsv,rgb);
    return RasterImage::Pixel(0xff&(int)(rgb[0]*255.0),
                              0xff&(int)(rgb[1]*255.0),
                              0xff&(int)(rgb[2]*255.0));
  } else {
    bool sg = false;
    Real rv=0;
    if(bin) {
      boost::mpl::for_each< boost::mpl::range_c<int,0,b> >( detail::bin_inner_loop<Real,b>(p,rv,d) );
      sg = rv<0.0;
      if(signcolor) rv = std::abs(rv);
      rv=n->Convert(rv/static_cast<Real>(b*b));
    } else {
      Real vv = intpol ? d.GetIntpolReal(v) : d.GetReal(p);
      sg = vv<0.0;
      if(signcolor) vv = std::abs(vv);
      rv = n->Convert(vv);
    }
    uchar v = (unsigned char)::round(rv);
    if(signcolor) {
      float hsv[3];
      hsv[0] = sg ? 0.0 : 120.0;
      hsv[1] = 1.0;
      hsv[2] = std::abs(rv/255.0);
      if(hsv[2]>1.0) hsv[2]=1.0;
      float rgb[3];
      HSVtoRGB(hsv,rgb);
      return RasterImage::Pixel(0xff&(int)(rgb[0]*255.0),
				0xff&(int)(rgb[1]*255.0),
				0xff&(int)(rgb[2]*255.0));
    } else {
      return RasterImage::Pixel(v,v,v);
    }
  }
  return RasterImage::Pixel(0,0,0); // should never get here
}

RasterImage::Pixel get_value(const Data& d, const Point& p,const NormalizerPtr& n) {
  return get_value_t<false,false,false,false,1>(d,p,geom::Vec3(),n);
}

template<unsigned int b>
RasterImage::Pixel get_binned_value(const Data& d, const Point& p,const NormalizerPtr& n) {
  return get_value_t<true,false,false,false,b>(d,p,geom::Vec3(),n);
}

RasterImage::Pixel get_signcolor_value(const Data& d, const Point& p,const NormalizerPtr& n) {
  return get_value_t<false,false,false,true,1>(d,p,geom::Vec3(),n);
}

template<unsigned int b>
RasterImage::Pixel get_binned_signcolor_value(const Data& d, const Point& p,const NormalizerPtr& n) {
  return get_value_t<true,false,false,true,b>(d,p,geom::Vec3(),n);
}

RasterImage::Pixel get_phasecolor_value(const Data& d, const Point& p, const NormalizerPtr& n)
{
  return get_value_t<false,false,true,false,1>(d,p,geom::Vec3(),n);
}

template<unsigned int b>
RasterImage::Pixel get_binned_phasecolor_value(const Data& d, const Point& p, const NormalizerPtr& n) 
{
  return get_value_t<true,false,true,false,b>(d,p,geom::Vec3(),n);
}

RasterImage::Pixel get_interpolated_value(const Data& d, const Vec3& v,const NormalizerPtr& n) 
{
  return get_value_t<false,true,false,false,1>(d,Point(),v,n);
}

RasterImage::Pixel get_interpolated_phasecolor_value(const Data& d, const Vec3& v, const NormalizerPtr& n) 
{
  return get_value_t<false,true,true,false,1>(d,Point(),v,n);
}

#ifdef USE_RASTER_IMAGE_STATE_ALG

template <typename T, class D, bool bin>
RasterImage::Pixel isi_to_pixel(const ImageStateImpl<T,D>& isi,
                                const Point& p,
                                const NormalizerPtr& n,
                                unsigned int b)
{
  Real rv=0;
  if(bin) {
    for(unsigned int i=0;i<b;++i){
      for(unsigned int j=0;j<b;++j){
        rv += n->Convert(Val2Val<T,Real>(isi.Value(p+Point(i,j))));
      }
    } 
    rv/=static_cast<Real>(b*b);  
  } else {
    rv = n->Convert(Val2Val<T,Real>(isi.Value(p)));
  }
  uchar v = (unsigned char)::round(rv);
  return RasterImage::Pixel(v,v,v);
}


class RasterImageFillerBase
{
public:
  RasterImageFillerBase() {}
  RasterImageFillerBase(RasterImage* ri,
                        unsigned int xs_start, unsigned int xs_end, unsigned int xd_start,
                        unsigned int ys_start, unsigned int ys_end, unsigned int yd_start,
                        int z,
                        unsigned int dstep, unsigned int binmax,
                        const NormalizerPtr& norm):
    ri_(ri),
    xs_start_(xs_start),xs_end_(xs_end),xd_start_(xd_start),
    ys_start_(ys_start),ys_end_(ys_end),yd_start_(yd_start),
    z_(z),
    dstep_(dstep),binmax_(binmax),
    norm_(norm)
  {}
    

  template <typename T, class D>
  void VisitState(const ImageStateImpl<T,D>& isi)
  {
    int yd = yd_start_;
    for(unsigned int ys=ys_start_;ys<ys_end_;++ys,yd+=dstep_) {
      int xd = xd_start_;
      for(unsigned int xs=xs_start_;xs<xs_end_;++xs,xd+=dstep_) {
        ri_->SetPixel(xs,ys, isi_to_pixel<T,D,false>(isi,Point(xd,yd,z_),norm_,binmax_));
      }
    }

  }


  static String GetAlgorithmName() {return "RasterImageFiller";}

private:
  RasterImage* ri_;
  unsigned int xs_start_, xs_end_, xd_start_;
  unsigned int ys_start_, ys_end_, yd_start_;
  int z_;
  unsigned int dstep_, binmax_;
  NormalizerPtr norm_;
};

typedef ImageStateNonModAlgorithm<RasterImageFillerBase> RasterImageFiller;

#endif

} // anon ns

RasterImage::RasterImage(unsigned int width, unsigned int height):
  width_(width), height_(height)
{
  data_ = new Pixel[width_*height_];

  external_data_=false;
}

// Initialization with size and already allocated memory area
RasterImage::RasterImage(unsigned int width, unsigned int height, uchar *dptr):
  width_(width), height_(height)
{
  data_ = (Pixel *)dptr;
  external_data_=true;
}

RasterImage::~RasterImage()
{
  if(!external_data_)  delete[] data_;
}

void RasterImage::SetPixel(unsigned int x, unsigned int y, const Pixel& p)
{
  assert(x<width_ &&  y<height_);
  data_[wh2i(x,y)]=p;
}

RasterImage::Pixel RasterImage::GetPixel(unsigned int x, unsigned int y) const
{
  if(x>=width_ || y>=height_) 
    return Pixel();
  else
    return data_[wh2i(x,y)];
}

void RasterImage::Fill(const Data& d, int logscale, const Vec3& offset, int z, const NormalizerPtr& norm,
                       Mode mode, bool fast_low_mag_flag, bool fast_high_mag_flag)
{
  Fill(d,logscale,offset,z,norm,mode,fast_low_mag_flag,fast_high_mag_flag,0,0,width_,height_);
}

#define BOOST_ADD_CASE(z, n, func) \
  case n+1: \
    get_value_fptr=func<n+1>; \
    break;

void RasterImage::Fill(const Data& d, int logscale, const Vec3& offset, int z, const NormalizerPtr& norm,
                       Mode mode, bool fast_low_mag_flag, bool fast_high_mag_flag,int x1,int y1,int x2,int y2)
{
  /*
  screen coordinates: xs,ys,xs_start,xs_end,ys_start,ys_end
  data coordinates:   xd,yd,xd_start,yd_start
  */
   
  Real zoom=ldexp(1.0,logscale);
  int dstep = static_cast<int>(1.0/zoom);
  int binmax=std::min(dstep,MAXBINNING);
  Extent e=d.GetExtent();
  RasterImage::Pixel (*get_value_fptr)(const Data&,const Point&, 
                                       const NormalizerPtr&)=0;

  // adjust x,y start,end to avoid drawing outside of image (speed tweak) + sanity check
  unsigned int xs_start = std::min<unsigned int>(width_ ,static_cast<unsigned int>(std::max<Real>(0.0,std::max<Real>(x1,(e.GetStart()[0]-offset[0])*zoom))));
  unsigned int xs_end   = std::min<unsigned int>(width_ ,static_cast<unsigned int>(std::max<Real>(0.0,std::min<Real>(x2,(e.GetEnd()[0]+1 - offset[0])*zoom))));
  unsigned int ys_start = std::min<unsigned int>(height_,static_cast<unsigned int>(std::max<Real>(0.0,std::max<Real>(y1,(e.GetStart()[1]-offset[1])*zoom))));
  unsigned int ys_end   = std::min<unsigned int>(height_,static_cast<unsigned int>(std::max<Real>(0.0,std::min<Real>(y2,(e.GetEnd()[1]+1 - offset[1])*zoom))));
  if(xs_start>=xs_end || ys_start>=ys_end) return;
  
  Real xd_start_fp=xs_start/zoom+offset[0];
  Real yd_start_fp=ys_start/zoom+offset[1];
  int xd_start=static_cast<int>(floor(xd_start_fp));
  int yd_start=static_cast<int>(floor(yd_start_fp));

  if(mode==PHASECOLOR){
    if(logscale<0 && !fast_low_mag_flag){
      switch(binmax){
        BOOST_PP_REPEAT(MAXBINNING, BOOST_ADD_CASE, get_binned_phasecolor_value)
      }
      //get_value_fptr=get_binned_phasecolor_value;
    }else{
      get_value_fptr=get_phasecolor_value;
    }
  }else{
    if(logscale<0 && !fast_low_mag_flag){
      if(mode==SIGNCOLOR){
        switch(binmax){
          BOOST_PP_REPEAT(MAXBINNING, BOOST_ADD_CASE, get_binned_signcolor_value)
        }
      }else{
        switch(binmax){
          BOOST_PP_REPEAT(MAXBINNING, BOOST_ADD_CASE, get_binned_value)
        }
      }
    }else{
      get_value_fptr= (mode==SIGNCOLOR) ? get_signcolor_value : get_value;
    }
  }

  if(logscale>0) {  // zoom 2:1 or higher
    if(!fast_high_mag_flag){
      RasterImage::Pixel (*get_interpolated_value_fptr)(const Data&,const Vec3&,const NormalizerPtr&);
      Real dstep_fp = 1.0/zoom;
      if(mode==PHASECOLOR){
        get_interpolated_value_fptr=get_interpolated_phasecolor_value;
      }else{
        get_interpolated_value_fptr=get_interpolated_value;
      }
      xd_start_fp+=dstep_fp*0.5-0.5;
      yd_start_fp+=dstep_fp*0.5-0.5;
      Real yd=yd_start_fp;
      for(unsigned int ys=ys_start;ys<ys_end;++ys,yd+=dstep_fp) {
        Real xd = xd_start_fp;
        for(unsigned int xs=xs_start;xs<xs_end;++xs,xd+=dstep_fp) {
          SetPixel(xs,ys, get_interpolated_value_fptr(d,Vec3(xd,yd,z),norm));
        }
      }
    }else{
      unsigned int int_zoom=static_cast<unsigned int>(zoom);
      unsigned int xd_subpixel_base=static_cast<unsigned int>((xd_start_fp-xd_start)*zoom);
      unsigned int yd_subpixel=static_cast<unsigned int>((yd_start_fp-yd_start)*zoom);
      int yd = yd_start;
      for(unsigned int ys=ys_start;ys<ys_end;++yd){
        for(;ys<ys_end && yd_subpixel<int_zoom;++ys,++yd_subpixel){
          unsigned int xd_subpixel=xd_subpixel_base;
          int xd = xd_start;
          for(unsigned int xs=xs_start;xs<xs_end;++xd){
            Pixel pval=get_value_fptr(d,Point(xd,yd,z),norm);
            for(;xs<xs_end && xd_subpixel<int_zoom;++xs,++xd_subpixel){
              SetPixel(xs,ys, pval);
            }
            xd_subpixel=0;
          }
        }
        yd_subpixel=0;
      }
    }
  } else { // zoom 1:1 or lower
#ifdef USE_RASTER_IMAGE_STATE_ALG
    RasterImageFiller rif(this,
                          xs_start,xs_end,xd_start,
                          ys_start,ys_end,yd_start,
                          z,
                          dstep, binmax,
                          norm);
    d.Apply(rif);
#else
    int yd = yd_start;
    for(unsigned int ys=ys_start;ys<ys_end;++ys,yd+=dstep) {
      int xd = xd_start;
      for(unsigned int xs=xs_start;xs<xs_end;++xs,xd+=dstep) {
        SetPixel(xs,ys, get_value_fptr(d,Point(xd,yd,z),norm));
      }
    }
#endif
  }
}


}} // ns
