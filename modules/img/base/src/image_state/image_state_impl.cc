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
  ImageStateImpl definition

  In order to separate implementation from declaration,
  the templated ImageStateImpl is defined here, and called
  from image_state_inst.hh

  Authors: Ansgar Philippsen, Andreas Schenk
*/


#ifndef IMG_IMAGE_STATE_INST_H
#error "image_state_impl.cc must be included from image_state_inst.hh"
#endif

#include <cassert>

#include <ost/base.hh>
#include "image_state_impl.hh"
#include "image_state_spatial_domain.hh"
#include "image_state_visitor.hh"

namespace ost { namespace img { namespace image_state {

template <typename T, class D>
ImageStateImpl<T,D>::ImageStateImpl():
  domain_(Extent()),
  data_(domain_.GetExtent().GetSize())
{}

template <typename T, class D>
ImageStateImpl<T,D>::ImageStateImpl(const Extent& e, const PixelSampling& s):
  domain_(e),
  data_(domain_.GetExtent().GetSize(),domain_.GetPhysicalExtent().GetSize()),
  sampling_(s),
  absolute_origin_(Vec3(0.0,0.0,0.0))
{
  sampling_.SetDomain(domain_.GetDomain());
}

template <typename T, class D>
ImageStateImpl<T,D>::ImageStateImpl(const Extent& e, const PixelSampling& s, const Vec3& c):
  domain_(e),
  data_(domain_.GetExtent().GetSize(),domain_.GetPhysicalExtent().GetSize()),
  sampling_(s),
  absolute_origin_(c)
{
  sampling_.SetDomain(domain_.GetDomain());
}

template <typename T, class D>
ImageStateImpl<T,D>::ImageStateImpl(const ImageStateImpl<T,D>& s):
  domain_(s.domain_),
  data_(s.data_),
  sampling_(s.sampling_),
  absolute_origin_(s.absolute_origin_)
{
  sampling_.SetDomain(domain_.GetDomain());
}


template <typename T, class D>
ImageStateImpl<T,D>& ImageStateImpl<T,D>::operator=(const ImageStateImpl<T,D>& s)
{
  if(this!=&s) {
    domain_=s.domain_;
    data_=s.data_; // copy
    sampling_=s.sampling_;
    absolute_origin_=s.absolute_origin_;
    sampling_.SetDomain(domain_.GetDomain());
  }
  return *this;
}


template <typename T, class D>
ImageStateImpl<T,D>::~ImageStateImpl()
{
  // cleanup is automatically done by the members
}
  
template <typename T, class D>
boost::shared_ptr<ImageStateImpl<T,D> > ImageStateImpl<T,D>::CloneState(bool cc) const
{
  if(cc) {
    return SharedPtrType(new ImageStateImpl<T,D>(*this));
  } else {
    return SharedPtrType(new ImageStateImpl<T,D>(domain_.GetLogicalExtent(),sampling_,this->absolute_origin_));
  }
}

template <typename T, class D>
ImageStateBasePtr ImageStateImpl<T,D>::Clone(bool cc) const
{
  return this->CloneState(cc);
}


template <typename T, class D>
size_t ImageStateImpl<T,D>::MemSize() const
{
  return data_.MemSize();
}

template <typename T, class D>
DataType ImageStateImpl<T,D>::GetType() const
{
  return data_.GetDataType();
}

template <typename T, class D>
DataDomain ImageStateImpl<T,D>::GetDomain() const 
{
  return domain_.GetDomain();
}

template <typename T, class D>
void ImageStateImpl<T,D>::SetSpatialOrigin(const Point& o) 
{
  domain_.SetSpatialOrigin(o);
}

template <typename T, class D>
void ImageStateImpl<T,D>::CenterSpatialOrigin() 
{
  SetSpatialOrigin(Extent(GetExtent().GetSize(),Point(0,0,0)).GetStart());
}

template <typename T, class D>
Point ImageStateImpl<T,D>::GetSpatialOrigin() const 
{
  return domain_.GetSpatialOrigin();
}

template <typename T, class D>
Extent ImageStateImpl<T,D>::GetExtent() const 
{
  return domain_.GetExtent();
}

template <typename T, class D>
Extent ImageStateImpl<T,D>::GetLogicalExtent() const 
{
  return domain_.GetLogicalExtent();
}

template <typename T, class D>
Real ImageStateImpl<T,D>::GetReal(const Point &p) const 
{
  return domain_.GetReal(p,data_);
}

template <typename T, class D>
void ImageStateImpl<T,D>::SetReal(const Point &p, const Real& r) 
{
  domain_.SetReal(p,r,data_);
}
 
template <typename T, class D>
Complex ImageStateImpl<T,D>::GetComplex(const Point &p) const 
{
  return domain_.GetComplex(p,data_);
}

template <typename T, class D>
void ImageStateImpl<T,D>::SetComplex(const Point &p, const Complex& c)
{
  domain_.SetComplex(p,c,data_);
}

template <typename T, class D>
Real ImageStateImpl<T,D>::GetIntpolReal(const Vec3 &v) const
{
  return Val2Val<T,Real>(CalcIntpolValue(v));
}
template <typename T, class D>
Real ImageStateImpl<T,D>::GetIntpolReal(const Vec2 &v) const
{
  return Val2Val<T,Real>(CalcIntpolValue(v));
}
template <typename T, class D>
Real ImageStateImpl<T,D>::GetIntpolReal(const Real &d) const
{
  return Val2Val<T,Real>(CalcIntpolValue(d));
}

template <typename T, class D>
Complex ImageStateImpl<T,D>::GetIntpolComplex(const Vec3 &v) const
{
  return Val2Val<T,Complex>(CalcIntpolValue(v));
}
template <typename T, class D>
Complex ImageStateImpl<T,D>::GetIntpolComplex(const Vec2 &v) const
{
  return Val2Val<T,Complex>(CalcIntpolValue(v));
}
template <typename T, class D>
Complex ImageStateImpl<T,D>::GetIntpolComplex(const Real &d) const
{
  return Val2Val<T,Complex>(CalcIntpolValue(d));
}


template <typename T, class D>
PixelSampling& ImageStateImpl<T,D>::GetSampling() 
{
  return sampling_;
}

template <typename T, class D>
const PixelSampling& ImageStateImpl<T,D>::GetSampling() const
{
  return sampling_;
}

template <typename T, class D>
void ImageStateImpl<T,D>::SetSampling(const PixelSampling& s)
{
  sampling_=s;
  sampling_.SetDomain(domain_.GetDomain());
}

template <typename T, class D>
Vec3 ImageStateImpl<T,D>::GetAbsoluteOrigin() const
{
  return absolute_origin_;
}

template <typename T, class D>
void ImageStateImpl<T,D>::SetAbsoluteOrigin(const Vec3& center)
{
  absolute_origin_ = center;
}

template <typename T, class D>
Vec3 ImageStateImpl<T,D>::IndexToCoord(const Point &i) const
{
  Vec3 offset(0.0,0.0,0.0);
  if (domain_.GetDomain() == SPATIAL)
  {
    offset = absolute_origin_;
  }
  return sampling_.Point2Coord(i)+offset;
}

template <typename T, class D>
Vec3 ImageStateImpl<T,D>::CoordToIndex(const Vec3& c) const
{
  Vec3 index;
  if (domain_.GetDomain() == SPATIAL)
  {
    index = sampling_.Coord2Point(c-absolute_origin_);
  } else {
    index = sampling_.Coord2Point(c);
  }
  return index;
}

template <typename T, class D>
Vec3 ImageStateImpl<T,D>::FractionalIndexToCoord(const Vec3 &p) const
{
  Vec3 offset(0.0,0.0,0.0);
  if (domain_.GetDomain() == SPATIAL)
  {
    offset = absolute_origin_;
  }
  return sampling_.Vec2Coord(p)+offset;
}


template <typename T, class D>
Size ImageStateImpl<T,D>::GetSize() const 
{
  return GetExtent().GetSize();
}

template <typename T, class D>
void ImageStateImpl<T,D>::AdjustPhaseOrigin(const Point& p)
{
  if(p==Point(0,0,0)) return;
  Size sz=this->GetLogicalExtent().GetSize();
  Real fu = 2.0*M_PI*Real(-p[0])/Real(sz[0]);
  Real fv = 2.0*M_PI*Real(-p[1])/Real(sz[1]);
  Real fw = 2.0*M_PI*Real(-p[2])/Real(sz[2]);

  for(IndexIterator it(this->GetSize()); !it.AtEnd(); ++it) {
    Index indx(it);
    Real ex = fu*static_cast<Real>(indx.u)+fv*static_cast<Real>(indx.v)+fw*static_cast<Real>(indx.w);
    Value(it)=Value(it)*Val2Val<Complex,T>(Complex(cos(ex),sin(ex)));
  }
}



template <typename T, class D>
T ImageStateImpl<T,D>::GetCheckedValue(const Point& p) const
{
  if(domain_.GetExtent().Contains(p)) return this->Value(p);
  else return T();
}

template <typename T, class D>
void ImageStateImpl<T,D>::SetCheckedValue(const Point& p, const T& v)
{
  if(domain_.GetExtent().Contains(p)) this->Value(p)=v;
}

namespace {

/*
  interpolation for 3D case

    /            / (1-r)
  -/------------/-
  /        |   / r
  --------------
  |        |  |         000 = (1-p) * (1-q) * (1-r)
  |   D    |C | (1-q)   100 =   p   * (1-q) * (1-r)
  |        |  |         010 = (1-p) *   q   * (1-r)
 -|--------|--|-        001 = (1-p) * (1-q) *   r
  |   A    |B | q       110 =   p   *   q   * (1-r)
  -------------         101 = ... etc
     p     |(1-p)
*/
static void calc_intpol_weights_3d(const Vec3& v, Real w[8], Point p[8])
{
  Real dix = floor(v[0]);
  Real diy = floor(v[1]);
  Real diz = floor(v[2]);

  Real dx = v[0]-dix;
  Real dy = v[1]-diy;
  Real dz = v[2]-diz;
  Real dxdy = dx*dy;
  Real dxdz = dx*dz;
  Real dydz = dy*dz;

  int ix = static_cast<int>(dix);
  int iy = static_cast<int>(diy);
  int iz = static_cast<int>(diz);

  p[0][0] = ix+0;  p[0][1] = iy+0;  p[0][2] = iz+0;
  p[1][0] = ix+1;  p[1][1] = iy+0;  p[1][2] = iz+0;
  p[2][0] = ix+0;  p[2][1] = iy+1;  p[2][2] = iz+0;
  p[3][0] = ix+1;  p[3][1] = iy+1;  p[3][2] = iz+0;
  p[4][0] = ix+0;  p[4][1] = iy+0;  p[4][2] = iz+1;
  p[5][0] = ix+1;  p[5][1] = iy+0;  p[5][2] = iz+1;
  p[6][0] = ix+0;  p[6][1] = iy+1;  p[6][2] = iz+1;
  p[7][0] = ix+1;  p[7][1] = iy+1;  p[7][2] = iz+1;

  w[7] = dxdy*dz;                            //   dx  *  dy  *  dz
  w[6] = dydz-w[7];                          // (1-dx)*  dy  *  dz
  w[5] = dxdz-w[7];                          //   dx  *(1-dy)*  dz
  w[4] = dz-dxdz-w[6];                       // (1-dx)*(1-dy)*  dz
  w[3] = dxdy-w[7];                          //   dx  *  dy  *(1-dz)
  w[2] = dy-dxdy-w[6];                       // (1-dx)*  dy  *(1-dz)
  w[1] = dx-dxdy-w[5];                       //   dx  *(1-dy)*(1-dz)
  w[0] = 1-w[4]-dx-dy+dxdy;                  // (1-dx)*(1-dy)*(1-dz)
}

template <typename T>
T do_calc_intpol_value_3d(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7, Real w[8]);

template <>
Word do_calc_intpol_value_3d(const Word& v0, const Word& v1, const Word& v2, const Word& v3, const Word& v4, const Word& v5, const Word& v6, const Word& v7, Real w[8])
{
  Word nrvo = static_cast<Word>(w[0]*static_cast<Real>(v0)+
				w[1]*static_cast<Real>(v1)+
				w[2]*static_cast<Real>(v2)+
				w[3]*static_cast<Real>(v3)+
				w[4]*static_cast<Real>(v4)+
				w[5]*static_cast<Real>(v5)+
				w[6]*static_cast<Real>(v6)+
				w[7]*static_cast<Real>(v7));
  return nrvo;
}

template <typename T>
T do_calc_intpol_value_3d(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4, const T& v5, const T& v6, const T& v7, Real w[8])
{
  T nrvo = w[0]*v0+w[1]*v1+w[2]*v2+w[3]*v3+w[4]*v4+w[5]*v5+w[6]*v6+w[7]*v7;
  return nrvo;
}

static void calc_intpol_weights_2d(const Vec2& v, int iz, Real w[4], Point p[4])
{
  Real dix = floor(v[0]);
  Real diy = floor(v[1]);

  Real dx = v[0]-dix;
  Real dy = v[1]-diy;

  int ix = static_cast<int>(dix);
  int iy = static_cast<int>(diy);

  p[0][0] = ix+0;  p[0][1] = iy+0;  p[0][2] = iz;
  p[1][0] = ix+1;  p[1][1] = iy+0;  p[1][2] = iz;
  p[2][0] = ix+0;  p[2][1] = iy+1;  p[2][2] = iz;
  p[3][0] = ix+1;  p[3][1] = iy+1;  p[3][2] = iz;

  w[3] = dx*dy;
  w[2] = dy-w[3];
  w[1] = dx-w[3];
  w[0] = Real(1.0)-dx-dy+w[3];
}


template <typename T>
T do_calc_intpol_value_2d(const T& v0, const T& v1, const T& v2, const T& v3, Real w[4]);

template <>
Word do_calc_intpol_value_2d(const Word& v0, const Word& v1, const Word& v2, const Word& v3, Real w[4])
{
  Word nrvo = static_cast<Word>(w[0]*static_cast<Real>(v0)+
              w[1]*static_cast<Real>(v1)+
              w[2]*static_cast<Real>(v2)+
              w[3]*static_cast<Real>(v3));
  return nrvo;
}

template <typename T>
T do_calc_intpol_value_2d(const T& v0, const T& v1, const T& v2, const T& v3, Real w[4])
{
  T nrvo = w[0]*v0+w[1]*v1+w[2]*v2+w[3]*v3;
  return nrvo;
}

static void calc_intpol_weights_1d(const Real& d,int iy, int iz, Real w[2], Point p[2])
{
  Real dix = floor(d);
  Real dx = d-dix;

  int ix = static_cast<int>(dix);

  p[0][0] = ix+0;  p[0][1] = iy;  p[0][2] = iz;
  p[1][0] = ix+1;  p[1][1] = iy;  p[1][2] = iz;

  w[1] = dx;
  w[0] = Real(1.0)-dx;
}

template <typename T>
T do_calc_intpol_value_1d(const T& v0, const T& v1, Real w[2]);

template <>
Word do_calc_intpol_value_1d(const Word& v0, const Word& v1, Real w[2])
{
  Word nrvo = static_cast<Word>(w[0]*static_cast<Real>(v0)+
              w[1]*static_cast<Real>(v1));
  return nrvo;
}

template <typename T>
T do_calc_intpol_value_1d(const T& v0, const T& v1, Real w[2])
{
  T nrvo = w[0]*v0+w[1]*v1;
  return nrvo;
}

} // anon ns


template <typename T, class D>
T ImageStateImpl<T,D>::CalcIntpolValue(const Vec3& v) const
{
  Real w[8];
  Point p[8];

  calc_intpol_weights_3d(v,w,p);
  if(domain_.GetExtent().Contains(p[0]) && domain_.GetExtent().Contains(p[7])){
    return do_calc_intpol_value_3d(Value(p[0]),Value(p[1]),
                                   Value(p[2]),Value(p[3]),
                                   Value(p[4]),Value(p[5]),
                                   Value(p[6]),Value(p[7]),
                                   w);
  }else{
    return do_calc_intpol_value_3d(GetCheckedValue(p[0]),
                                   GetCheckedValue(p[1]),
                                   GetCheckedValue(p[2]),
                                   GetCheckedValue(p[3]),
                                   GetCheckedValue(p[4]),
                                   GetCheckedValue(p[5]),
                                   GetCheckedValue(p[6]),
                                   GetCheckedValue(p[7]),
                                   w);
  }
}

template <typename T, class D>
T ImageStateImpl<T,D>::CalcIntpolValue(const Vec2& v) const
{
  Real w[4];
  Point p[4];

  calc_intpol_weights_2d(v,domain_.GetExtent().GetStart()[2],w,p);
  if(domain_.GetExtent().Contains(p[0]) && domain_.GetExtent().Contains(p[3])){
    return do_calc_intpol_value_2d(Value(p[0]),Value(p[1]),
                                     Value(p[2]),Value(p[3]),
                                   w);
  }else{
    return do_calc_intpol_value_2d(GetCheckedValue(p[0]),
                                   GetCheckedValue(p[1]),
                                   GetCheckedValue(p[2]),
                                   GetCheckedValue(p[3]),
                                   w);
  }
}

template <typename T, class D>
T ImageStateImpl<T,D>::CalcIntpolValue(const Real& d) const
{
  Real w[2];
  Point p[2];

  calc_intpol_weights_1d(d,domain_.GetExtent().GetStart()[1],domain_.GetExtent().GetStart()[2],w,p);
  return do_calc_intpol_value_1d(GetCheckedValue(p[0]),
                                 GetCheckedValue(p[1]),
                                 w);
}

template <typename T, class D>
T& ImageStateImpl<T,D>::Value(const Index& i)
{
  return data_.Value(i);
}

template <typename T, class D>
const T& ImageStateImpl<T,D>::Value(const Index& i) const
{
  return data_.Value(i);
}

template <typename T, class D>
T& ImageStateImpl<T,D>::Value(size_t i)
{
    return data_.GetData()[i];
}

template <typename T, class D>
const T& ImageStateImpl<T,D>::Value(size_t i) const
{
    return data_.GetData()[i];
}

template <typename T, class D>
ValueHolder<T>& ImageStateImpl<T,D>::Data()
{
  return data_;
}

template <typename T, class D>
const ValueHolder<T>& ImageStateImpl<T,D>::Data() const
{
  return data_;
}


template <typename T, class D>
void ImageStateImpl<T,D>::ApplyIP(ImageStateNonModVisitorBase& v) const
{
  v.Visit_(*this);
}

template <typename T, class D>
void ImageStateImpl<T,D>::Apply(ImageStateNonModVisitorBase& v) const
{
  v.Visit_(*this);
}

template <typename T, class D>
void ImageStateImpl<T,D>::ApplyIP(ImageStateModIPVisitorBase& v) 
{
  v.Visit_(*this);
}

template <typename T, class D>
void ImageStateImpl<T,D>::ApplyIP(const ImageStateConstModIPVisitorBase& v) 
{
  v.Visit_(*this);
}

template <typename T, class D>
ImageStateBasePtr ImageStateImpl<T,D>::Apply(ImageStateModIPVisitorBase& v) const
{
  SharedPtrType s=this->CloneState();
  v.Visit_(*s);
  return s;
}

template <typename T, class D>
ImageStateBasePtr ImageStateImpl<T,D>::Apply(const ImageStateConstModIPVisitorBase& v) const
{
  SharedPtrType s=this->CloneState();
  v.Visit_(*s);
  return s;
}

template <typename T, class D>
ImageStateBasePtr ImageStateImpl<T,D>::Apply(ImageStateModOPVisitorBase& v) const 
{
  return v.Visit_(*this);
}

template <typename T, class D>
ImageStateBasePtr ImageStateImpl<T,D>::Apply(const ImageStateConstModOPVisitorBase& v) const 
{
  return v.Visit_(*this);
}

template <typename T, class D>
ImageStateBasePtr ImageStateImpl<T,D>::Apply(ImageStateMorphVisitorBase& v) 
{
  return v.Visit_(*this);
}

template <typename T, class D>
ImageStateBase& ImageStateImpl<T,D>::operator+=(Real v) 
{
  add_scalar(v);
  return *this;
}

template <typename T, class D>
ImageStateBase& ImageStateImpl<T,D>::operator+=(const Complex& v)
{
  add_scalar(v);
  return *this;
}

template <typename T, class D>
ImageStateBase& ImageStateImpl<T,D>::operator-=(Real v) 
{
  add_scalar(-v);
  return *this;
};

template <typename T, class D>
ImageStateBase& ImageStateImpl<T,D>::operator-=(const Complex& v) 
{
  add_scalar(-v);
  return *this;
}

template <typename T, class D>
ImageStateBase& ImageStateImpl<T,D>::operator*=(Real v) 
{
  mul_scalar(v);
  return *this;
}

template <typename T, class D>
ImageStateBase& ImageStateImpl<T,D>::operator*=(const Complex& v) 
{
  mul_scalar(v);
  return *this;
}

template <typename T, class D>
ImageStateBase& ImageStateImpl<T,D>::operator/=(Real v) 
{
  mul_scalar(1.0/v);
  return *this;
}

template <typename T, class D>
ImageStateBase& ImageStateImpl<T,D>::operator/=(const Complex& v) 
{
  mul_scalar(Real(1.0)/v);
  return *this;
}


template <typename T, class D>
template <typename U>
void ImageStateImpl<T,D>::add_scalar(const U& v) 
{
  T v2=Val2Val<U,T>(v);
  for(T* ptr=data_.GetData(); ptr<data_.GetEnd(); ++ptr)
    (*ptr)+=v2;
}

template <typename T, class D>
template <typename U> void ImageStateImpl<T,D>::mul_scalar(const U& v) 
{
  T v2=Val2Val<U,T>(v);
  for(T* ptr=data_.GetData(); ptr<data_.GetEnd(); ++ptr)
    (*ptr)*=v2;
}

}}} // ns
