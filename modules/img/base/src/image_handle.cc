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
  image handle and const handle

  Author: Ansgar Philippsen, Andreas Schenk
*/

#include <ost/message.hh>
#include <ost/log.hh>
#include <ost/img/module_config.hh>

#include "image_handle.hh"
#include "point.hh"
#include "size.hh"
#include "extent.hh"
#include "data.hh"
#include "data_algorithm.hh"
#include "image_impl.hh"
#include "observable.hh"
#include "paste_impl.hh"
#include "data_observer.hh"

#include "image_state/image_state_factory.hh"
#include "image_state/image_state_base.hh"

namespace ost { namespace img {

// ImageHandle

// construction/deconstruction/assignement

ImageHandle ImageHandle::Create(const Extent& e, DataType type, DataDomain dom)
{
  return ImageHandle(StatePtrPtr((new ImageStateBasePtr(image_state::CreateState(e,type,dom)))));
}

ImageHandle::ImageHandle():
  Data(),
  impl_(),
  obs_()
{
}

ImageHandle::ImageHandle(const ImageHandle& h):
  Data(h),
  impl_(h.impl_),
  obs_(h.obs_)
{}

ImageHandle& ImageHandle::operator=(const ImageHandle& h)
{
  Data::operator=(h);
  impl_=h.impl_;
  obs_=h.obs_;
  Notify();
  return *this;
}

ImageHandle::ImageHandle(const StatePtrPtr& spp):
  Data(),
  impl_(spp),
  obs_(new DataObservable())
{
}

ImageHandle ImageHandle::Copy(bool cc) const
{
  /*
    While this looks complicated, it isnt. First, a clone of the current
    image state is created, which is returned as an ImageStateBasePtr. 
    Then, a new ImageStateBasePtr is created on the heap, using the
    copy constructor,with the cloned ImageStateBasePtr as a parameter.
    This new ptr on the heap is then in turn used to initialize the StatePtrPtr.
  */
  return ImageHandle(StatePtrPtr(new ImageStateBasePtr(ImageStatePtr()->Clone(cc))));
}

ImageHandle ImageHandle::Extract(const Extent& e) const
{
#if 0
  ImageHandle ih = CreateImage(e,this->GetType(),this->GetDomain());
  
  if(this->GetType()==REAL) {
    for(ExtentIterator it(e); !it.AtEnd();++it) {
      ih.SetReal(it,this->GetReal(it));
    }
  } else {
    for(ExtentIterator it(e); !it.AtEnd();++it) {
      ih.SetComplex(it,this->GetComplex(it));
    }
  }
#else
  // hardcoded spatial on purpose
  ImageHandle ih = CreateImage(e,this->GetType(),SPATIAL);
  // carry along pixel sampling as well, but don't convert from reci to spat
  ih.SetPixelSampling(this->GetPixelSampling());

  ih.Paste(*this);
#endif
  return ih;
}

void ImageHandle::Paste(const Data& d)
{
  detail::PasteImpl pst(*this);
  try {
    d.Apply(pst);
  } catch(InvalidExtentException&) {
    // ignore
  }
  Notify();
}

void ImageHandle::Set(const ImageHandle& h)
{
  impl_=h.impl_;
  obs_=h.obs_;
}


void ImageHandle::Swap(ImageHandle& h)
{
  impl_.swap(h.impl_);
  obs_.swap(h.obs_);
  Notify();
}

void ImageHandle::Reset(const Extent &e, DataType type,DataDomain dom)
{
  if(impl_) impl_=StatePtrPtr(); // clears memory
  impl_=StatePtrPtr(new ImageStateBasePtr(image_state::CreateState(e,type,dom)));
  Notify();
}

bool ImageHandle::IsValid() const
{
  return static_cast<bool>(impl_);
}

long ImageHandle::MemSize() const
{
  if(impl_) {
    return ImageStatePtr()->MemSize()+(*obs_).MemSize()+sizeof(*this);
  } else {
    return sizeof(*this);
  }
}

void ImageHandle::Apply(NonModAlgorithm& a) const 
{
  a.Visit(*this);
}

void ImageHandle::ApplyIP(NonModAlgorithm& a) const 
{
  this->Apply(a);
  Notify();
}

void ImageHandle::ApplyIP(ModIPAlgorithm& a) 
{
  a.Visit(*this);
  Notify();
}

ImageHandle ImageHandle::Apply(ModIPAlgorithm& a) const 
{
  ImageHandle nh = this->Copy(true);
  nh.ApplyIP(a);
  return nh;
}

void ImageHandle::ApplyIP(const ConstModIPAlgorithm& a) 
{
  a.Visit(*this);
  Notify();
}

ImageHandle ImageHandle::Apply(const ConstModIPAlgorithm& a) const 
{
  ImageHandle nh = this->Copy(true);
  nh.ApplyIP(a);
  return nh;
}

void ImageHandle::ApplyIP(ModOPAlgorithm& a) 
{
  ImageHandle nh = a.Visit(*this);
  ImageStatePtr().swap(nh.ImageStatePtr());
  Notify();
}

ImageHandle ImageHandle::Apply(ModOPAlgorithm& a) const 
{
  ImageHandle nh = a.Visit(*this);
  return nh;
}

void ImageHandle::ApplyIP(const ConstModOPAlgorithm& a) 
{
  ImageHandle nh = a.Visit(*this);
  ImageStatePtr().swap(nh.ImageStatePtr());
  Notify();
}

ImageHandle ImageHandle::Apply(const ConstModOPAlgorithm& a) const 
{
  ImageHandle nh = a.Visit(*this);
  return nh;
}

// image state visitor interface

void ImageHandle::StateApply(image_state::ImageStateNonModVisitorBase& v) const
{
  ImageStatePtr()->ApplyIP(v);
}

void ImageHandle::StateApplyIP(image_state::ImageStateModIPVisitorBase& v)
{
  ImageStatePtr()->ApplyIP(v);
  Notify();
}

ImageHandle ImageHandle::StateApply(image_state::ImageStateModIPVisitorBase& v) const
{
  ImageHandle nh=this->Copy();
  nh.StateApplyIP(v);
  return nh;
}

void ImageHandle::StateApplyIP(const image_state::ImageStateConstModIPVisitorBase& v)
{
  ImageStatePtr()->ApplyIP(v);
  Notify();
}

ImageHandle ImageHandle::StateApply(const image_state::ImageStateConstModIPVisitorBase& v) const
{
  ImageHandle nh=this->Copy();
  nh.StateApplyIP(v);
  return nh;
}

void ImageHandle::StateApplyIP(image_state::ImageStateModOPVisitorBase& v)
{
  ImageHandle nh=this->StateApply(v);
  ImageStatePtr().swap(nh.ImageStatePtr());
  Notify();
}

ImageHandle ImageHandle::StateApply(image_state::ImageStateModOPVisitorBase& v) const
{
  return ImageHandle(StatePtrPtr(new ImageStateBasePtr(ImageStatePtr()->Apply(v))));
}

void ImageHandle::StateApplyIP(const image_state::ImageStateConstModOPVisitorBase& v)
{
  ImageHandle nh=this->StateApply(v);
  ImageStatePtr().swap(nh.ImageStatePtr());
  Notify();
}

ImageHandle ImageHandle::StateApply(const image_state::ImageStateConstModOPVisitorBase& v) const
{
  return ImageHandle(StatePtrPtr(new ImageStateBasePtr(ImageStatePtr()->Apply(v))));
}

void ImageHandle::StateApplyIP(image_state::ImageStateMorphVisitorBase& v)
{
  ImageStateBasePtr nb(ImageStatePtr()->Apply(v));
  ImageStatePtr().swap(nb);
}

ImageHandle ImageHandle::StateApply(image_state::ImageStateMorphVisitorBase& v) const
{
  ImageStateBasePtr nb(ImageStatePtr()->Apply(v));
  return ImageHandle(StatePtrPtr(new ImageStateBasePtr(nb)));
}

void ImageHandle::Attach(DataObserver *o) const
{
  ImageStatePtr(); // assert valid
  obs_.get()->Attach(o);
}

void ImageHandle::Detach(DataObserver *o) const
{
  ImageStatePtr(); // assert valid
  obs_.get()->Detach(o);
}

void ImageHandle::Notify() const
{
  ImageStatePtr(); // assert valid
  obs_.get()->Notify();
}
void ImageHandle::Notify(const Extent& e) const
{
  ImageStatePtr(); // assert valid
  obs_.get()->Notify(e);
}
void ImageHandle::Notify(const Point& p) const
{
  ImageStatePtr(); // assert valid
  obs_.get()->Notify(p);
}

DataType ImageHandle::GetType() const 
{
  return ImageStatePtr()->GetType();
}

Extent ImageHandle::GetExtent() const 
{
  return ImageStatePtr()->GetLogicalExtent();
}

Point ImageHandle::GetSpatialOrigin() const 
{
  return ImageStatePtr()->GetSpatialOrigin();
}

void ImageHandle::SetSpatialOrigin(const Point& p) 
{
  // TODO: adjust phase origin
  ImageStatePtr()->SetSpatialOrigin(p);
  Notify();
}

void ImageHandle::CenterSpatialOrigin()
{
  ImageStatePtr()->SetSpatialOrigin(Extent(GetExtent().GetSize(),Point(0,0,0)).GetStart());
  Notify();
}

Vec3 ImageHandle::GetAbsoluteOrigin() const
{
  return ImageStatePtr()->GetAbsoluteOrigin();
}

void ImageHandle::SetAbsoluteOrigin(const Vec3& c)
{
  ImageStatePtr()->SetAbsoluteOrigin(c);
  Notify();
}

Vec3 ImageHandle::IndexToCoord(const Point &p) const
{
  return ImageStatePtr()->IndexToCoord(p);
}

Vec3 ImageHandle::CoordToIndex(const Vec3& c) const
{
  return ImageStatePtr()->CoordToIndex(c);
}

Vec3 ImageHandle::FractionalIndexToCoord(const Vec3 &p) const
{
  return ImageStatePtr()->FractionalIndexToCoord(p);
}

DataDomain ImageHandle::GetDomain() const 
{
  return ImageStatePtr()->GetDomain();
}

Real ImageHandle::GetReal(const Point& p) const 
{
  return ImageStatePtr()->GetReal(p);
}

Complex ImageHandle::GetComplex(const Point& p) const 
{
  return ImageStatePtr()->GetComplex(p);
}

void ImageHandle::SetReal(const Point &p, Real r) 
{
  ImageStatePtr()->SetReal(p,r);
  Notify(p);
}

void ImageHandle::SetComplex(const Point &p, const Complex& c) 
{
  ImageStatePtr()->SetComplex(p,c);
  Notify(p);
}
Real ImageHandle::GetIntpolReal(const Vec3 &v) const
{
  return ImageStatePtr()->GetIntpolReal(v);
}
Real ImageHandle::GetIntpolReal(const Vec2 &v) const
{
  return ImageStatePtr()->GetIntpolReal(v);
}
Real ImageHandle::GetIntpolReal(const Real &d) const
{
  return ImageStatePtr()->GetIntpolReal(d);
}
Complex ImageHandle::GetIntpolComplex(const Vec3 &v) const
{
  return ImageStatePtr()->GetIntpolComplex(v);
}
Complex ImageHandle::GetIntpolComplex(const Vec2 &v) const
{
  return ImageStatePtr()->GetIntpolComplex(v);
}
Complex ImageHandle::GetIntpolComplex(const Real &d) const
{
  return ImageStatePtr()->GetIntpolComplex(d);
}

PixelSampling& ImageHandle::Sampling()
{
  Notify();  
  return ImageStatePtr()->GetSampling();
}

const PixelSampling& ImageHandle::Sampling() const
{
  return ImageStatePtr()->GetSampling();
}

ImageStateBasePtr& ImageHandle::ImageStatePtr()
{
  if(impl_) {
    return *impl_;
  } else {
    throw InvalidImageHandle();
  }
}

const ImageStateBasePtr& ImageHandle::ImageStatePtr() const
{
  if(impl_) {
    return *impl_;
  } else {
    throw InvalidImageHandle();
  }
}

// operators

bool ImageHandle::operator==(const ImageHandle& ih) const
{
  return equal(ih);
}

bool ImageHandle::operator!=(const ImageHandle& ih) const
{
  return !equal(ih);
}

bool ImageHandle::equal(const ImageHandle& ih) const
{
  return impl_ == ih.impl_;
}


ImageHandle& ImageHandle::operator+=(Real v)
{
  ImageStatePtr()->operator+=(v);
  Notify();
  return *this;
}

ImageHandle& ImageHandle::operator+=(const Complex& v)
{
  ImageStatePtr()->operator+=(v);
  Notify();
  return *this;
}

ImageHandle& ImageHandle::operator-=(Real v)
{
  ImageStatePtr()->operator-=(v);
  Notify();
  return *this;
}

ImageHandle& ImageHandle::operator-=(const Complex& v)
{
  ImageStatePtr()->operator-=(v);
  Notify();
  return *this;
}

ImageHandle& ImageHandle::operator*=(Real v)
{
  ImageStatePtr()->operator*=(v);
  Notify();
  return *this;
}

ImageHandle& ImageHandle::operator*=(const Complex& v)
{
  ImageStatePtr()->operator*=(v);
  Notify();
  return *this;
}

ImageHandle& ImageHandle::operator/=(Real v)
{
  ImageStatePtr()->operator/=(v);
  Notify();
  return *this;
}

ImageHandle& ImageHandle::operator/=(const Complex& v)
{
  ImageStatePtr()->operator/=(v);
  Notify();
  return *this;
}

ImageHandle& ImageHandle::operator+=(const ConstImageHandle& h)
{
  try {
    ImageStatePtr()->operator+=(*h.ImageStatePtr());
  } catch(InvalidExtentException&) {
    // ignore
  }
  return *this;
}

ImageHandle& ImageHandle::operator-=(const ConstImageHandle& h)
{
  try {
    ImageStatePtr()->operator-=(*h.ImageStatePtr());
  } catch(InvalidExtentException&) {
    // ignore
  }
  return *this;
}

ImageHandle& ImageHandle::operator*=(const ConstImageHandle& h)
{
  try {
    ImageStatePtr()->operator*=(*h.ImageStatePtr());
  } catch(InvalidExtentException&) {
    // ignore
  }
  return *this;
}

ImageHandle& ImageHandle::operator/=(const ConstImageHandle& h)
{
  try {
    ImageStatePtr()->operator/=(*h.ImageStatePtr());
  } catch(InvalidExtentException&) {
    // ignore
  }
  return *this;
}

void ImageHandle::Reciproce()
{
  throw Error("Reciproce not implemented");
  //ImageStatePtr()->Reciproce();
}

ExtentIterator ImageHandle::GetIterator() const {
   return ExtentIterator(GetExtent(), GetDomain());
}

// ConstImageHandle
ConstImageHandle::ConstImageHandle():
  handle_()
{
}

ConstImageHandle::ConstImageHandle(const ConstImageHandle& h):
  handle_(h.handle_)
{}

ConstImageHandle::ConstImageHandle(const ImageHandle& h):
  handle_(h)
{}

ConstImageHandle& ConstImageHandle::operator=(const ConstImageHandle& h)
{
  handle_=h.handle_;
  Notify();
  return *this;
}

ConstImageHandle& ConstImageHandle::operator=(const ImageHandle& h)
{
  handle_=h;
  Notify();
  return *this;
}

ImageHandle ConstImageHandle::Copy(bool cc) const
{
  return handle_.Copy(cc);
}

ImageHandle ConstImageHandle::Extract(const Extent& e) const
{
  return handle_.Extract(e);
}

bool ConstImageHandle::IsValid() const
{
  return handle_.IsValid();
}

long ConstImageHandle::MemSize() const
{
  return handle_.MemSize();
}

void ConstImageHandle::Apply(NonModAlgorithm& a) const 
{
  handle_.Apply(a);
}

void ConstImageHandle::ApplyIP(NonModAlgorithm& a) const 
{
  Apply(a);
}

ImageHandle ConstImageHandle::Apply(ModIPAlgorithm& a) const 
{
  return handle_.Apply(a);
}

ImageHandle ConstImageHandle::Apply(const ConstModIPAlgorithm& a) const 
{
  return handle_.Apply(a);
}

ImageHandle ConstImageHandle::Apply(ModOPAlgorithm& a) const 
{
  return handle_.Apply(a);
}

ImageHandle ConstImageHandle::Apply(const ConstModOPAlgorithm& a) const 
{
  return handle_.Apply(a);
}

void ConstImageHandle::StateApply(image_state::ImageStateNonModVisitorBase& v) const
{
  handle_.StateApply(v);
}

ImageHandle ConstImageHandle::StateApply(image_state::ImageStateModIPVisitorBase& v) const
{
  return handle_.StateApply(v);
}

ImageHandle ConstImageHandle::StateApply(const image_state::ImageStateConstModIPVisitorBase& v) const
{
  return handle_.StateApply(v);
}

ImageHandle ConstImageHandle::StateApply(image_state::ImageStateModOPVisitorBase& v) const
{
  return handle_.StateApply(v);
}

ImageHandle ConstImageHandle::StateApply(const image_state::ImageStateConstModOPVisitorBase& v) const
{
  return handle_.StateApply(v);
}

ImageHandle ConstImageHandle::StateApply(image_state::ImageStateMorphVisitorBase& v) const
{
  return handle_.StateApply(v);
}

void ConstImageHandle::Attach(DataObserver *o) const
{
  handle_.Attach(o);
}

void ConstImageHandle::Detach(DataObserver *o) const
{
  handle_.Detach(o);
}

void ConstImageHandle::Notify() const
{
  handle_.Notify();
}
void ConstImageHandle::Notify(const Extent& e) const
{
  handle_.Notify(e);
}
void ConstImageHandle::Notify(const Point& p) const
{
  handle_.Notify(p);
}

DataType ConstImageHandle::GetType() const 
{
  return handle_.GetType();
}

Extent ConstImageHandle::GetExtent() const 
{
  return handle_.GetExtent();
}

Point ConstImageHandle::GetSpatialOrigin() const 
{
  return handle_.GetSpatialOrigin();
}

Vec3 ConstImageHandle::GetAbsoluteOrigin() const
{
  return ImageStatePtr()->GetAbsoluteOrigin();
}

Vec3 ConstImageHandle::IndexToCoord(const Point &p) const
{
  return ImageStatePtr()->IndexToCoord(p);
}

Vec3 ConstImageHandle::CoordToIndex(const Vec3& c) const
{
  return ImageStatePtr()->CoordToIndex(c);
}

Vec3 ConstImageHandle::FractionalIndexToCoord(const Vec3 &p) const
{
  return ImageStatePtr()->FractionalIndexToCoord(p);
}

DataDomain ConstImageHandle::GetDomain() const 
{
  return handle_.GetDomain();
}

Real ConstImageHandle::GetReal(const Point& p) const 
{
  return handle_.GetReal(p);
}

Complex ConstImageHandle::GetComplex(const Point& p) const 
{
  return handle_.GetComplex(p);
}

Real ConstImageHandle::GetIntpolReal(const Vec3& v) const 
{
  return handle_.GetIntpolReal(v);
}
Real ConstImageHandle::GetIntpolReal(const Vec2& v) const 
{
  return handle_.GetIntpolReal(v);
}
Real ConstImageHandle::GetIntpolReal(const Real& d) const 
{
  return handle_.GetIntpolReal(d);
}

Complex ConstImageHandle::GetIntpolComplex(const Vec3& v) const 
{
  return handle_.GetIntpolComplex(v);
}
Complex ConstImageHandle::GetIntpolComplex(const Vec2& v) const 
{
  return handle_.GetIntpolComplex(v);
}
Complex ConstImageHandle::GetIntpolComplex(const Real& d) const 
{
  return handle_.GetIntpolComplex(d);
}

PixelSampling& ConstImageHandle::Sampling()
{
  return handle_.Sampling();
}

const PixelSampling& ConstImageHandle::Sampling() const
{
  return handle_.Sampling();
}

const ImageStateBasePtr& ConstImageHandle::ImageStatePtr() const
{
  return handle_.ImageStatePtr();
}

bool ConstImageHandle::operator==(const ConstImageHandle& ih) const
{
  return equal(ih);
}

bool ConstImageHandle::operator==(const ImageHandle& ih) const
{
  return equal(ih);
}

bool ConstImageHandle::operator!=(const ConstImageHandle& ih) const
{
  return !equal(ih);
}

bool ConstImageHandle::operator!=(const ImageHandle& ih) const
{
  return !equal(ih);
}


bool ConstImageHandle::equal(const ConstImageHandle& h) const
{
  return handle_==h.handle_;
}

bool ConstImageHandle::equal(const ImageHandle& h) const
{
  return handle_==h;
}

ExtentIterator ConstImageHandle::GetIterator() const {
   return handle_.GetIterator();
}

// global ops

bool operator==(const ImageHandle& lhs, const ConstImageHandle& rhs)
{
  return rhs.operator==(lhs);
}
bool operator!=(const ImageHandle& lhs, const ConstImageHandle& rhs)
{
  return rhs.operator!=(lhs);
}

ImageHandle operator+(const ConstImageHandle& h, Real v)
{
  return h.Copy()+=v;
}

ImageHandle operator+(const ConstImageHandle& h, const Complex& v)
{
  return h.Copy()+=v;
}

ImageHandle operator+(Real v, const ConstImageHandle& h)
{
  return h.Copy()+=v;
}

ImageHandle operator+(const Complex& v, const ConstImageHandle& h)
{
  return h.Copy()+=v;
}

ImageHandle operator-(const ConstImageHandle& h, Real v)
{
  return h.Copy()-=v;
}

ImageHandle operator-(const ConstImageHandle& h, const Complex& v)
{
  return h.Copy()-=v;
}

ImageHandle operator-(Real v, const ConstImageHandle& h)
{
  return h.Copy()-=v;
}

ImageHandle operator-(const Complex& v, const ConstImageHandle& h)
{
  return h.Copy()-=v;
}

ImageHandle operator*(const ConstImageHandle& h, Real v)
{
  return h.Copy()*=v;
}

ImageHandle operator*(const ConstImageHandle& h, const Complex& v)
{
  return h.Copy()*=v;
}

ImageHandle operator*(Real v, const ConstImageHandle& h)
{
  return h.Copy()*=v;
}

ImageHandle operator*(const Complex& v, const ConstImageHandle& h)
{
  return h.Copy()*=v;
}

ImageHandle operator/(const ConstImageHandle& h, Real v)
{
  return h.Copy()/=v;
}

ImageHandle operator/(const ConstImageHandle& h, const Complex& v)
{
  return h.Copy()/=v;
}

ImageHandle operator+(const ConstImageHandle& h1, const ConstImageHandle& h2)
{
  ImageHandle nh=h1.Copy();
  nh+=h2;
  return nh;
}

ImageHandle operator-(const ConstImageHandle& h1, const ConstImageHandle& h2)
{
  ImageHandle nh=h1.Copy();
  nh-=h2;
  return nh;
}

ImageHandle operator*(const ConstImageHandle& h1, const ConstImageHandle& h2)
{
  ImageHandle nh=h1.Copy();
  nh*=h2;
  return nh;
}

ImageHandle operator/(const ConstImageHandle& h1, const ConstImageHandle& h2)
{
  ImageHandle nh=h1.Copy();
  nh/=h2;
  return nh;
}



}} // namespace
