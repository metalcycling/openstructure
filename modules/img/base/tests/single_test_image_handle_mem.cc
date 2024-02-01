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

#include <iostream>

#include <boost/shared_ptr.hpp>

#include <ost/img/data.hh>
#include <ost/img/image_handle.hh>
#include <ost/img/image_impl.hh>
#include <ost/img/data_observer.hh>
#include <ost/img/observable.hh>

using namespace ::img;

typedef boost::shared_ptr<ImageStateBase> isbp;
typedef Observable<DataObserver> DataObservable;
typedef boost::shared_ptr<DataObservable> ObsPtr;

class IWrap: public Data {
public:
  IWrap(const detail::ImageImplPtr& iip):
    Data(),
    iip_(iip),
    obs_(new DataObservable())
  {}

  virtual DataType GetType() const {return REAL;}
  virtual DataDomain GetDomain() const {return SPATIAL;}
  virtual Extent GetExtent() const {return Extent();}
  virtual Point GetSpatialOrigin() const {return Point();}
  virtual Real GetReal(const Point& p) const {return 0.0;}
  virtual Complex GetComplex(const Point& p) const {return Complex(0.0,0.0);}
  virtual void Apply(NonModAlgorithm& a) const {}
  virtual void Attach(DataObserver* o) const {}
  virtual void Detach(DataObserver* o) const {}
  virtual void Notify() const {}
  virtual ImageHandle Crop(const Extent& e) const {return ImageHandle();}
  virtual PixelSampling& Sampling() {static PixelSampling ps; return ps;}
  virtual const PixelSampling& Sampling() const {static PixelSampling ps; return ps;}
  virtual void SetSpatialOrigin(const Point& p) {}
  

  void R() {
    iip_->RandomizeIP();
  }

  static IWrap Create(const Extent &e) {
    return IWrap(detail::ImageImpl::Create(e,REAL,SPATIAL));
  }

private:
  detail::ImageImplPtr iip_;
  ObsPtr obs_;
};

IWrap Create(const Extent& e)
{
  return IWrap::Create(e);
}


int main(int argc, char* argv[])
{
  for(int i=0;i<10000;++i) {
    std::cerr << i << std::endl;
    // variant 0
    /*
    isbp p(new RealSpatialImageState(Size(1000,1000)));
    p->RandomizeIP();
    */
    // variant 1 
    /*
    detail::ImageImplPtr iip(detail::ImageImpl::Create(Size(1000,1000),REAL,SPATIAL));
    iip->RandomizeIP();
    */
    // variant 2
    /*
    ImageHandle img = CreateImage(Size(1000,1000));
    img.RandomizeIP();
    */
    // variant 3

    IWrap im=Create(Size(1000,1000));
    im.R();
  }
}




