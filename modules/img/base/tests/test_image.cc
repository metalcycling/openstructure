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
  Author: Ansgar Philippsen
*/

#include <iostream>

#include "test_image.hh"
#include <ost/invalid_handle.hh>
#include <ost/img/image.hh>
#include <ost/img/data_observer.hh>
#include <ost/img/observable.hh>
#include <ost/img/algorithm.hh>

#include <ost/img/alg/randomize.hh>

using namespace ost::img;

/*
  black-box image tests
*/

namespace test_image {

template <typename T>
bool check_close(const T& t1, const T& t2, Real eps);

template <>
bool check_close<Vec3>(const Vec3& t1, const Vec3& t2, Real eps)
{
  return 
    std::abs(t1[0]-t2[0])<eps &&
    std::abs(t1[1]-t2[1])<eps &&
    std::abs(t1[2]-t2[2])<eps;
}

template <typename T>
bool check_close(const T& t1, const T& t2, Real eps)
{
  return std::abs(t1-t2)<eps;
}

void test_Creation()
{
  // all default
  ImageHandle ih0=CreateImage();
  BOOST_CHECK(ih0.GetExtent()==Extent(Point(0,0,0),Point(0,0,0)));
  BOOST_CHECK(ih0.GetType()==REAL);
  BOOST_CHECK(ih0.GetDomain()==SPATIAL);

  // extent
  ImageHandle ih1=CreateImage(Size(5,6,7));
  BOOST_CHECK(ih1.GetExtent()==Extent(Point(0,0,0),Point(4,5,6)));
  BOOST_CHECK(ih1.GetType()==REAL);
  BOOST_CHECK(ih1.GetDomain()==SPATIAL);

  ih1=CreateImage( (Extent(Point(-1,-2,-3),Size(5,6,7))) );
  BOOST_CHECK(ih1.GetExtent()==Extent(Point(-1,-2,-3),Point(3,3,3)));
  BOOST_CHECK(ih1.GetType()==REAL);
  BOOST_CHECK(ih1.GetDomain()==SPATIAL);

  // extent and type
  ImageHandle ih2=CreateImage(Size(3,1,4), COMPLEX);
  BOOST_CHECK(ih2.GetExtent()==Extent(Point(0,0,0),Point(2,0,3)));
  BOOST_CHECK(ih2.GetType()==COMPLEX);
  BOOST_CHECK(ih2.GetDomain()==SPATIAL);

  // extent, type and domain
  ImageHandle ih3=CreateImage(Size(6,6,6), COMPLEX, FREQUENCY);
  BOOST_CHECK(ih3.GetExtent()==Extent(Point(-2,-2,-2),Point(3,3,3)));
  BOOST_CHECK(ih3.GetType()==COMPLEX);
  BOOST_CHECK(ih3.GetDomain()==FREQUENCY);
}

void test_Cctor_Assignement()
{
  ImageHandle ih1=CreateImage(Size(10,10,10));
  ih1.ApplyIP(alg::Randomize());

  ImageHandle ih2(ih1);
  BOOST_CHECK(ih1==ih2);
  for(ExtentIterator it(ih1.GetExtent());!it.AtEnd(); ++it) {
    BOOST_CHECK(ih1.GetReal(it)==ih2.GetReal(it));
  }

  ImageHandle ih3=CreateImage(Size(10));
  ih3.ApplyIP(alg::Randomize());
  ih3=ih1;
  BOOST_CHECK(ih1==ih3);
  for(ExtentIterator it(ih1.GetExtent());!it.AtEnd(); ++it) {
    BOOST_CHECK(ih1.GetReal(it)==ih3.GetReal(it));
  }
}

void test_ConstHandle()
{
  ImageHandle i1 = CreateImage(Size(10));
  i1.ApplyIP(alg::Randomize());
  ConstImageHandle c1(i1);
  BOOST_CHECK(c1.GetExtent()==i1.GetExtent());
  BOOST_CHECK(c1.GetReal(Point(0))==i1.GetReal(Point(0)));
  i1.SetReal(Point(1),2.0);
  BOOST_CHECK(c1.GetReal(Point(1))==2.0);

  ImageHandle i2 = CreateImage(Size(4));
  i2.ApplyIP(alg::Randomize());
  c1=i2;
  BOOST_CHECK(c1.GetExtent()==i2.GetExtent());
  BOOST_CHECK(c1.GetReal(Point(0))==i2.GetReal(Point(0)));
  i2.SetReal(Point(1),3.0);
  BOOST_CHECK(c1.GetReal(Point(1))==3.0);
}


void test_Reset()
{
  ImageHandle ih=CreateImage();
  ih.Reset(Extent(Size(5,6,7)));
  BOOST_CHECK(ih.GetExtent()==Extent(Point(0,0,0),Point(4,5,6)));
  BOOST_CHECK(ih.GetType()==REAL);
  BOOST_CHECK(ih.GetDomain()==SPATIAL);

  ih.Reset(Extent(Size(6,6,6)), COMPLEX, FREQUENCY);
  BOOST_CHECK(ih.GetExtent()==Extent(Point(-2,-2,-2),Point(3,3,3)));
  BOOST_CHECK(ih.GetType()==COMPLEX);
  BOOST_CHECK(ih.GetDomain()==FREQUENCY);
}

class MyObserver: public DataObserver {
public:
  MyObserver(const ImageHandle& d):
    DataObserver(d),
    update_count(0),
    release_count(0)
  {}

  virtual void ObserverUpdate() {
    update_count++;
  }

  virtual void ObserverRelease() {
    release_count++;
  }

  void CheckData() {
    GetObservedData();
  }

  int update_count;
  int release_count;
};

void test_Observer()
{
  MyObserver* to1=0;
  MyObserver* to2=0;
  
  {
    ImageHandle ih = CreateImage();
    BOOST_CHECK(ih.obs_);

    BOOST_CHECK(ih.obs_.get()->GetListSize()==0);

    to1 = new MyObserver(ih);
    BOOST_CHECK(ih.obs_.get()->GetListSize()==1);
    to2 = new MyObserver(ih);
    BOOST_CHECK(ih.obs_.get()->GetListSize()==2);

    BOOST_CHECK(to1->update_count==0);
    BOOST_CHECK(to1->release_count==0);
    BOOST_CHECK(to2->update_count==0);
    BOOST_CHECK(to2->release_count==0);
    
    ih.Notify();

    BOOST_CHECK(to1->update_count==1);
    BOOST_CHECK(to2->update_count==1);

    delete to2;

    BOOST_CHECK(ih.obs_.get()->GetListSize()==1);

    ImageHandle ih2 = ih;

    BOOST_CHECK(ih.obs_.get() == ih2.obs_.get());

    ih2.Notify();

    BOOST_CHECK(to1->update_count==2);
  }
  // ih and ih2 go out-of-scope

  // observer retains internal reference, therefore no call to release should happen
  BOOST_CHECK(to1->release_count==0);

  BOOST_CHECK_NO_THROW(to1->CheckData());

  delete to1;
}

void test_Interpolation()
{
  Real dx=0.25,dy=0.75,dz=0.125;
  Real one=1.0;
  ImageHandle ih3d=CreateImage(Extent(Point(0,0,0),Point(1,1,1)));
  
  ih3d.SetReal(Point(0,0,0), Real(0.000));
  ih3d.SetReal(Point(1,0,0), Real(0.500));
  ih3d.SetReal(Point(0,1,0), Real(0.250));
  ih3d.SetReal(Point(0,0,1), Real(0.125));
  ih3d.SetReal(Point(1,0,1), Real(0.500));
  ih3d.SetReal(Point(0,1,1), Real(0.250));
  ih3d.SetReal(Point(1,1,0), Real(1.000));  
  ih3d.SetReal(Point(1,1,1), Real(2.000));

  ih3d.ApplyIP(alg::Randomize());
  
  Real val3= (one-dx)*(one-dy)*(one-dz)*ih3d.GetReal(Point(0,0,0))
              +(  dx  )*(one-dy)*(one-dz)*ih3d.GetReal(Point(1,0,0))  
              +(one-dx)*(  dy  )*(one-dz)*ih3d.GetReal(Point(0,1,0))  
              +(  dx  )*(  dy  )*(one-dz)*ih3d.GetReal(Point(1,1,0))  
              +(one-dx)*(one-dy)*(  dz  )*ih3d.GetReal(Point(0,0,1))
              +(  dx  )*(one-dy)*(  dz  )*ih3d.GetReal(Point(1,0,1))  
              +(one-dx)*(  dy  )*(  dz  )*ih3d.GetReal(Point(0,1,1))  
              +(  dx  )*(  dy  )*(  dz  )*ih3d.GetReal(Point(1,1,1));
  BOOST_CHECK(check_close(val3,ih3d.GetIntpolReal(Vec3(dx,dy,dz)),1.0e-6));


  ImageHandle ih2d=CreateImage(Extent(Point(0,0),Point(1,1)));
  ih3d.SetReal(Point(0,0,0), Real(0.000));
  ih3d.SetReal(Point(1,0,0), Real(0.500));
  ih3d.SetReal(Point(0,1,0), Real(0.250));
  ih3d.SetReal(Point(1,1,0), Real(0.125));
  Real val2= (one-dx)*(one-dy)*ih2d.GetReal(Point(0,0))
              +(  dx  )*(one-dy)*ih2d.GetReal(Point(1,0))  
              +(one-dx)*(  dy  )*ih2d.GetReal(Point(0,1))  
              +(  dx  )*(  dy  )*ih2d.GetReal(Point(1,1));  
  
  BOOST_CHECK(check_close(val2,ih2d.GetIntpolReal(Vec2(dx,dy)),1.0e-6));

  ImageHandle ih1d=CreateImage(Extent(Point(0),Point(1)));
  ih3d.SetReal(Point(0,0,0), Real(0.000));
  ih3d.SetReal(Point(1,0,0), Real(0.500));
  Real val1= (one-dx)*ih1d.GetReal(Point(0,0))
              +(  dx  )*ih1d.GetReal(Point(1,0));  
  
  BOOST_CHECK(check_close(val1,ih1d.GetIntpolReal(dx),1.0e-6));

}

void test_ScalarOps()
{
  ImageHandle ih1=CreateImage(Size(10,10));
  ImageHandle ih2 = ih1.Copy();
  ih2+=1.0;

  for(ExtentIterator it(ih2.GetExtent()); !it.AtEnd(); ++it) {
    BOOST_CHECK(ih2.GetReal(it) == ih1.GetReal(it) + 1.0);
  }
    
  ih2 = ih1.Copy();
  ih2-=2.5;

  for(ExtentIterator it(ih2.GetExtent()); !it.AtEnd(); ++it) {
    BOOST_CHECK(ih2.GetReal(it) == ih1.GetReal(it) - 2.5);
  }
    
  ih2 = ih1.Copy();
  ih2*=0.9;

  for(ExtentIterator it(ih2.GetExtent()); !it.AtEnd(); ++it) {
    BOOST_CHECK(std::abs(ih2.GetReal(it) - ih1.GetReal(it) * 0.9)<1e-10);
  }

  ih2 = ih1.Copy();
  ih2/=2.1;

  for(ExtentIterator it(ih2.GetExtent()); !it.AtEnd(); ++it) {
    BOOST_CHECK(std::abs(ih2.GetReal(it) - ih1.GetReal(it) / 2.1)<1e-10);
  }
}

void test_ImageOps()
{
  Extent ex1(Point(-3,-2),Point(2,2));
  ImageHandle h1=CreateImage(ex1);
  Real count=1.0;
  for(ExtentIterator it(ex1); !it.AtEnd(); ++it) {
    h1.SetReal(it,count);
    count+=1.0;
  }
  Extent ex2(Point(-1,-2),Point(3,4));
  ImageHandle h2=CreateImage(ex2, COMPLEX);
  for(ExtentIterator it(ex2); !it.AtEnd(); ++it) {
    h2.SetReal(it,count);
    count+=1.0;
  }

  // addition
  ImageHandle h3=h1+h2;
  for(ExtentIterator it(Extent(Point(-4,-3),Point(4,5))); !it.AtEnd(); ++it) {
    if(ex1.Contains(it) && ex2.Contains(it)) {
      BOOST_CHECK_SMALL(h3.GetReal(it)-(h1.GetReal(it)+h2.GetReal(it)),Real(1e-6));
    } else if (ex1.Contains(it)) {
      BOOST_CHECK(h3.GetReal(it)==h1.GetReal(it));
    } else {
      BOOST_CHECK(h3.GetReal(it)==0.0);
    }
  }

  // subtraction
  h3=h1-h2;
  for(ExtentIterator it(Extent(Point(-2,-1),Point(1,2))); !it.AtEnd(); ++it) {
    if(ex1.Contains(it) && ex2.Contains(it)) {
      BOOST_CHECK_SMALL(h3.GetReal(it)-(h1.GetReal(it)-h2.GetReal(it)),Real(1e-6));
    } else if (ex1.Contains(it)) {
      BOOST_CHECK(h3.GetReal(it)==h1.GetReal(it));
    } else {
      BOOST_CHECK(h3.GetReal(it)==0.0);
    }
  }

  // half frequency tests
  h1=CreateImage(Size(4,5),COMPLEX,HALF_FREQUENCY);
  h2=CreateImage(Size(4,5),COMPLEX,HALF_FREQUENCY);
  for(ExtentIterator it(h1.GetExtent()); !it.AtEnd(); ++it) {
    Complex c(count,count+2);
    h1.SetComplex(it,c);
    count+=1.0;
  }
  for(ExtentIterator it(h2.GetExtent()); !it.AtEnd(); ++it) {
    Complex c(count,count+2);
    h2.SetComplex(it,c);
    count+=1.0;
  }

  h3=h1+h2;
  for(ExtentIterator it(h1.GetExtent());!it.AtEnd();++it) {
    BOOST_CHECK_SMALL(std::abs(h1.GetComplex(it)+h2.GetComplex(it)-h3.GetComplex(it)),Real(1e-10));
  }

  h3=h1-h2;
  for(ExtentIterator it(Extent(Point(-1,0),Point(2,2)));!it.AtEnd();++it) {
    BOOST_CHECK_SMALL(std::abs(h1.GetComplex(it)-h2.GetComplex(it)-h3.GetComplex(it)),Real(1e-10));
  }
}

void test_Empty()
{
  ImageHandle ih;

  BOOST_CHECK_THROW(ih.GetExtent(),InvalidImageHandle);
}


void test_Origin() 
{
  ImageHandle rs = CreateImage(Extent(Size(6,7),Point(0,0)));
  std::ostringstream msg;

  msg << rs.GetExtent() << " != " << Extent(Point(-2,-3),Size(6,7));
  BOOST_CHECK_MESSAGE(rs.GetExtent()==Extent(Point(-2,-3),Point(3,3)),msg.str());

  msg << rs.GetSpatialOrigin() << " != " << Point(-2,-3);
  BOOST_CHECK_MESSAGE(rs.GetSpatialOrigin()==Point(-2,-3),msg.str());
  
  rs.SetSpatialOrigin(Point(0,0));
  msg.str("");
  msg << rs.GetSpatialOrigin() << " != " << Point(0,0);
  BOOST_CHECK_MESSAGE(rs.GetSpatialOrigin()==Point(0,0),msg.str());
  msg.str("");
  msg << rs.GetExtent() << " != " << Extent(Point(0,0),Size(6,7));
  BOOST_CHECK_MESSAGE(rs.GetExtent()==Extent(Point(0,0),Size(6,7)),msg.str());

  rs.CenterSpatialOrigin();
  msg.str("");
  msg << rs.GetSpatialOrigin() << " != " << Point(-2,-3);
  BOOST_CHECK_MESSAGE(rs.GetSpatialOrigin()==Point(-2,-3),msg.str());
}

void test_DataAlgorithms()
{
  class TestIPAlg: public ModIPAlgorithm {
  public:
    TestIPAlg():
      ModIPAlgorithm("") {}
    virtual void Visit(ImageHandle& h) {
      h.SetReal(Point(0),1.0);
    }
  };

  class TestOPAlg: public ModOPAlgorithm {
  public:
    TestOPAlg():
      ModOPAlgorithm("") {}
    virtual ImageHandle Visit(const ConstImageHandle& h) {
      ImageHandle nh=h.Copy(true);
      nh.SetReal(Point(0),-1.0);
      return nh;
    }
  };

  TestIPAlg ip_alg;
  TestOPAlg op_alg;

  ImageHandle ih1=CreateImage(Size(1));
  BOOST_CHECK(ih1.GetReal(Point(0))==0.0);

  // in-place alg applied in-place
  ih1.ApplyIP(ip_alg);
  BOOST_CHECK(ih1.GetReal(Point(0))==1.0);

  // out-of-place alg applied out-of-place
  ImageHandle ih2=ih1.Apply(op_alg);
  BOOST_CHECK(ih1.GetReal(Point(0))==1.0);
  BOOST_CHECK(ih2.GetReal(Point(0))==-1.0);

  // out-of-place alg applied in-place
  ih1.ApplyIP(op_alg);
  BOOST_CHECK(ih1.GetReal(Point(0))==-1.0);

  // in-place alg applied out-of-place
  ImageHandle ih3=ih2.Apply(ip_alg);
  BOOST_CHECK(ih2.GetReal(Point(0))==-1.0);
  BOOST_CHECK(ih3.GetReal(Point(0))==1.0);
}

void test_Huge()
{
  ImageHandle h1=CreateImage(Size(200,200,200));
  {
    ImageHandle h2;
    h2=h1;
  }
}

void test_Mult() 
{
  ImageHandle im1=CreateImage(Extent(Point(-2,-1,-3),Point(1,2,1)),REAL);
  ImageHandle im2=CreateImage(Extent(Point(-1,-2,-1),Point(2,1,3)),COMPLEX);
  Real count=0;
  for(ExtentIterator it(im1.GetExtent()); !it.AtEnd(); ++it) {
    im1.SetReal(it,count);
    count+=1.0;
  }
  for(ExtentIterator it(im2.GetExtent()); !it.AtEnd(); ++it) {
    Complex c(count,count+2);
    im2.SetComplex(it,c);
    count+=1.0;
  }

  ImageHandle im12=im1*im2;
  BOOST_CHECK(im12.GetExtent()==im1.GetExtent());

  for(ExtentIterator it(im12.GetExtent());!it.AtEnd();++it) {
    if(im1.GetExtent().Contains(it) && im2.GetExtent().Contains(it)) {
      BOOST_CHECK_CLOSE(im12.GetReal(it),im1.GetReal(it)*im2.GetReal(it),Real(1e-10));
    } else {
      //BOOST_CHECK(im21.GetReal(it)==0.0);
    }
  }

  ImageHandle im21=im2*im1;
  BOOST_CHECK(im21.GetExtent()==im2.GetExtent());

  for(ExtentIterator it(im21.GetExtent());!it.AtEnd();++it) {
    if(im1.GetExtent().Contains(it) && im2.GetExtent().Contains(it)) {
        BOOST_CHECK_SMALL(std::abs(im21.GetComplex(it)-im2.GetComplex(it)*im1.GetComplex(it)),Real(1e-10));
    } else {
      //BOOST_CHECK(im21.GetReal(it)==0.0);
    }
  }
}

void test_AddSub() 
{
  ImageHandle im1=CreateImage(Extent(Point(-2,-1,-3),Point(1,2,1)),REAL);
  ImageHandle im2=CreateImage(Extent(Point(-1,-2,-1),Point(2,1,3)),COMPLEX);
  im1.ApplyIP(alg::Randomize());
  im2.ApplyIP(alg::Randomize());

  ImageHandle im12=im1+im2;
  BOOST_CHECK(im12.GetExtent()==im1.GetExtent());

  for(ExtentIterator it(im12.GetExtent());!it.AtEnd();++it) {
    if(im2.GetExtent().Contains(it)) {
      BOOST_CHECK_CLOSE(im12.GetReal(it),im1.GetReal(it)+im2.GetReal(it),Real(1e-10));
    } else {
      BOOST_CHECK(im12.GetReal(it)==im1.GetReal(it));
    }
  }

  im12=im1-im2;
  BOOST_CHECK(im12.GetExtent()==im1.GetExtent());

  for(ExtentIterator it(im12.GetExtent());!it.AtEnd();++it) {
    if(im2.GetExtent().Contains(it)) {
      BOOST_CHECK_CLOSE(im12.GetReal(it),im1.GetReal(it)-im2.GetReal(it),Real(1e-10));
    } else {
      BOOST_CHECK(im12.GetReal(it)==im1.GetReal(it));
    }
  }

  ImageHandle im21=im2+im1;
  BOOST_CHECK(im21.GetExtent()==im2.GetExtent());

  for(ExtentIterator it(im21.GetExtent());!it.AtEnd();++it) {
    if(im2.GetExtent().Contains(it)) {
      BOOST_CHECK_CLOSE(real(im21.GetComplex(it)),real(im1.GetComplex(it)+im2.GetComplex(it)),Real(1e-10));
      BOOST_CHECK_CLOSE(imag(im21.GetComplex(it)),imag(im1.GetComplex(it)+im2.GetComplex(it)),Real(1e-10));
    } else {
      BOOST_CHECK_CLOSE(im21.GetReal(it),im2.GetReal(it),Real(1e-10));
    }
  }

  im21=im2-im1;
  BOOST_CHECK(im21.GetExtent()==im2.GetExtent());

  for(ExtentIterator it(im21.GetExtent());!it.AtEnd();++it) {
    if(im2.GetExtent().Contains(it)) {
      BOOST_CHECK_CLOSE(real(im21.GetComplex(it)),real(im2.GetComplex(it)-im1.GetComplex(it)),Real(1e-10));
      BOOST_CHECK_CLOSE(imag(im21.GetComplex(it)),imag(im2.GetComplex(it)-im1.GetComplex(it)),Real(1e-10));
    } else {
      BOOST_CHECK_CLOSE(im21.GetReal(it),im2.GetReal(it),Real(1e-10));
    }
  }
}


void test_Copy()
{
  // real spatial
  ImageHandle i1=CreateImage(Extent(Point(0,0),Point(9,9)));
  i1.SetSpatialOrigin(Point(7,3,-2));
  i1.ApplyIP(alg::Randomize());
  ImageHandle i2=i1.Copy();
  BOOST_CHECK(i2.GetSpatialOrigin()==Point(7,3,-2));
  ImageHandle i3=i1.Copy(false);
  BOOST_CHECK(i3.GetSpatialOrigin()==Point(7,3,-2));
  for(ExtentIterator it(i1.GetExtent());!it.AtEnd();++it) {
    BOOST_CHECK(i1.GetReal(it)==i2.GetReal(it));
    BOOST_CHECK(i3.GetReal(it)==0.0);
  }
  // complex half frequency
  i1=CreateImage(Size(6,7),COMPLEX,HALF_FREQUENCY);
  i1.SetSpatialOrigin(Point(-3,4));
  i2=i1.Copy();
  BOOST_CHECK(i1.GetExtent()==i2.GetExtent());
  BOOST_CHECK(i2.GetSpatialOrigin()==i2.GetSpatialOrigin());
}

void test_Extract()
{
  ImageHandle i1=CreateImage(Extent(Point(0,0),Point(9,9)));
  Vec3 ps(1.0,2.0,3.0);
  i1.SetPixelSampling(ps);
  i1.ApplyIP(alg::Randomize());
  ImageHandle i2 = i1.Extract(Extent(Point(2,3),Point(7,8)));
  BOOST_CHECK(i1.GetPixelSampling()==i2.GetPixelSampling());
  BOOST_CHECK(i2.GetDomain()==SPATIAL);
  BOOST_CHECK(i2.GetExtent()==Extent(Point(2,3),Point(7,8)));
  for(ExtentIterator it(i2.GetExtent());!it.AtEnd();++it) {
    BOOST_CHECK(i2.GetReal(it)==i1.GetReal(it));
  }

  i1=CreateImage(Size(13,20),COMPLEX,FREQUENCY);
  i1.SetPixelSampling(ps);
  i1.ApplyIP(alg::Randomize());
  i2 = i1.Extract(Extent(Point(2,3),Point(7,8)));
  BOOST_CHECK(i1.GetPixelSampling()==i2.GetPixelSampling());
  BOOST_CHECK(i2.GetDomain()==SPATIAL);
  BOOST_CHECK(i2.GetExtent()==Extent(Point(2,3),Point(7,8)));
  for(ExtentIterator it(i2.GetExtent());!it.AtEnd();++it) {
    BOOST_CHECK(i2.GetComplex(it)==i1.GetComplex(it));
  }
  
  i1=CreateImage(Size(13,30),COMPLEX,HALF_FREQUENCY);
  i1.SetPixelSampling(ps);
  i1.ApplyIP(alg::Randomize());
  i2 = i1.Extract(Extent(Point(-2,-3),Point(7,8)));
  BOOST_CHECK(i1.GetPixelSampling()==i2.GetPixelSampling());
  BOOST_CHECK(i2.GetDomain()==SPATIAL);
  BOOST_CHECK(i2.GetExtent()==Extent(Point(-2,-3),Point(7,8)));
  for(ExtentIterator it(i2.GetExtent());!it.AtEnd();++it) {
    BOOST_CHECK(i2.GetComplex(it)==i1.GetComplex(it));
  }
  
}

void test_Paste()
{
  ImageHandle i1=CreateImage(Extent(Point(-3,-4),Point(2,5)));
  ImageHandle i2=CreateImage(Extent(Point(0,0),Point(2,2)));
  i1.ApplyIP(alg::Randomize());
  i2.ApplyIP(alg::Randomize());
  ImageHandle i3=i1.Copy();
  i3.Paste(i2);
  
  for(ExtentIterator it(i3.GetExtent());!it.AtEnd();++it) {
    if(i2.GetExtent().Contains(it)) {
      BOOST_CHECK(i3.GetReal(it)==i2.GetReal(it));
    } else {
      BOOST_CHECK(i3.GetReal(it)==i1.GetReal(it));
    }
  }
}

void test_PointToCoord() {
  ImageHandle spatial = CreateImage(Extent(Size(10,10,10)));
  spatial.SetSpatialSampling(2.0);
  ImageHandle frequency = CreateImage(Extent(Size(10,10,10)),COMPLEX,FREQUENCY);
  frequency.SetSpatialSampling(0.1);
  spatial.SetAbsoluteOrigin(Vec3(5.0,5.0,5.0));
  frequency.SetAbsoluteOrigin(Vec3(5.0,5.0,5.0));
  BOOST_CHECK(spatial.IndexToCoord(Point(2.0,2.0,2.0)) == Vec3(9.0,9.0,9.0));
  BOOST_CHECK(frequency.IndexToCoord(Point(3.0,3.0,3.0)) == Vec3(3.0,3.0,3.0));
  BOOST_CHECK(frequency.CoordToIndex(Vec3(2.0,2.0,2.0)) == Vec3(2.0,2.0,2.0));
  BOOST_CHECK(spatial.CoordToIndex(Vec3(11.0,11.0,11.0)) == Vec3(3.0,3.0,3.0));
}

void test_ThrowInvalidImageHandle() {
  ImageHandle image;
  BOOST_CHECK_THROW(ost::CheckHandleValidity(image), ost::InvalidHandle);
  image=CreateImage(Size(1,1));
  BOOST_CHECK_NO_THROW(ost::CheckHandleValidity(image));
}

} // ns

test_suite* CreateImageTest()
{
  using namespace test_image;
  test_suite* ts=BOOST_TEST_SUITE("Image Test");

  ts->add(BOOST_TEST_CASE(&test_Creation));
  ts->add(BOOST_TEST_CASE(&test_Cctor_Assignement));
  ts->add(BOOST_TEST_CASE(&test_ConstHandle));
  ts->add(BOOST_TEST_CASE(&test_Reset));
  ts->add(BOOST_TEST_CASE(&test_Observer));
  ts->add(BOOST_TEST_CASE(&test_ScalarOps));
  ts->add(BOOST_TEST_CASE(&test_ImageOps));
  ts->add(BOOST_TEST_CASE(&test_Empty));
  ts->add(BOOST_TEST_CASE(&test_Origin));
  ts->add(BOOST_TEST_CASE(&test_DataAlgorithms));
  ts->add(BOOST_TEST_CASE(&test_Huge));
  ts->add(BOOST_TEST_CASE(&test_Copy));
  ts->add(BOOST_TEST_CASE(&test_Extract));
  ts->add(BOOST_TEST_CASE(&test_Paste));
  ts->add(BOOST_TEST_CASE(&test_ThrowInvalidImageHandle));
  ts->add(BOOST_TEST_CASE(&test_Mult));
  ts->add(BOOST_TEST_CASE(&test_AddSub));
  ts->add(BOOST_TEST_CASE(&test_Interpolation));
  ts->add(BOOST_TEST_CASE(&test_PointToCoord));
  return ts;
}
