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
#include "polar.hh"
#include <ost/img/image_state.hh>
#include <ost/img/algorithm.hh>
#include <ost/geom/vecmat2_op.hh>

namespace ost { namespace  img { namespace alg {

namespace detail {

class PolarImplFnc {
public:
  PolarImplFnc(Real s=1.0): sampling_(s) {}

  static String GetAlgorithmName() {return "Polar";}

  template <typename T, class D>
  ImageStateBasePtr VisitState(const img::ImageStateImpl<T,D>& isi) {
    Size s=isi.GetExtent().GetSize();
//    unsigned int rmax=static_cast<unsigned int>(std::min(s[0],s[1])*0.5);
    unsigned int rmax=static_cast<unsigned int>(round(sqrt(s[0]*s[0]+s[1]*s[1])*0.5));
    Extent e(Size(rmax,static_cast<unsigned int>(rmax*2.0*M_PI*sampling_)));
    std::cerr << e << std::endl;
    boost::shared_ptr<ImageStateImpl<T,image_state::SpatialDomain> > result(new ImageStateImpl<T,image_state::SpatialDomain>(e,isi.GetSampling()));

    Real prefac=1.0/(rmax*sampling_);

    for(ExtentIterator it(e);!it.AtEnd();++it) {
      Vec2 v(it[0],0);
      v=Rotate(v,it[1]*prefac);
      result->Value(it)=isi.CalcIntpolValue(Vec3(v));
    }
    return result;
  }
private:
  Real sampling_;
};

typedef image_state::ImageStateModOPAlgorithm<PolarImplFnc> PolarImpl;

} // ns

Polar::Polar():
  ModOPAlgorithm("Polar"),
  sampling_(1.0)
{}

Polar::Polar(Real sampling):
  ModOPAlgorithm("Polar"),
  sampling_(sampling)
{}

ImageHandle Polar::Visit(const ConstImageHandle& i)
{
#if 0
  PolarImpl impl(sampling_);
  return i.Apply(impl);
#else
  Size s=i.GetExtent().GetSize();
//  unsigned int rmax=static_cast<unsigned int>(std::min(s[0],s[1])*0.5);
  unsigned int rmax=static_cast<unsigned int>(round(sqrt(static_cast<double>(s[0]*s[0]+s[1]*s[1]))*0.5));
  Extent e(Size(rmax,static_cast<unsigned int>(rmax*2.0*M_PI*sampling_)));
  ImageHandle result = CreateImage(e,i.GetType(),SPATIAL);
  Real prefac=1/(rmax*sampling_);
  for(ExtentIterator it(e);!it.AtEnd();++it) {
    Vec2 v(it[0],0);
    v=Rotate(v,it[1]*prefac);
    result.SetComplex(it,i.GetIntpolComplex(Vec3(v)));
  }
  return result;
#endif
}

EllipticPolar::EllipticPolar(Real ratio,Real angle,Real sampling):
  ModOPAlgorithm("EllipticPolar"),
  ratio_(ratio),
  angle_(angle),
  sampling_(sampling)
{}


ImageHandle EllipticPolar::Visit(const ConstImageHandle& i)
{
  Ellipse2 ellipse(Vec2(0,0),ratio_,1,angle_);
  Size s=i.GetExtent().GetSize();
  Real ang=atan(static_cast<double>(s[1]/s[0]));
  Real mdist=std::min(Length(ellipse.AtAngle(ang)),Length(ellipse.AtAngle(-ang)));
  unsigned int rmax=static_cast<unsigned int>(round(sqrt(static_cast<double>(s[0]*s[0]+s[1]*s[1]))*0.5));
  Extent e(Size(rmax,static_cast<unsigned int>(rmax*2.0*M_PI*sampling_)));
  ImageHandle result = CreateImage(e,i.GetType(),SPATIAL);
  Real prefac=1/(rmax*sampling_);
  for(ExtentIterator it(e);!it.AtEnd();++it) {
    Ellipse2 ellipse(Vec2(0,0),ratio_*it[0]/mdist,it[0]/mdist,angle_);
    result.SetComplex(it,i.GetIntpolComplex(Vec3(ellipse.At(it[1]*prefac))));
  }
  return result;
}

Cartesian::Cartesian():
  ModOPAlgorithm("Cartesian")
{}


ImageHandle Cartesian::Visit(const ConstImageHandle& i)
{
  Size s=i.GetExtent().GetSize();
  Extent e(Size(2*s[0],2*s[0]),Point(0,0));
  ImageHandle result = CreateImage(e,i.GetType(),SPATIAL);
  Real prefac=s[1]*0.5/M_PI;
  for(ExtentIterator it(e);!it.AtEnd();++it) {
    Vec2 v(sqrt(static_cast<double>(it[0]*it[0]+it[1]*it[1])),
           (fmod(atan2(static_cast<double>(it[1]),
                       static_cast<double>(it[0]))
                 +2.0*M_PI,2.0*M_PI))*prefac);

    result.SetComplex(it,i.GetIntpolComplex(Vec3(v)));
  }
  return result;
}

}}} // ns
