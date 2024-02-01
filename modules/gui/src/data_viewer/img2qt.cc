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

#include <ost/img/data.hh>
#include <ost/img/algorithm.hh>
#include <ost/img/image_state.hh>
#include <ost/stdint.hh>
#include "img2qt.hh"

#include <QImage>
using namespace ost::img;


namespace ost { namespace img { namespace gui {

namespace {

template <typename T>
T clamp(const T& v, const T& a, const T& b)
{
  return std::min<T>(b,std::max<T>(v,a));
}

/*
  this is an attempt to speed up the data viewer as much as possible;
  qualitiative speed measurements point to isi.Value(Point()) as the 
  bottleneck though, and that cannot be sped up without caching the
  values upon translation
*/
class ConvBase0
{
public:
  ConvBase0(QImage* img, const geom::Vec3& offset, int z, 
            const NormalizerPtr& norm):
    img_(img), foff_(offset), z_(z), norm_(norm), off_(offset)
  {
  }

  // image state algorithm interface
  template <typename T, class D>
  void VisitState(const ImageStateImpl<T,D>& isi)
  {
    Extent iex(off_,Size(img_->size().width(),img_->size().height()));
    if(!HasOverlap(iex,isi.GetExtent())) return;
    Extent area=Overlap(iex,isi.GetExtent());
    
    for(int y=area.GetStart()[1];y<=area.GetEnd()[1];++y) {
      for(int x=area.GetStart()[0];x<=area.GetEnd()[0];++x) {
        int cv=clamp(static_cast<int>(norm_->Convert(Val2Val<T,Real>(isi.Value(Point(x,y,z_))))),0,255);
        img_->setPixel(x-off_[0],y-off_[1],qRgb(cv,cv,cv));
      }
    }

  }


  static String GetAlgorithmName() {return "Conv0";}


private:

  QImage* img_;
  geom::Vec3 foff_;
  int z_;
  NormalizerPtr norm_;
  Point off_;
};

} // anon ns


void img2qt(const Data& d, QImage* img,
             int logscale, const geom::Vec3& offset, int z,
             const NormalizerPtr& norm)
{
  if(logscale==0) {
    ImageStateNonModAlgorithm<ConvBase0> conv(img, offset, z, norm);
    d.Apply(conv);
  } else if(logscale<0) {
  } else {
  }
  
}


}}}  //ns
