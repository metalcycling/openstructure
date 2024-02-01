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

#ifndef IMG_ALG_TRANSFORM_ALG_TRANSFORM_HH
#define IMG_ALG_TRANSFORM_ALG_TRANSFORM_HH

#include <ost/img/image_state.hh>
#include <ost/img/value_util.hh>
#include "transformations.hh"
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

struct DLLEXPORT_IMG_ALG TransformFnc {

  TransformFnc(): tf_() {}
  TransformFnc(const Transformation& tf, const Vec3 o=Vec3(0.0,0.0,0.0)): tf_(tf), offset_(o) {}

  template <typename T, class D>
  ImageStateBasePtr VisitState(const ImageStateImpl<T,D>& isi) const {
    /*
      old extent code, needs re-implementation
      Extent src_extent = isi.GetExtent();
      Extent dst_extent = (r==Transformer::ORIGINAL) ? src_extent : tf::Transform(src_extent,t);
      ImageStateImpl<V,DP> *nisi = new ImageStateImpl<V,DP>(dst_extent);

      Extent dst_extent = src_extent;
    */

    boost::shared_ptr<ImageStateImpl<T,D> > nisi = isi.CloneState(false);

    Mat4 imat = tf_.InverseMatrix();

    if(offset_==Vec3(0.0,0.0,0.0)) {
      for(ExtentIterator it(nisi->GetExtent()); !it.AtEnd(); ++it) {
        Point p(it);
        Vec3 vold_norm=Vec3(imat * Vec4(p.ToVec3()));
        nisi->Value(it) = isi.CalcIntpolValue(vold_norm);
      }
    } else {
      for(ExtentIterator it(nisi->GetExtent()); !it.AtEnd(); ++it) {
        Point p(it);
        Vec3 vold_norm=Vec3(imat * Vec4(p.ToVec3()-offset_))+offset_;
        nisi->Value(it) = isi.CalcIntpolValue(vold_norm);
      }
    }
    
    return nisi;
  }

  static String GetAlgorithmName() {return "Transform";}

private:
  Transformation tf_;
  Vec3 offset_;
};

typedef ImageStateConstModOPAlgorithm<TransformFnc> Transform;

}}} // ns

#endif

