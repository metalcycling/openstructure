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
  Authors: Ansgar Philippsen, Andreas Schenk
*/


#include <iostream>
#include <limits>

#include <ost/message.hh>
#include <ost/img/alg/stat.hh>

#include "viewer_normalizer.hh"

namespace ost { namespace img { namespace gui {

ViewerNormalizer::ViewerNormalizer():
  NormalizerImpl(),
  RangeHandler(),
  p1_(1.0),
  p2_(0.0),
  p3_(0.0),
  gamma_(1.0),
  invert_sign_(1.0),
  clipmax_(1.0)
{
}
ViewerNormalizer::ViewerNormalizer(Real mininput, Real maxinput, Real minoutput, Real maxoutput, Real gamma, bool invert):
  NormalizerImpl(),
  RangeHandler(mininput,maxinput),
  p1_(pow(maxoutput-minoutput,gamma)*1.0/(maxinput-mininput)),
  p2_(pow(maxoutput-minoutput,gamma)*mininput/(mininput-maxinput)),
  p3_(invert?maxoutput:minoutput),
  gamma_(gamma),
  invert_sign_(invert?-1.0:1.0),
  clipmax_(pow(maxoutput-minoutput,gamma))
{
}
ViewerNormalizer::ViewerNormalizer(const ViewerNormalizer& n):
  NormalizerImpl(n),
  RangeHandler(n),
  p1_(n.p1_),
  p2_(n.p2_),
  p3_(n.p3_),
  gamma_(n.gamma_),
  invert_sign_(n.invert_sign_),
  clipmax_(n.clipmax_)
{
}
ViewerNormalizer::~ViewerNormalizer()
{
}

Real ViewerNormalizer::BackConvert(Real v) const
{
  return Clip((pow(invert_sign_*(v-p3_),gamma_)-p2_)/p1_);
}
Complex ViewerNormalizer::BackConvert(Complex v) const
{
  return Clip((pow(invert_sign_*(v-p3_),gamma_)-p2_)/p1_);
}
Real ViewerNormalizer::Convert(Real v) const
{
  return invert_sign_*pow(ClipMinMax(p1_*v+p2_,0.0,clipmax_),Real(1.0)/gamma_)+p3_;
}
Complex ViewerNormalizer::Convert(Complex v) const
{
  return invert_sign_*pow(ClipMinMax(p1_*v+p2_,0.0,clipmax_),Real(1.0)/gamma_)+p3_;
}

Real ViewerNormalizer::GetGamma()
{
  return gamma_;
}

bool ViewerNormalizer::GetInvert()
{
  return invert_sign_<0.0;
}

}}}  //ns
