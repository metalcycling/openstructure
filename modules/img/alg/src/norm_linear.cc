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

#include "norm_linear.hh"

namespace ost { namespace img { namespace alg {

LinearNormalizer::LinearNormalizer():
  NormalizerImpl(),
 scale_(),
 offset_()
{
}
LinearNormalizer::LinearNormalizer(Real scale, Real offset):
  NormalizerImpl(),
 scale_(scale),
 offset_(offset)
{
}
LinearNormalizer::LinearNormalizer(const LinearNormalizer& n):
  NormalizerImpl(),
 scale_(n.scale_),
 offset_(n.offset_)
{
}


Real LinearNormalizer::BackConvert(Real v) const
{
  return (v-offset_)/scale_;
}

Complex LinearNormalizer::BackConvert(Complex v) const
{
  return (v-offset_)/scale_;
}

Real LinearNormalizer::Convert(Real v) const
{
  return v*scale_+offset_;
}

Complex LinearNormalizer::Convert(Complex v) const
{
  return v*scale_+offset_;
}



LinearRangeNormalizer::LinearRangeNormalizer():
 LinearNormalizer(),
 RangeHandler()
{
}
LinearRangeNormalizer::LinearRangeNormalizer(Real mininput, Real maxinput, Real minoutput, Real maxoutput):
  LinearNormalizer(maxinput!=mininput ? (maxoutput-minoutput)/(maxinput-mininput) : 1.0,
                   maxinput!=mininput ? (minoutput*maxinput-maxoutput*mininput)/(maxinput-mininput) : 0.0),
  RangeHandler(mininput,maxinput)
{
}
LinearRangeNormalizer::LinearRangeNormalizer(const LinearRangeNormalizer& n):
  LinearNormalizer(n),
  RangeHandler(n)
{
}

Real LinearRangeNormalizer::BackConvert(Real v) const
{
  return Clip(LinearNormalizer::BackConvert(v));
}

Complex LinearRangeNormalizer::BackConvert(Complex v) const
{
  return Clip(LinearNormalizer::BackConvert(v));
}

Real LinearRangeNormalizer::Convert(Real v) const
{
  return LinearNormalizer::Convert(Clip(v));
}

Complex LinearRangeNormalizer::Convert(Complex v) const
{
  return LinearNormalizer::Convert(Clip(v));
}

}}} // namespaces

