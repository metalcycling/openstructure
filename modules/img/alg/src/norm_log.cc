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

#include "norm_log.hh"

namespace ost { namespace img { namespace alg {

LogNormalizer::LogNormalizer():
  NormalizerImpl(),
 scale_(1.0),
 min_(0.0),
 offset_(0.0)
{
}
LogNormalizer::LogNormalizer(Real scale,Real min, Real offset):
  NormalizerImpl(),
 scale_(scale),
 min_(min),
 offset_(offset)
{
}
LogNormalizer::LogNormalizer(const LogNormalizer& n):
  NormalizerImpl(),
 scale_(n.scale_),
 min_(n.min_),
 offset_(n.offset_)
{
}

Real LogNormalizer::BackConvert(Real v) const
{
  return exp((v-offset_)/scale_)+min_-Real(1.0);
}

Complex LogNormalizer::BackConvert(Complex v) const
{
  return exp((v-offset_)/scale_)+min_-Real(1.0);
}

Real LogNormalizer::Convert(Real v) const
{
  return log(v-min_+Real(1.0))*scale_+offset_;
}

Complex LogNormalizer::Convert(Complex v) const
{
  return log(v-min_+Real(1.0))*scale_+offset_;
}


LogRangeNormalizer::LogRangeNormalizer():
 LogNormalizer(),
 RangeHandler()
{
}
LogRangeNormalizer::LogRangeNormalizer(Real mininput, Real maxinput, Real minoutput, Real maxoutput):
  LogNormalizer(maxinput!=mininput ? (maxoutput-minoutput)/log(1-mininput+maxinput) : Real(1.0) ,mininput,minoutput),
  RangeHandler(mininput,maxinput)
{
}
LogRangeNormalizer::LogRangeNormalizer(const LogRangeNormalizer& n):
  LogNormalizer(n),
  RangeHandler(n)
{
}

Real LogRangeNormalizer::BackConvert(Real v) const
{
  return Clip(LogNormalizer::BackConvert(v));
}

Complex LogRangeNormalizer::BackConvert(Complex v) const
{
  return Clip(LogNormalizer::BackConvert(v));
}

Real LogRangeNormalizer::Convert(Real v) const
{
  return LogNormalizer::Convert(Clip(v));
}

Complex LogRangeNormalizer::Convert(Complex v) const
{
  return LogNormalizer::Convert(Clip(v));
}


}}} // namespaces
