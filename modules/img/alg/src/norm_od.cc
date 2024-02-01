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
  Author: Andreas Schenk
*/
#include "norm_od.hh"

namespace ost { namespace img { namespace alg {

TransmissionToODNormalizer::TransmissionToODNormalizer():
  NormalizerImpl(),
 scale_(Real(1.0)),
 min_(Real(0.0)),
 min_tr_(Real(1.0e-10))
{
}
TransmissionToODNormalizer::TransmissionToODNormalizer(Real mininput,Real maxinput, Real maxOD):
  NormalizerImpl(),
 scale_(Real(1.0)/(maxinput-mininput)),
 min_(mininput),
 min_tr_(Real(pow(Real(10.0),-maxOD)))
{
}
TransmissionToODNormalizer::TransmissionToODNormalizer(const TransmissionToODNormalizer& n):
  NormalizerImpl(),
 scale_(n.scale_),
 min_(n.min_),
 min_tr_(n.min_tr_)
{
}

Real TransmissionToODNormalizer::BackConvert(Real v) const
{
  return std::max<Real>(Real(pow(Real(10.0),-v))/scale_+min_,min_tr_);
}

Complex TransmissionToODNormalizer::BackConvert(Complex v) const
{
  Real val=Real(pow(Real(10.0),-abs(v)))/scale_+min_;
  return std::max<Real>(val,min_tr_)/abs(v)*v;
}

Real TransmissionToODNormalizer::Convert(Real v) const
{
  Real val=std::max(min_tr_,(v-min_)*scale_);
  return -log10(val);
}

Complex TransmissionToODNormalizer::Convert(Complex v) const
{
  Real val=std::max<Real>(min_tr_,(abs(v)-min_)*scale_);
  return Real(-log10(val))*v/(abs(v));
}



}}} // namespaces
