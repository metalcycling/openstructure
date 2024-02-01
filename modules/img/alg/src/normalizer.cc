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

  Author:Andreas Schenk
*/

#include "normalizer.hh"

namespace ost { namespace img {namespace alg{

Normalizer::Normalizer(const NormalizerPtr& ptr):
  ConstModIPAlgorithm("Normalizer"),
  impl_(ptr)
{
}

Normalizer::~Normalizer()
{
}
Real Normalizer::BackConvert(Real v) const
{
  return impl_->BackConvert(v);
}
Complex Normalizer::BackConvert(Complex v) const
{
  return impl_->BackConvert(v);
}
Real Normalizer::Convert(Real v) const
{
  return impl_->Convert(v);
}
Complex Normalizer::Convert(Complex v) const
{
  return impl_->Convert(v);
}
void Normalizer::Visit(ImageHandle& ih) const
{
  impl_->Visit(ih);
}
NormalizerPtr Normalizer::GetImplementation()
{
  return impl_;
}

}}}//ns
