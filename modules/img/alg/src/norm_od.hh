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
#ifndef NORM_OD_HH_
#define NORM_OD_HH_


#include <ost/img/normalizer_impl.hh>
#include <ost/img/alg/module_config.hh>

namespace ost { namespace img { namespace alg {

class DLLEXPORT_IMG_ALG TransmissionToODNormalizer: public NormalizerImpl {
public:
  // ctor
  TransmissionToODNormalizer();
  TransmissionToODNormalizer(Real mininput,Real maxinput, Real maxOD=10.0);
  // copy ctor
  TransmissionToODNormalizer(const TransmissionToODNormalizer& n);
  virtual ~TransmissionToODNormalizer() {}

  // normalizer interface
  virtual Real BackConvert(Real v) const;
  virtual Complex BackConvert(Complex v) const;
  virtual Real Convert(Real v) const;
  virtual Complex Convert(Complex v) const;
  virtual NormalizerImpl* Clone() const {return new TransmissionToODNormalizer((*this));}

private:
  Real scale_;
  Real min_;
  Real min_tr_;
};


}}} // namespaces

#endif /*NORM_OD_HH_*/
