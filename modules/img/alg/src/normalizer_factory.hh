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
  Normalizer factory methods 
   
  Author: Andreas Schenk
*/

#include <ost/img/image.hh>
#include "norm_linear.hh"
#include "norm_log.hh"
#include "norm_od.hh"
#include "normalizer.hh"

#ifndef NORMALIZER_FACTORY_HH_
#define NORMALIZER_FACTORY_HH_

namespace ost { namespace img { namespace alg {
  
DLLEXPORT_IMG_ALG Normalizer CreateNoOpNormalizer();
DLLEXPORT_IMG_ALG Normalizer CreateLinearNormalizer(Real scale, Real offset);
DLLEXPORT_IMG_ALG Normalizer CreateLinearRangeNormalizer(const ConstImageHandle& ih, Real minoutput, Real maxoutput);
DLLEXPORT_IMG_ALG Normalizer CreateLinearRangeNormalizer(Real mininput,Real maxinput, Real minoutput, Real maxoutput);
DLLEXPORT_IMG_ALG Normalizer CreateLogNormalizer(Real scale, Real min, Real offset);
DLLEXPORT_IMG_ALG Normalizer CreateLogRangeNormalizer(const ConstImageHandle& ih, Real minoutput, Real maxoutput);
DLLEXPORT_IMG_ALG Normalizer CreateLogRangeNormalizer(Real mininput,Real maxinput, Real minoutput, Real maxoutput);
DLLEXPORT_IMG_ALG Normalizer CreateTransmissionToODNormalizer(Real mininput,Real maxinput, Real maxOD=10.0);

}}}//ns

#endif /*NORMALIZER_FACTORY_HH_*/
