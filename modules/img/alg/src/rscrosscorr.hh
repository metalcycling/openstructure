//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
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
#ifndef IMG_ALG_REAL_SPACE_CROSSCORR_HH
#define IMG_ALG_REAL_SPACE_CROSSCORR_HH

#include <ost/img/image.hh>
#include <ost/img/alg/module_config.hh>
namespace ost { namespace img { namespace alg {

/// \brief Computes cross correlation of images
/// Both input images must be Real Spatial their extent should be the same.
Real DLLEXPORT_IMG_ALG RealSpatialCrossCorrelation(const ConstImageHandle& image1,
                                                   const ConstImageHandle& image2);

/// \brief Computes cross correlation of images
/// the cross correlation is computed only within the given extent.
///
/// Both input images must be Real Spatial and their overlap with the given
/// extent should be the same.
Real DLLEXPORT_IMG_ALG RealSpatialCrossCorrelation(const ConstImageHandle& image1,
                                                   const ConstImageHandle& image2,
                                                   const img::Extent& extent);


}}} // ns

#endif
