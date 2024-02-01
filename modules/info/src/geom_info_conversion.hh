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
#ifndef OST_INFO_GEOM_INFO_CONVERSION_HH
#define OST_INFO_GEOM_INFO_CONVERSION_HH

/*
  A variety of info conversion routines for classes of the geom module.

  Author: Marco Biasini
 */

#include <ost/geom/geom.hh>

#include <ost/info/module_config.hh>
#include <ost/info/info_fw.hh>


namespace ost { namespace info {


/// \name matrix conversion routines
//@{
/// \brief store Mat2 in info group
void DLLEXPORT_OST_INFO Mat2ToInfo(const geom::Mat2& mat, InfoGroup& group);

/// \brief load Mat2 from info group
geom::Mat2 DLLEXPORT_OST_INFO Mat2FromInfo(const InfoGroup& group);

/// \brief store Mat3 in info group
void DLLEXPORT_OST_INFO Mat3ToInfo(const geom::Mat3& mat, InfoGroup& group);

///\brief load Mat3 from info
geom::Mat3 DLLEXPORT_OST_INFO Mat3FromInfo(const InfoGroup& group);

///\brief store Mat4 in info group
void DLLEXPORT_OST_INFO Mat4ToInfo(const geom::Mat4& mat, InfoGroup& group);

///\brief load Mat4 from info
geom::Mat4 DLLEXPORT_OST_INFO Mat4FromInfo(const InfoGroup& group);

/// \brief read transformation from info group
/// \relates Transform
geom::Transform DLLEXPORT_OST_INFO TransformFromInfo(const info::InfoGroup& group);
/// \brief store transformation in info group
/// \relates Transform
void DLLEXPORT_OST_INFO TransformToInfo(const geom::Transform& transform,
                                           info::InfoGroup& group);

//@}
}}

#endif
