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
  Authors: Andreas Schenk, Ansgar Philippsen
*/

#ifndef MASK_INFO_CONVERT_HH_
#define MASK_INFO_CONVERT_HH_

#include <ost/info/info_group.hh>
#include "mask_base_fw.hh"

namespace ost { namespace img {

MaskPtr DLLEXPORT InfoToMask(const info::InfoGroup& g);
void DLLEXPORT MaskToInfo(const MaskPtr& mptr, info::InfoGroup& g);

}} //ns

#endif /*MASK_INFO_CONVERT_HH_*/
