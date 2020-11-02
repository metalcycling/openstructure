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
#ifndef OST_MOL_ALG_SEC_STRUCTURE_SEGMENTS_HH
#define OST_MOL_ALG_SEC_STRUCTURE_SEGMENTS_HH

/*
  Author: Marco Biasini
 */
#include <ost/mol/entity_handle.hh>
#include <ost/mol/alg/module_config.hh>

namespace ost { namespace mol { namespace alg {

/// \brief a consecutive secondary structure element
struct DLLEXPORT_OST_MOL_ALG SecStructureSegment {
  SecStructureSegment(int f, int l, SecStructure s):
    first(f), last(l), ss_type(s) 
  { }
  SecStructureSegment():
    first(0), last(0), ss_type(SecStructure::COIL)
  { }
  bool operator==(const SecStructureSegment& rhs) const
  {
    return first==rhs.first && last==rhs.last && ss_type==rhs.ss_type;
  }
  bool operator!=(const SecStructureSegment& rhs) const
  {
    return !(*this==rhs);
  }

  /// \brief index of first residue part of the segment
  int                first;
  /// \brief index of last residue part of the segment
  int                last;
  /// \brief secondandary structure type of the segment
  SecStructure       ss_type;
};

typedef std::vector<SecStructureSegment> SecStructureSegments;

/// \name extract secondary structure segments
/// 
/// These functions can be used to get a list of consecutive secondary structure 
/// elements, i.e. all helices of a chain or all extended strands.
/// \relates mol::SecStructure
//@{
/// \brief get list of helical segments in chain handle
SecStructureSegments DLLEXPORT_OST_MOL_ALG 
ExtractHelicalSegments(const ChainHandle& chain);
/// \brief get list of helical segments in chain view
SecStructureSegments DLLEXPORT_OST_MOL_ALG 
ExtractHelicalSegments(const ChainView& chain);
/// \brief get list of extended segments in chain handle
SecStructureSegments DLLEXPORT_OST_MOL_ALG 
ExtractExtendedSegments(const ChainHandle& chain);
/// \brief get list of extended segments in chain view
SecStructureSegments DLLEXPORT_OST_MOL_ALG 
ExtractExtendedSegments(const ChainView& chain);
/// \brief get list of all secondary structure segments in chain handle
SecStructureSegments DLLEXPORT_OST_MOL_ALG 
ExtractSecStructureSegments(const ChainHandle& chain);

/// \brief get list of all secondary structure segments in chain view
SecStructureSegments DLLEXPORT_OST_MOL_ALG 
ExtractSecStructureSegments(const ChainView& chain);
//@}

}}}

#endif
