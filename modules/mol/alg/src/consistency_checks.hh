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
#ifndef OST_MOL_CONSISTENCY_CHECKS_HH
#define OST_MOL_CONSISTENCY_CHECKS_HH

#include <ost/mol/entity_view.hh>
#include <ost/mol/alg/module_config.hh>

namespace ost { namespace mol { namespace alg {
  
/// \brief Checks that residue types with the same ResNum in the two structures match
///
/// Requires a reference structure and a probe structure. The function checks that all the 
/// residues in the reference structure that appear in the probe structure (i.e., that have the 
/// same ResNum) are of the same residue type. Chains are comapred by order, not by chain name 
/// (i.e.: the first chain of the reference will be compared with the first chain of the probe 
/// structure, etc.)
bool DLLEXPORT_OST_MOL_ALG ResidueNamesMatch(const EntityView& probe, 
                                             const EntityView& reference, 
                                             bool log_as_error=false);
  
}}}

#endif


