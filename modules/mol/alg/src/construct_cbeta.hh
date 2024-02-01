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
#ifndef OST_MOL_ALG_CONSTRUCT_CBETA_HH
#define OST_MOL_ALG_CONSTRUCT_CBETA_HH

#include <ost/mol/entity_view.hh>
#include <ost/mol/alg/module_config.hh>

namespace ost { namespace mol { namespace alg {


geom::Vec3 DLLEXPORT_OST_MOL_ALG CBetaPosition(const ResidueHandle& residue, 
                                               Real bond_length=1.5);

geom::Vec3 DLLEXPORT_OST_MOL_ALG CBetaPosition(const geom::Vec3& n_pos, 
                                               const geom::Vec3& ca_pos, 
                                               const geom::Vec3& c_pos, 
                                               Real bond_length=1.5);

void DLLEXPORT_OST_MOL_ALG ConstructCBetas(EntityHandle& entity_handle, 
                                           bool include_gly=false);


}}} // ns

#endif // OST_MOL_ALG_CONSTRUCT_CBETA_HH
