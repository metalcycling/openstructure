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

/*
 * Niklaus Johner
 */
#ifndef OST_STRUCTURE_ANALYSIS_HH
#define OST_STRUCTURE_ANALYSIS_HH

#include <ost/mol/alg/module_config.hh>

#include <ost/mol/entity_view.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/img/map.hh>

namespace ost { namespace mol { namespace alg {
  geom::Vec3List DLLEXPORT_OST_MOL_ALG GetPosListFromView(const EntityView& view);
  std::vector<Real> DLLEXPORT_OST_MOL_ALG CalculateAgreementWithDensityMap(const geom::Vec3List& vl, img::MapHandle& density_map);
  Real DLLEXPORT_OST_MOL_ALG CalculateAverageAgreementWithDensityMap(const geom::Vec3List& vl, img::MapHandle& density_map);
  void DLLEXPORT_OST_MOL_ALG WrapEntityInPeriodicCell(EntityHandle eh, const geom::Vec3 cell_center, const geom::Vec3 ucell_size, const geom::Vec3 ucell_angles=geom::Vec3(), bool group_res=true, bool follow_bond=true);
  std::vector< std::vector<Real> > DLLEXPORT_OST_MOL_ALG PariwiseDistanceMatrix(const EntityView& view1, const EntityView& view2);
  std::vector< std::vector<Real> > DLLEXPORT_OST_MOL_ALG PariwiseDistanceMatrix(const EntityView& view);
}}}//ns
#endif
