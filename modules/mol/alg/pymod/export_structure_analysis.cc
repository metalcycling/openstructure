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
#include <boost/python.hpp>
using namespace boost::python;

#include <ost/mol/alg/structure_analysis.hh>

using namespace ost;
using namespace ost::mol::alg;

std::vector< std::vector<Real> > (*pair_dist1)(const mol::EntityView&) = PariwiseDistanceMatrix;
std::vector< std::vector<Real> > (*pair_dist2)(const mol::EntityView&,const mol::EntityView&) = PariwiseDistanceMatrix;

void export_StructureAnalysis()
{
  def("GetPosListFromView",&GetPosListFromView, (arg("view")));
  def("CalculateAverageAgreementWithDensityMap",&CalculateAverageAgreementWithDensityMap,(arg("pos_list"),arg("density_map")));
  def("CalculateAgreementWithDensityMap",&CalculateAgreementWithDensityMap,(arg("pos_list"),arg("density_map")));
  def("WrapEntityInPeriodicCell",&WrapEntityInPeriodicCell,(arg("Entity"),arg("cell_center"),arg("ucell_size"),arg("ucell_angles")=geom::Vec3(),arg("group_res")=true,arg("follow_bonds")=false));
  

  def("PariwiseDistanceMatrix",pair_dist2,(arg("EntityView1"),arg("EntityView2")));
  def("PariwiseDistanceMatrix",pair_dist1,(arg("EntityView")));
}
