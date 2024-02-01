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
  Authors: Marco Biasini, Niklaus Johner, Ansgar Philippsen
*/

#include <boost/python.hpp>
using namespace boost::python;

#include <ost/mol/coord_frame.hh>
#include <ost/geom/vec3.hh>
#include <ost/mol/entity_handle.hh>

using namespace ost;
using namespace ost::mol;

geom::Vec3 (CoordFrame::*get_atom_pos)(const AtomHandle&) const = &CoordFrame::GetAtomPos;
Real (CoordFrame::*get_dist_atom)(const AtomHandle&, const AtomHandle&) const = &CoordFrame::GetDistanceBetwAtoms;
Real (CoordFrame::*get_angle)(const AtomHandle&, const AtomHandle&, const AtomHandle&) const = &CoordFrame::GetAngle;
Real (CoordFrame::*get_dihedral)(const AtomHandle&, const AtomHandle&, const AtomHandle&, const AtomHandle&) const = &CoordFrame::GetDihedralAngle;
geom::Vec3 (CoordFrame::*get_cm)(const mol::EntityView&) const = &CoordFrame::GetCenterOfMassPos;
Real (CoordFrame::*get_dist_cm)(const mol::EntityView&, const mol::EntityView&) const = &CoordFrame::GetDistanceBetwCenterOfMass;
Real (CoordFrame::*get_min_dist_cm_v)(const mol::EntityView&, const mol::EntityView&) const = &CoordFrame::GetMinDistBetwCenterOfMassAndView;
Real (CoordFrame::*get_rmsd)(const mol::EntityView&, const mol::EntityView&) const = &CoordFrame::GetRMSD;
Real (CoordFrame::*get_min_dist)(const mol::EntityView&, const mol::EntityView&) const = &CoordFrame::GetMinDistance;
Real (CoordFrame::*get_alpha)(const mol::EntityView&) const = &CoordFrame::GetAlphaHelixContent;
geom::Line3 (CoordFrame::*get_odr_line)(const mol::EntityView&) const = &CoordFrame::GetODRLine;
geom::Plane (CoordFrame::*get_odr_plane)(const mol::EntityView&) const = &CoordFrame::GetODRPlane;
boost::python::tuple wrap_FitCylinder(const CoordFrame& f,const mol::EntityView& ev) {
  std::pair<geom::Line3,Real> pair=f.FitCylinder(ev);
  return boost::python::make_tuple<geom::Line3,Real>(pair.first,pair.second);
}  

// TODO: move to geom
geom::Line3 (CoordFrame::*get_odr_line2)() const = &geom::Vec3List::GetODRLine;

CoordFrame create_coord_frame1(const geom::Vec3List& atom_pos)
{
  return CreateCoordFrame(atom_pos);
}

CoordFrame create_coord_frame2(const geom::Vec3List& atom_pos, 
                               const geom::Vec3& cell_size,
                               const geom::Vec3& cell_angles)
{
  return CreateCoordFrame(atom_pos,cell_size,cell_angles);
}

void export_CoordFrame()
{
  // TODO: add ctors or factory functions that take python sequences
  class_<CoordFrame>("CoordFrame",init<size_t, optional<geom::Vec3> >())
    .def("SetCellSize",&CoordFrame::SetCellSize)
    .def("GetCellSize",&CoordFrame::GetCellSize)
    .add_property("cell_size",&CoordFrame::GetCellSize,&CoordFrame::SetCellSize)
    .def("SetCellAngles",&CoordFrame::SetCellAngles)
    .def("GetCellAngles",&CoordFrame::GetCellAngles)
    .add_property("cell_angles",&CoordFrame::GetCellAngles,&CoordFrame::SetCellAngles)
    .def("GetCellVectors",&CoordFrame::GetCellVectors)
    .def("GetAtomPos", get_atom_pos)
    .def("GetDistanceBetwAtoms", get_dist_atom)
    .def("GetAngle", get_angle)
    .def("GetDihedralAngle", get_dihedral)
    .def("GetCenterOfMassPos", get_cm)
    .def("GetDistanceBetwCenterOfMass", get_dist_cm)
    .def("GetMinDistBetwCenterOfMassAndView", get_min_dist_cm_v)
    .def("GetRMSD",get_rmsd)
    .def("GetMinDistance",get_min_dist)
    .def("GetODRPlane",get_odr_plane)
    .def("GetODRLine",get_odr_line)
    .def("GetODRLine",get_odr_line2)
    .def("GetAlphaHelixContent",get_alpha)
    .def("FitCylinder",wrap_FitCylinder)
  ;
  def("CreateCoordFrame",create_coord_frame1);
  def("CreateCoordFrame",create_coord_frame2);
}
