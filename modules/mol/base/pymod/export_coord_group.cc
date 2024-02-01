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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

#include <ost/mol/coord_group.hh>
#include <ost/mol/entity_handle.hh>

using namespace ost;
using namespace ost::mol;

namespace {
  geom::Vec3 cg_getitem(const CoordGroupHandle& cg, tuple i)
  {
    return (*cg.GetFrame(extract<int>(i[0])))[extract<int>(i[1])];
  }
  void cg_setitem(CoordGroupHandle& cg, tuple i, const geom::Vec3& v)
  {
    (*cg.GetFrame(extract<int>(i[0])))[extract<int>(i[1])]=v;
  }

  void (CoordGroupHandle::*capture1)()=&CoordGroupHandle::Capture;
  void (CoordGroupHandle::*capture2)(uint)=&CoordGroupHandle::Capture;
  void (CoordGroupHandle::*add_frame1)(const geom::Vec3List&)=&CoordGroupHandle::AddFrame;
  void (CoordGroupHandle::*add_frame2)(const geom::Vec3List&,const geom::Vec3&,const geom::Vec3&)=&CoordGroupHandle::AddFrame;
}

void export_CoordGroup()
{
  class_<CoordGroupHandle>("CoordGroupHandle",no_init)
    .def("IsValid",&CoordGroupHandle::IsValid)
    .def("GetEntity",&CoordGroupHandle::GetEntity)
    .add_property("entity",&CoordGroupHandle::GetEntity)
    .def("GetAtomCount",&CoordGroupHandle::GetAtomCount)
    .add_property("atom_count",&CoordGroupHandle::GetAtomCount)
    .def("GetFrame",&CoordGroupHandle::GetFrame2)
    .def("AddFrames", &CoordGroupHandle::AddFrames)
    .def("AddFrame", add_frame1)
    .def("AddFrame", add_frame2)
    .def("GetFrameCount",&CoordGroupHandle::GetFrameCount)
    .add_property("frame_count",&CoordGroupHandle::GetFrameCount)
    .def("GetFramePositions",&CoordGroupHandle::GetFramePositions)
    .def("SetFramePositions",&CoordGroupHandle::SetFramePositions)
    .def("SetFrameCellSize",&CoordGroupHandle::SetFrameCellSize)
    .def("SetAtomPos",&CoordGroupHandle::SetAtomPos)
    .def("GetAtomPos",&CoordGroupHandle::GetAtomPos)
    .def("CopyFrame",&CoordGroupHandle::CopyFrame)
    .def("IsValid", &CoordGroupHandle::IsValid)
    .def("Capture", capture1)
    .def("Capture", capture2)
    .def("CaptureInto",&CoordGroupHandle::CaptureInto)
    .def("GetAtomList",&CoordGroupHandle::GetAtomList)
    .add_property("atoms", &CoordGroupHandle::GetAtomList)
    .def("GetDelta",&CoordGroupHandle::GetDelta)
    .def("SetDelta",&CoordGroupHandle::SetDelta)
    .add_property("delta",&CoordGroupHandle::GetDelta,&CoordGroupHandle::SetDelta)
    .def("GetStartTime",&CoordGroupHandle::GetStartTime)
    .def("SetStartTime",&CoordGroupHandle::SetStartTime)
    .add_property("start_time",&CoordGroupHandle::GetStartTime,&CoordGroupHandle::SetStartTime)
    .def("__getitem__",cg_getitem)
    .def("__setitem__",cg_setitem)
    .def("Filter", &CoordGroupHandle::Filter, (arg("selected"),arg("first")=0,arg("last")=-1,arg("stride")=1))
  ;

  def("CreateCoordGroup",CreateCoordGroup);
}
