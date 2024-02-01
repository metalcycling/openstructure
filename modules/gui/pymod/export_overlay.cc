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
  Author: Ansgar Philippsen
*/

#include <boost/python.hpp>
#include <ost/gui/data_viewer/overlay.hh>
#include <ost/gui/data_viewer/null_overlay.hh>
#include <ost/gui/data_viewer/pointlist_overlay.hh>
#include <ost/gui/data_viewer/draw_overlay.hh>
#include <ost/gui/data_viewer/mask_overlay.hh>

using namespace boost::python;
using namespace ost;
using namespace ost::img;
using namespace ost::img::gui;
namespace {

void plo_add1(PointlistOverlay& plo, const Point& p)
{
  plo.Add(p);
}


void plo_add2(PointlistOverlay& plo, const PointList& p)
{
  plo.Add(p);
}

void plo_add3(PointlistOverlay& plo, const Point& p, double scale)
{
  plo.Add(p,scale);
}

void plo_add4(PointlistOverlay& plo, const PointList& p, double scale)
{
  plo.Add(p,scale);
}

QColor qcolor_from_pyobj(boost::python::object& obj)
{
return QColor(boost::python::extract<int>(obj.attr("red")()),
              boost::python::extract<int>(obj.attr("green")()),
              boost::python::extract<int>(obj.attr("blue")()));
}

void set_active_color(PointlistOverlayBase& plo,boost::python::object& obj)
{
  QColor color=qcolor_from_pyobj(obj);
  plo.SetActiveColor(color);
}
void set_passive_color(PointlistOverlayBase& plo,boost::python::object& obj)
{
  QColor color=qcolor_from_pyobj(obj);
  plo.SetPassiveColor(color);
}

}//ns

void export_overlay()
{
  class_<Overlay, OverlayPtr, boost::noncopyable>("Overlay",no_init)
    .def("GetName",&Overlay::GetName)
    .def("SetName",&Overlay::SetName)
    ;

   class_<NullOverlay,bases<Overlay>,boost::noncopyable>("NullOverlay",init<>())
     ;

  class_<PointlistOverlayBase,bases<Overlay>,boost::noncopyable>("PointlistOverlayBase",init<const String&>())
    .def("SetSymbolSize",&PointlistOverlayBase::SetSymbolSize)
    .def("GetSymbolSize",&PointlistOverlayBase::GetSymbolSize)
    .def("SetSymbolShape",&PointlistOverlayBase::SetSymbolShape)
    .def("GetSymbolShape",&PointlistOverlayBase::GetSymbolShape)
    .def("SetCrosshair",&PointlistOverlayBase::SetCrosshair)
    .def("GetCrosshair",&PointlistOverlayBase::GetCrosshair)
    .def("SetActiveColor",set_active_color)
    .def("SetPassiveColor",set_passive_color)
    ;

  class_<PointlistOverlay,bases<PointlistOverlayBase>,boost::noncopyable>("PointlistOverlay",init<optional<const String&> >())
    .def(init<const PointList&,optional<const String&> >())
    .def("Add",plo_add1)
    .def("Add",plo_add2)
    .def("Add",plo_add3)
    .def("Add",plo_add4)
    .def("Remove",&PointlistOverlay::Remove)
    .def("Clear",&PointlistOverlay::Clear)
    ;

   class_<DrawOverlay,bases<Overlay>,boost::noncopyable>("DrawOverlay",init<ImageHandle>())
     ;

  class_<MaskOverlay,bases<Overlay>,boost::noncopyable>("MaskOverlay",init<>())
    .def(init<const MaskPtr&>())
    .def("SetMask",&MaskOverlay::SetMask)
    .def("GetMask",&MaskOverlay::GetMask)
    .def("ClearMask",&MaskOverlay::ClearMask)
    .def("GetShift",&MaskOverlay::GetShift)
    .def("ClearShift",&MaskOverlay::ClearShift)
    .def("SetShift",&MaskOverlay::SetShift)
    .def("ApplyShiftToMask",&MaskOverlay::ApplyShiftToMask)
    ;

/*  class_<Gauss2DOverlay,bases<Overlay>,boost::noncopyable>("Gauss2DOverlay",init<const alg::ParamsGauss2D&>())
    ;*/
}
