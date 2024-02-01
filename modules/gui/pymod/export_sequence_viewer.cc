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
#include <vector>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <ost/gui/sequence_viewer/sequence_viewer.hh>

#include "sip_handler.hh"

using namespace boost::python;
using namespace ost;
using namespace ost::gui;

namespace {

void change_display_mode_a(SequenceViewer* seq_viewer, const QString& mode)
{
  seq_viewer->ChangeDisplayMode(mode);
}

void change_display_mode_b(SequenceViewer* seq_viewer, const gfx::EntityP& entity, const QString& mode)
{
  seq_viewer->ChangeDisplayMode(entity, mode);
}

void change_display_mode_c(SequenceViewer* seq_viewer, const seq::AlignmentHandle& alignment, const QString& mode)
{
  seq_viewer->ChangeDisplayMode(alignment, mode);
}

String get_current_display_mode_a(SequenceViewer* seq_viewer)
{
  return seq_viewer->GetCurrentDisplayMode().toStdString();
}

String get_current_display_mode_b(SequenceViewer* seq_viewer, const gfx::EntityP& entity)
{
  return seq_viewer->GetCurrentDisplayMode(entity).toStdString();
}

String get_current_display_mode_c(SequenceViewer* seq_viewer, const seq::AlignmentHandle& alignment)
{
  return seq_viewer->GetCurrentDisplayMode(alignment).toStdString();
}

std::vector<String> get_display_modes_a(SequenceViewer* seq_viewer)
{
  std::vector<String> modes;
  const QStringList& list = seq_viewer->GetDisplayModes();
  for (int i=0; i<list.size(); i++){
    modes.push_back(list.at(i).toStdString());
  }
  return modes;
}

std::vector<String> get_display_modes_b(SequenceViewer* seq_viewer, const gfx::EntityP& entity)
{
  std::vector<String> modes;
  const QStringList& list = seq_viewer->GetDisplayModes(entity);
  for (int i=0; i<list.size(); i++){
    modes.push_back(list.at(i).toStdString());
  }
  return modes;
}

std::vector<String> get_display_modes_c(SequenceViewer* seq_viewer, const seq::AlignmentHandle& alignment)
{
  std::vector<String> modes;
  const QStringList& list = seq_viewer->GetDisplayModes(alignment);
  for (int i=0; i<list.size(); i++){
    modes.push_back(list.at(i).toStdString());
  }
  return modes;
}

}


void export_SequenceViewer()
{
  class_<SequenceViewer, boost::noncopyable >("SequenceViewer")
    .def(init<optional<bool, bool, const QString&, 
              QWidget*> >())
    .def("Show", &SequenceViewer::show)
    .def("Hide", &SequenceViewer::hide)
    .def("AddEntity", &SequenceViewer::AddEntity)
    .def("RemoveEntity", &SequenceViewer::RemoveEntity)
    .def("SetAlignment", &SequenceViewer::SetAlignment)
    .def("AddAlignment", &SequenceViewer::AddAlignment)
    .def("RemoveAlignment", &SequenceViewer::RemoveAlignment)
    .def("GetDisplayModes", &get_display_modes_a)
    .def("GetDisplayModes", &get_display_modes_b)
    .def("GetDisplayModes", &get_display_modes_c)
    .def("GetCurrentDisplayMode", &get_current_display_mode_a)
    .def("GetCurrentDisplayMode", &get_current_display_mode_b)
    .def("GetCurrentDisplayMode", &get_current_display_mode_c)
    .def("ChangeDisplayMode",&change_display_mode_a)
    .def("ChangeDisplayMode",&change_display_mode_b)
    .def("ChangeDisplayMode",&change_display_mode_c)
    .def("GetQObject",&get_py_qobject<SequenceViewer>)
    .add_property("qobject", &get_py_qobject<SequenceViewer>)
  ;
}

