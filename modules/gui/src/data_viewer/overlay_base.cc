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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include "overlay_base.hh"

namespace ost { namespace img { namespace gui {


Overlay::Overlay(const String& name):
  name_(name)
{}

void Overlay::SetName(const String& name)
{
  name_ = name;
}

String Overlay::GetName() const
{
  return name_;
}

void Overlay::OnDraw(QPainter& pnt,  DataViewerPanel* dvp, bool is_active)
{
  // default is to do nothing
}

void Overlay::OnMenuEventSlot(QAction* e) 
{
  this->OnMenuEvent(e); // polymorphic call
  emit SettingsChanged();
}

void Overlay::OnMenuEvent(QAction* e) 
{
}

bool Overlay::OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, const QPoint& lastmouse) 
{
  return false;
}

bool Overlay::OnKeyEvent(QKeyEvent* e,  DataViewerPanel* dvp) 
{
  return false;
}

QMenu* Overlay::GetMenu() 
{
  return 0;
}

}}}  //ns
