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
  Authors: Ansgar Philippsen
*/


#include <iostream>

#include <ost/message.hh>

#include "null_overlay.hh"

namespace ost { namespace img { namespace gui {

NullOverlay::NullOverlay():
  Overlay("Null")
{}

void NullOverlay::OnDraw(QPainter& pnt,  DataViewerPanel* dvp, bool is_active)
{
  std::cerr << "NullOverlay::OnDraw" << std::endl;
}

bool NullOverlay::OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, 
                               const QPoint& lastmouse)
{
  std::cerr << "NullOverlay::OnMouseEvent" << std::endl;
  return false;
}

bool NullOverlay::OnKeyEvent(QKeyEvent* e,  DataViewerPanel* dvp)
{
  std::cerr << "NullOverlay::OnKeyEvent" << std::endl;
  return false;
}

void NullOverlay::OnMenuEvent(QAction* e)
{
  std::cerr << "NullOverlay::OnMenuEvent: " << e->text().toStdString() << std::endl;
}

QMenu* NullOverlay::GetMenu()
{
  QMenu* menu=new QMenu("Null Menu");
  menu->addAction("One");
  menu->addAction("Two");
  menu->addAction("Three");
  return menu;
}

}}}  //ns
