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

#include "data_viewer_panel.hh"
#include "draw_overlay.hh"

namespace ost { namespace img { namespace gui {

DrawOverlay::DrawOverlay(ImageHandle imagehandle):
  Overlay("DrawOverlay"),
  imagehandle_(imagehandle),
  pixel_value_(1.0)
{
  menu_=new QMenu;
}

void DrawOverlay::OnMenuEvent(QAction* e) 
{
}

bool DrawOverlay::OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, 
                               const QPoint& lastmouse) 
{
  if(e->buttons()==Qt::LeftButton){
    imagehandle_.SetReal(dvp->WinToPoint(e->x(),e->y()),pixel_value_);
    return true;
  } else if(e->buttons()==Qt::RightButton) {
    imagehandle_.SetReal(dvp->WinToPoint(e->x(),e->y()),0.0);
    return true;
  }
  return false;
}

QMenu* DrawOverlay::GetMenu() 
{
  return menu_;
}


}}}  //ns
