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

#ifndef IMG_GUI_DRAW_OVERLAY_H
#define IMG_GUI_DRAW_OVERLAY_H

#include <ost/img/image_handle.hh>

#include <ost/gui/module_config.hh>
#include "overlay_base.hh"

namespace ost { namespace img { namespace gui {

// fw decl
class DataViewerPanel;

// overlay base class
class DLLEXPORT_OST_GUI DrawOverlay: public Overlay
{
public:
  //! initialization requires name
  DrawOverlay(ImageHandle imagehandle);

  virtual bool OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, 
                            const QPoint& lastmouse);
  virtual void  OnMenuEvent(QAction* e);
  virtual QMenu* GetMenu();

protected:
  ImageHandle imagehandle_;
  Real pixel_value_;
  QMenu* menu_;
};

}}}  //ns

#endif
