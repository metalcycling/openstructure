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

#ifndef IMG_GUI_MASK_OVERLAY_HH_
#define IMG_GUI_MASK_OVERLAY_HH_

#include <ost/geom/geom.hh>
#include <ost/img/image_handle.hh>
#include <ost/img/mask.hh>

#include <ost/gui/module_config.hh>

#include "overlay_base.hh"

namespace ost { namespace img { namespace gui {

// fw decl
class DataViewerPanel;

// overlay base class
class DLLEXPORT_OST_GUI MaskOverlay: public Overlay
{
public:
  MaskOverlay();
  MaskOverlay(const MaskPtr& m);

  virtual bool OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, 
                            const QPoint& lastmouse);
  virtual bool OnKeyEvent(QKeyEvent* e,  DataViewerPanel* dvp);
  virtual void OnDraw(QPainter& pnt,  DataViewerPanel* dvp, bool is_active);
  virtual QMenu* GetMenu();
  
  MaskPtr GetMask() const;
  void SetMask(const MaskPtr& m);
  void ClearMask();
  void ClearShift(){shift_=geom::Vec2();}
  geom::Vec2 GetShift(){return shift_;}
  void SetShift(geom::Vec2 shift);
  void ApplyShiftToMask();

protected:
  std::vector<geom::Polygon2> polygons_;
  QMenu* menu_;
  int active_;
  bool add_mode_;
  geom::Polygon2 new_poly_;
  geom::Vec2 old_mouse_pos_;
  geom::Vec2 shift_;
};

}}}  //ns

#endif
