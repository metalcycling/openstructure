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
  Authors: Andreas Schenk, Ansgar Philippsen
*/

#ifndef IMG_POINTLIST_OVERLAY_HH
#define IMG_POINTLIST_OVERLAY_HH

#include <ost/img/point_list.hh>

#include <ost/gui/module_config.hh>

#include "pointlist_overlay_base.hh"

namespace ost { namespace img { namespace gui {

class DLLEXPORT_OST_GUI PointlistOverlay: public PointlistOverlayBase
{
public:
  PointlistOverlay(const String& name="Pointlist");
  // initialize with pointlist
  PointlistOverlay(const PointList& pl, const String& name="Pointlist");
  
  virtual void OnMenuEvent(QAction* e);
  virtual void OnDraw(QPainter& pnt,  DataViewerPanel* dvp, bool is_active);
  virtual bool OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, 
                            const QPoint& lastmouse);
  
  void Add(const Point& pnt, double scale=1.0);
  void Add(const PointList& pl, double scale=1.0);
  
  //! remove all points on the given point
  void Remove(const Point& pnt);
  
  // clear all points
  void Clear();
  
private:
  std::vector<std::pair<Point,double> > pointlist_;

  QAction* a_clr_;

  void build_();
};

}}}  //ns


#endif
