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

#include <ost/message.hh>

#include "data_viewer_panel.hh"
#include "pointlist_overlay.hh"

namespace ost { namespace img { namespace gui {

PointlistOverlay::PointlistOverlay(const String& name):
  PointlistOverlayBase(name),
  pointlist_()
{
  build_();
}

PointlistOverlay::PointlistOverlay(const PointList& pl, const String& name):
  PointlistOverlayBase(name),
  pointlist_()
{
  build_();
  Add(pl);
}

void PointlistOverlay::OnMenuEvent(QAction* e)
{
  if(e==a_clr_) {
    Clear();
  } else {
    PointlistOverlayBase::OnMenuEvent(e);
  }
}

void PointlistOverlay::Add(const Point& pnt, double scale)
{
  pointlist_.push_back(std::pair<Point,double>(pnt,scale));
}

void PointlistOverlay::Add(const PointList& pl, double scale)
{
  for(PointList::const_iterator it=pl.begin();it!=pl.end();++it) {
    Add(*it,scale);
  }
}

void PointlistOverlay::Remove(const Point& pnt)
{
  for(std::vector<std::pair<Point,double> >::iterator it=pointlist_.begin();it!=pointlist_.end();++it){
    if (it->first==pnt) {
      pointlist_.erase(it);
    }
  }
}

void PointlistOverlay::Clear()
{
  pointlist_.clear();
}

void PointlistOverlay::OnDraw(QPainter& pnt, DataViewerPanel* dvp, bool is_active)
{
  std::vector<std::pair<QPoint,double> > qpointlist;
  for(std::vector<std::pair<Point,double> >::iterator it=pointlist_.begin();it!=pointlist_.end();++it){
    qpointlist.push_back(std::pair<QPoint,double> (dvp->FracPointToWinCenter((it->first).ToVec2()),it->second));
  }
  
  DrawVariableSizePointList(pnt,dvp,
		is_active? active_color_ : passive_color_,
		qpointlist);
}

bool PointlistOverlay::OnMouseEvent(QMouseEvent* e, 
                                    DataViewerPanel* dvp,
                                    const QPoint& lastmouse)
{
  if(e->buttons()==Qt::LeftButton && e->modifiers()&Qt::ShiftModifier) {
    // toggle point at mouse position
    Point mpoint = dvp->WinToPoint(e->x(),e->y());
    for(std::vector<std::pair<Point,double> >::iterator it=pointlist_.begin();it!=pointlist_.end();++it){
      if (it->first==mpoint) {
        pointlist_.erase(it);
        return true;
      }
    }
    pointlist_.push_back(std::pair<Point,double>(mpoint,1.0));
    return true;
  }
  return false;
}


void PointlistOverlay::build_()
{
  menu_->addSeparator();
  a_clr_ = menu_->addAction("Clear");
}

}}}  //ns
