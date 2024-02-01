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

#include <ost/img/mask_visitor.hh>

#include "data_viewer_panel.hh"
#include "mask_overlay.hh"

namespace ost { namespace img { namespace gui {

MaskOverlay::MaskOverlay():
  Overlay("MaskOverlay"),
  polygons_(),
  menu_(0),
  active_(-1),
  add_mode_(false),
  new_poly_(),
  old_mouse_pos_(),
  shift_()
{
}

MaskOverlay::MaskOverlay(const MaskPtr& m):
  Overlay("MaskOverlay"),
  polygons_(),
  menu_(0),
  active_(-1),
  add_mode_(false),
  new_poly_(),
  old_mouse_pos_(),
  shift_()
{
  SetMask(m);
}

void MaskOverlay::OnDraw(QPainter& pnt, DataViewerPanel* dvp, bool is_active)
{
  QPainter::RenderHints render_hints=pnt.renderHints();
  pnt.setRenderHints(render_hints | QPainter::Antialiasing);
  for (int i=0;i<static_cast<int>(polygons_.size());++i){
    if(is_active){
      if(i==active_){
        pnt.setPen(QPen(QColor(255,255,0),1));
        pnt.setBrush(QColor(255,255,0,30));
      }else{
        pnt.setPen(QPen(QColor(255,0,0),1));
        pnt.setBrush(QColor(255,0,0,30));
      }
    }else{
      pnt.setPen(QPen(QColor(100,100,100),1));
      pnt.setBrush(QColor(100,100,100,30));
    }
    geom::Polygon2 pol=polygons_[i];
    QPolygon qpol;
    for(int j=0;j<static_cast<int>(pol.GetNodeCount());++j){
      qpol << dvp->FracPointToWinCenter(pol.GetNode(j)+shift_);
      pnt.drawEllipse(qpol.back(),4,4);
    }
    pnt.drawPolygon(qpol);
  }

  if(add_mode_) {
    if(is_active){
      pnt.setPen(QPen(QColor(255,0,255),1));
      pnt.setBrush(QColor(255,0,255,30));
    }else{
      pnt.setPen(QPen(QColor(100,100,100),1));
      pnt.setBrush(QColor(100,100,100,30));
    }
    QPolygon qpol;
    for(int j=0;j<static_cast<int>(new_poly_.GetNodeCount());++j){
      qpol << dvp->FracPointToWinCenter(new_poly_.GetNode(j)+shift_);
      pnt.drawEllipse(qpol.back(),4,4);
    }
    pnt.drawPolygon(qpol);
  }
  pnt.setRenderHints(render_hints);
}

bool MaskOverlay::OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, 
                               const QPoint& lastmouse) 
{
  static int active_node_=-1;

  geom::Vec2 mousepos=dvp->WinToFracPoint(e->pos());
  geom::Vec2 diffpos=mousepos-old_mouse_pos_;
  if(e->buttons() & Qt::RightButton && e->button()==Qt::NoButton){
    shift_+=diffpos;
  }
  old_mouse_pos_=mousepos;
  if(e->button() == Qt::LeftButton && e->type()==QEvent::MouseButtonPress){
    if(e->modifiers() & Qt::ShiftModifier) {
      if(!add_mode_) {
        new_poly_=geom::Polygon2();
        add_mode_=true;
      }
      new_poly_.AddNode(mousepos-shift_);
    } else {
      if(active_>=0){
        geom::Polygon2 pol=polygons_[active_];
        for(unsigned int j=0;j<pol.GetNodeCount();++j){
          if(Length(mousepos-(pol.GetNode(j)+shift_))<4){
            active_node_=j;
            return true;
          }
        }
      }
      for(int i=0;i<static_cast<int>(polygons_.size());++i){
        if(i!=active_){
          geom::Polygon2 pol=polygons_[i];
          for(unsigned int j=0;j<pol.GetNodeCount();++j){
            if(Length(mousepos-(pol.GetNode(j)+shift_))<4){
              active_=i;
              active_node_=j;
              return true;
            }
          }
        }
      }
      active_=-1;
      active_node_=-1;
      return false;
    }
  }
  if(e->buttons() & Qt::LeftButton && active_>=0 && active_node_>=0 && ! (e->modifiers() & Qt::ShiftModifier)){
    polygons_[active_].SetNode(active_node_,mousepos-shift_);
  }
  return false;
}

bool MaskOverlay::OnKeyEvent(QKeyEvent* e,  DataViewerPanel* dvp) 
{
  if(e->key()=='R') {
    // reset new polygon
    add_mode_=false;
    return true;
  } else if(e->key()=='A') {
    // accept new polygon
    if(add_mode_) {
      polygons_.push_back(new_poly_);
      add_mode_=false;
      return true;
    }
  }

  return false;
}

QMenu* MaskOverlay::GetMenu() 
{
  return menu_;
}

MaskPtr MaskOverlay::GetMask() const 
{
  MaskPtr nrvo=Mask(geom::Polygon2());
  if (!polygons_.empty()){
    nrvo=Mask(polygons_[0]);
    for(unsigned int i=1;i<polygons_.size();++i){
      nrvo=nrvo | Mask(polygons_[i]);
    }
  }
  return nrvo;
}

namespace {

class MaskOverlayMaskVisitor: public MaskVisitor {
public:
  virtual void VisitPolygonMask(PolygonMask& m) {
    pl_.push_back(geom::Polygon2(m));
  }
  std::vector<geom::Polygon2> GetPolygonList() const {return pl_;}
private:
  std::vector<geom::Polygon2> pl_;
};

} // ns

void MaskOverlay::SetMask(const MaskPtr& m)
{
  MaskOverlayMaskVisitor mv;
  m->Apply(mv);
  polygons_=mv.GetPolygonList();
}

void MaskOverlay::ClearMask()
{
  polygons_.clear();
  active_=-1;
  add_mode_=false;
}

void MaskOverlay::SetShift(geom::Vec2 shift)
{
  shift_=shift;
}

void MaskOverlay::ApplyShiftToMask()
{
  for(std::vector<geom::Polygon2>::iterator it=polygons_.begin();it!=polygons_.end();++it){
    (*it)=(*it)+shift_;
  }
  new_poly_=new_poly_+shift_;
  shift_=geom::Vec2(0.0,0.0);
}


}}}  //ns
