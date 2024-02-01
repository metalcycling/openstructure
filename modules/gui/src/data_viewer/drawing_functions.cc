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

#include <algorithm>
#include <deque>
#include "drawing_functions.hh"

namespace ost { namespace img { namespace gui {

using namespace geom;

void DrawDiagonalCross(QPainter& pnt, const QPoint& p, int s)
{
  pnt.drawPoint(p);
  for(int i=1;i<s;i++) {
    pnt.drawPoint(p+QPoint(i,i));
    pnt.drawPoint(p+QPoint(-i,-i));
    pnt.drawPoint(p+QPoint(-i,i));
    pnt.drawPoint(p+QPoint(i,-i));
  }
}
void DrawCross(QPainter& pnt, const QPoint& p, int s)
{
  pnt.drawPoint(p);
  for(int i=1;i<s;i++) {
    pnt.drawPoint(p+QPoint(0,i));
    pnt.drawPoint(p+QPoint(0,-i));
    pnt.drawPoint(p+QPoint(-i,0));
    pnt.drawPoint(p+QPoint(i,0));
  }
}

void DrawRectangle(QPainter& pnt, const QPoint& p1, const QPoint& p2) 
{
  pnt.drawRect(QRect(p1,p2));
}

void DrawEllipse(QPainter& pnt, const QPoint& center, Real rx, Real ry, Real angle)
{
  static int numpoints=72;
  QPolygon qpoly;
  for(int j=0; j<numpoints;++j){
    Real ang=static_cast<Real>(j)/static_cast<Real>(numpoints)*2.0*M_PI;
    Real dx=rx*cos(ang);
    Real dy=ry*sin(ang);
    qpoly << QPoint(static_cast<int>(round(center.x()+dx*cos(-angle)+dy*sin(-angle))),static_cast<int>(round(center.y()-dx*sin(-angle)+dy*cos(-angle))));
  }
  pnt.drawPolyline(qpoly);
}

void DrawHyperbola(QPainter& pnt, const QSize& size, const geom::Vec2& center, Real rx, Real ry, Real angle,bool orientation)
{
  bool flagp1=true;
  bool flagp2=true;
  bool flagp3=true;
  bool flagp4=true;
  QPolygon points1,points2;
  Real x=rx;
  Real y=ry;
  QPoint p1,p2,p3,p4;
  do{
    if(orientation){
      y=sqrt((x*x/rx/rx-1.0)*ry*ry);
    }else{
      x=sqrt((y*y/ry/ry-1.0)*rx*rx);
    }
    geom::Vec2 pos=Rotate(Vec2(x,y),-angle)+center;
    p1=QPoint(static_cast<int>(round(pos[0])),static_cast<int>(round(pos[1])));
    pos=Rotate(Vec2(x,-y),-angle)+center;
    p2=QPoint(static_cast<int>(round(pos[0])),static_cast<int>(round(pos[1])));
    pos=Rotate(Vec2(-x,y),-angle)+center;
    p3=QPoint(static_cast<int>(round(pos[0])),static_cast<int>(round(pos[1])));
    pos=Rotate(Vec2(-x,-y),-angle)+center;
    p4=QPoint(static_cast<int>(round(pos[0])),static_cast<int>(round(pos[1])));
    if(orientation){
      if(flagp1){
        points1.push_back(p1);
      }
      if(flagp2){
        points1.push_front(p2);
      }
      if(flagp3){
        points2.push_back(p3);
      }
      if(flagp4){
        points2.push_front(p4);
      }
    }else{
      if(flagp1){
        points1.push_back(p1);
      }
      if(flagp2){
        points2.push_back(p2);
      }
      if(flagp3){
        points1.push_front(p3);
      }
      if(flagp4){
        points2.push_front(p4);
      }
    }
    if(!(p1.x()>=0 && p1.x() <=size.width() && p1.y()>=0 && p1.y() <=size.height())){
      flagp1=false;
    }
    if(!(p2.x()>=0 && p2.x() <=size.width() && p2.y()>=0 && p2.y() <=size.height())){
      flagp2=false;
    }
    if(!(p3.x()>=0 && p3.x() <=size.width() && p3.y()>=0 && p3.y() <=size.height())){
      flagp3=false;
    }
    if(!(p4.x()>=0 && p4.x() <=size.width() && p4.y()>=0 && p4.y() <=size.height())){
      flagp4=false;
    }
    if(orientation){
      x++;
    }else{
      y++;
    }
  }while(flagp1 | flagp2 | flagp3 | flagp4);
  pnt.drawPolyline(points1);
  pnt.drawPolyline(points2);
}

void DrawHalfEllipse(QPainter& pnt, const QPoint& center, Real rx, Real ry, Real angle)
{
  static int numpoints=72;
  QPolygon qpoly;
  for(int j=0; j<numpoints;++j){
    Real ang=j/static_cast<Real>(numpoints)*M_PI+M_PI;
    Real dx=rx*cos(ang);
    Real dy=ry*sin(ang);
    qpoly << QPoint(static_cast<int>(round(center.x()+dx*cos(-angle)+dy*sin(-angle))),static_cast<int>(round(center.y()-dx*sin(-angle)+dy*cos(-angle))));
  }
  pnt.drawPolyline(qpoly);
}

}}}  //ns
