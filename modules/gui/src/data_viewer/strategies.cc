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

#include "strategies.hh"

namespace ost { namespace img { namespace gui {

//////////////////////////////
// base

SymbolDrawingStrategy::SymbolDrawingStrategy():
  symbolsize_(10),
  symbolstrength_(1),
  pen_(QPen(QColor(255,0,0),1)),
  brush_(QBrush(QColor(255,0,0),Qt::NoBrush))
{
}

void SymbolDrawingStrategy::Draw(QPainter& pnt, const QPoint& center)
{
  pen_.setWidth(symbolstrength_);
  pnt.setPen(pen_);
  pnt.setBrush(brush_);
  DrawSymbol(pnt,center);
}

//////////////////////////////
// circle

void CircleDrawingStrategy::DrawSymbol(QPainter& pnt, const QPoint& center)
{
  pnt.drawEllipse(center, symbolsize_, symbolsize_);
}

//////////////////////////////
// crosshair
    
CrosshairBaseDrawingStrategy::CrosshairBaseDrawingStrategy():
  SymbolDrawingStrategy(),
  pen2_(QPen())
{
  SetCrosshairColor();	
}

CrosshairBaseDrawingStrategy::CrosshairBaseDrawingStrategy(const SymbolDrawingStrategy& strategy):
  SymbolDrawingStrategy(strategy),
  pen2_(QPen())
{
  SetCrosshairColor();	
}

void CrosshairBaseDrawingStrategy::SetPenColor(const QColor& color)
{
  pen_.setColor(color);
  SetCrosshairColor();
}

void CrosshairBaseDrawingStrategy::SetPen(const QPen& pen)
{
  pen_=pen;
  SetCrosshairColor();
}

void CrosshairBaseDrawingStrategy::SetCrosshairColor()
{
  pen2_.setColor(pen_.color());
}

void CrosshairBaseDrawingStrategy::DrawCrosshair(QPainter& pnt, const QPoint& center)
{
  pnt.drawLine(center-QPoint(symbolsize_,0),center+QPoint(symbolsize_,0));
  pnt.drawLine(center-QPoint(0,symbolsize_),center+QPoint(0,symbolsize_));
}

void CrosshairBaseDrawingStrategy::Draw(QPainter& pnt, const QPoint& center)
{
  pnt.setBrush(brush_);
  if(symbolsize_>=3){
    pen2_.setWidth(1);
    pnt.setPen(pen2_);
    DrawCrosshair(pnt,center);
  }
  pen_.setWidth(symbolstrength_);
  pnt.setPen(pen_);
  DrawSymbol(pnt,center);
}

//////////////////////////////
// crosshair square

void CrosshairSquareDrawingStrategy::DrawSymbol(QPainter& pnt, const QPoint& center)
{
  pnt.drawRect(QRect(center-QPoint(symbolsize_,symbolsize_),
                     center+QPoint(symbolsize_,symbolsize_)));
}

//////////////////////////////
// square

void SquareDrawingStrategy::DrawSymbol(QPainter& pnt, const QPoint& center)
{
  pnt.drawRect(QRect(center-QPoint(symbolsize_,symbolsize_),
                     center+QPoint(symbolsize_,symbolsize_)));
}

//////////////////////////////
// crosshair circle

void CrosshairCircleDrawingStrategy::DrawSymbol(QPainter& pnt, const QPoint& center)
{
  pnt.drawEllipse(center, symbolsize_, symbolsize_);
}

}}}  //ns
