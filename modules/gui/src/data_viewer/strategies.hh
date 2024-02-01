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

#ifndef IMG_GUI_STRATEGIES_HH_
#define IMG_GUI_STRATEGIES_HH_


#include <ost/geom/geom.hh>

#include "strategies_fw.hh"

#include <QPainter>
#include <QPoint>
#include <QColor>
#include <QPen>
#include <QBrush>
namespace ost { namespace img { namespace gui {


class SymbolDrawingStrategy
{
public:
  SymbolDrawingStrategy();
  virtual ~SymbolDrawingStrategy() {};

  virtual void Draw(QPainter& pnt, const QPoint& center);

  virtual void SetPenColor(const QColor& color){pen_.setColor(color);}
  virtual QColor GetPenColor() const {return pen_.color();}
  void SetPen(const QPen& pen) {pen_=pen;}
  virtual QPen GetPen() const {return pen_;}
  virtual void SetBrush(const QBrush& brush) {brush_=brush;}
  virtual QBrush GetBrush() const {return brush_;}
  virtual unsigned int GetSymbolSize() const {return symbolsize_;}
  virtual void SetSymbolSize(unsigned int symbolsize) {symbolsize_=symbolsize;}
  virtual unsigned int GetSymbolStrength() const {return symbolstrength_;}
  virtual void SetSymbolStrength(unsigned int s) {symbolstrength_=s;}
  virtual bool HasCrosshair(){return false;}
  virtual String GetShape()=0;
protected:
  virtual void DrawSymbol(QPainter& pnt, const QPoint& center)=0;
  unsigned int symbolsize_;
  unsigned int symbolstrength_;
  QPen pen_;
  QBrush brush_;
};


class SquareDrawingStrategy: public SymbolDrawingStrategy
{
public:
  SquareDrawingStrategy():SymbolDrawingStrategy()
  {}
  SquareDrawingStrategy(const SymbolDrawingStrategy& strategy):
    SymbolDrawingStrategy(strategy)
  {}
  virtual String GetShape(){return "Square";}
protected:
  virtual void DrawSymbol(QPainter& pnt, const QPoint& center);
};

class CrosshairBaseDrawingStrategy: public SymbolDrawingStrategy
{
public:
  CrosshairBaseDrawingStrategy();
  CrosshairBaseDrawingStrategy(const SymbolDrawingStrategy& strategy);
  void SetPen(const QPen& pen);
  void SetPenColor(const QColor& color);
  virtual void Draw(QPainter& pnt, const QPoint& center);
  virtual bool HasCrosshair(){return true;}
protected:
  virtual void DrawCrosshair(QPainter& pnt, const QPoint& center);
  void SetCrosshairColor();
  QPen pen2_;
};

class CrosshairCircleDrawingStrategy: public CrosshairBaseDrawingStrategy
{
public:
  CrosshairCircleDrawingStrategy():
    CrosshairBaseDrawingStrategy()
  {}
  CrosshairCircleDrawingStrategy(const SymbolDrawingStrategy& strategy):
    CrosshairBaseDrawingStrategy(strategy){}
  virtual String GetShape(){return "Circle";}
protected:
  virtual void DrawSymbol(QPainter& pnt, const QPoint& center);
};

class CrosshairSquareDrawingStrategy: public CrosshairBaseDrawingStrategy
{
public:
  CrosshairSquareDrawingStrategy():
    CrosshairBaseDrawingStrategy()
  {}
  CrosshairSquareDrawingStrategy(const SymbolDrawingStrategy& strategy):
    CrosshairBaseDrawingStrategy(strategy)
  {}
  virtual String GetShape(){return "Square";}
protected:
  virtual void DrawSymbol(QPainter& pnt, const QPoint& center);
};

class CircleDrawingStrategy: public SymbolDrawingStrategy
{
public:
  CircleDrawingStrategy():
    SymbolDrawingStrategy()
  {}
  CircleDrawingStrategy(const SymbolDrawingStrategy& strategy):
    SymbolDrawingStrategy(strategy)
  {}
  virtual String GetShape(){return "Circle";}
protected:
  virtual void DrawSymbol(QPainter& pnt, const QPoint& center);
};

}}}  //ns

#endif /*STRATEGIES_HH_*/
