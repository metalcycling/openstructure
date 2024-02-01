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

#ifndef IMG_GUI_DRAWING_FUNCTIONS_HH_
#define IMG_GUI_DRAWING_FUNCTIONS_HH_

#include <ost/geom/geom.hh>
#include <QPainter>
#include <QPoint>
#include <QSize>


namespace ost { namespace img { namespace gui {

DLLEXPORT void DrawDiagonalCross(QPainter& pnt, const QPoint& p, int s);
DLLEXPORT void DrawCross(QPainter& pnt, const QPoint& p, int s);
DLLEXPORT void DrawRectangle(QPainter& pnt, const QPoint& p1, const QPoint& p2);
DLLEXPORT void DrawEllipse(QPainter& pnt, const QPoint& center, Real rx, Real ry, Real angle);
DLLEXPORT void DrawHalfEllipse(QPainter& pnt, const QPoint& center, Real rx, Real ry, Real angle);
DLLEXPORT void DrawHyperbola(QPainter& pnt, const QSize& s, const geom::Vec2& center, Real rx, Real ry, Real angle,bool orientation=false);

}}}  //ns

#endif

