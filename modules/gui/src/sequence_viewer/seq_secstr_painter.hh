//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
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
#ifndef OST_SEQUENCE_VIEWER_SEQ_SECSTR_PAINTER
#define OST_SEQUENCE_VIEWER_SEQ_SECSTR_PAINTER

/*
  Author: Stefan Scheuber
 */


#include "painter.hh"
#include <QObject>
#include <QVarLengthArray>

namespace ost { namespace gui {

class SeqSecStrPainter : public Painter
{
  Q_OBJECT
public:
  SeqSecStrPainter(QObject* parent = 0);
  void Paint(QPainter *painter, const QStyleOptionViewItem &option,
      const QModelIndex &index);
private:
  void PaintRect(QPainter* painter, const QRect& rect, int center, int diff);
  void DrawArrow(QPainter* painter, const QRect& rect, int center, int std_diff, int diff, int max_diff);
  QBrush brush_;
  QPen border_pen_;
  QPen brush_pen_;
};

}}

#endif
