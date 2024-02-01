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
/*
  Authors: Marco Biasini, Andreas Schenk
 */
#include "python_shell_widget.hh"
#include "gutter.hh"

#include <QPainter>
#include <QTextBlock>

namespace ost { namespace gui {

Gutter::Gutter(PythonShellWidget* shell):
  QWidget(shell)
{  
}

void Gutter::Update(QRect area, int)
{
  update(0, area.y(), width(), area.height());
}

QSize Gutter::sizeHint() const
{
  return QSize(GUTTER_WIDTH, 0);
}

void Gutter::paintEvent(QPaintEvent* paint_event)
{
  PythonShellWidget* py_shell=dynamic_cast<PythonShellWidget*>(parent());
  GutterBlockList blocks=py_shell->GetGutterBlocks(paint_event->rect());
  QPainter painter(this);
  painter.setPen(QPen(Qt::black));
  painter.fillRect(paint_event->rect(), QBrush(Qt::lightGray));
  for (GutterBlockList::const_iterator it=blocks.begin();
       it!=blocks.end(); ++it) {
    if (it->type & BLOCKTYPE_ACTIVE) {
      QColor green(Qt::green);
      painter.setPen(QPen(green.darker(), 3));
    }else  if (it->type & BLOCKTYPE_BLOCKEDIT) {
        painter.setPen(QPen(Qt::yellow, 3));
    }else  if (it->type & BLOCKTYPE_CODE) {
        QColor green(Qt::green);
        painter.setPen(QPen(green.darker(), 1));
    } else if (it->type & BLOCKTYPE_OUTPUT) {
      painter.setPen(QPen(Qt::blue, 1));
    } else if (it->type & BLOCKTYPE_ERROR) {
      painter.setPen(QPen(Qt::red, 1));
    } else {
      painter.setPen(QPen(Qt::black, 1));
    }
    painter.drawLine(QLineF(Gutter::GUTTER_WIDTH/2, it->start,
                            Gutter::GUTTER_WIDTH/2, it->end));      
  }
}


}}//ns
