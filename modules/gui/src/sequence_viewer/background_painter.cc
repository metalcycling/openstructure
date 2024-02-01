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
  Author: Stefan Scheuber
 */



#include "background_painter.hh"

#include <QtGui>
namespace ost { namespace gui {

BackgroundPainter::BackgroundPainter(QObject* parent)
    : Painter(parent)
{}

void BackgroundPainter::Paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index){
  painter->save();
  if (index.column()>0 && (index.column()-1)%10 < 5){
    painter->fillRect(option.rect, QColor(240,240,240));

  }
  else{
    painter->fillRect(option.rect, QColor(255,255,255));
  }
  if(index.row()>0 && (index.column())%10 == 0){
    painter->setPen(QPen(QColor(135,135,135)));
    painter->drawLine(option.rect.topRight(),option.rect.bottomRight());
  }
  painter->restore();
}

}}
