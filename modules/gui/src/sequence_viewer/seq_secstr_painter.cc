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



#include "seq_secstr_painter.hh"

#include <QtGui>
namespace ost { namespace gui {

SeqSecStrPainter::SeqSecStrPainter(QObject* parent)
    : Painter(parent), brush_(QColor(255,255,177,192)), border_pen_(Qt::lightGray), brush_pen_(QColor(255,255,177,192))
{}

void SeqSecStrPainter::Paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index){
  QVariant var = index.data(Qt::UserRole);
  painter->save();
  painter->setPen(QPen(Qt::lightGray));
  const QVarLengthArray<mol::SecStructure>& sec_str = index.data(Qt::UserRole).value<QVarLengthArray<mol::SecStructure> >();
  int column = index.column()-1;
  if(! sec_str.isEmpty() && column < sec_str.size()){
    QSize size = index.data(Qt::UserRole+1).toSize();
    mol::SecStructure sec_element = sec_str[column];

    int std_diff = (size.height()/2);
    int diff = size.height()/1.25;
    int center = option.rect.top()+(option.rect.height()/2);

    if (sec_element.IsCoil()) {
      this->PaintRect(painter,option.rect,center,std_diff);

    } else if (sec_element.IsHelical()) {
      this->PaintRect(painter,option.rect,center,diff);

      if(column - 1 >= 0 && !sec_str[column-1].IsHelical()){
        painter->drawLine(option.rect.left(),center+std_diff,option.rect.left(),center+diff);
        painter->drawLine(option.rect.left(),center-std_diff,option.rect.left(),center-diff);
      }
      if(column + 1 < sec_str.size() && !sec_str[column+1].IsHelical()){
        painter->drawLine(option.rect.right(),center+std_diff,option.rect.right(),center+diff);
        painter->drawLine(option.rect.right(),center-std_diff,option.rect.right(),center-diff);
      }

    } else if (sec_element.IsExtended()) {
      if(column + 1 < sec_str.size() && !sec_str[column+1].IsExtended()){
        this->DrawArrow(painter,option.rect,center,std_diff,diff,size.height());
      }
      else{
        this->PaintRect(painter,option.rect,center,diff);
      }
      if(column - 1 >= 0 && !sec_str[column-1].IsExtended()){
        painter->drawLine(option.rect.left(),center+std_diff,option.rect.left(),center+diff);
        painter->drawLine(option.rect.left(),center-std_diff,option.rect.left(),center-diff);
      }
    }
    if(!(column+1 < sec_str.size())){
      painter->drawLine(option.rect.right(),center+std_diff,option.rect.right(),center-std_diff);
    }
  }
  painter->restore();
}

void SeqSecStrPainter::PaintRect(QPainter* painter, const QRect& rect, int center, int diff)
{
  QPoint top_left(rect.left(),center+diff);
  QPoint bottom_right(rect.right(),center-diff);
  QRect draw_rect = QRect(top_left,bottom_right);
  painter->setBrush(brush_);
  painter->setPen(brush_pen_);
  painter->drawRect(draw_rect);
  painter->setPen(border_pen_);
  painter->drawLine(draw_rect.topLeft(),draw_rect.topRight());
  painter->drawLine(draw_rect.bottomLeft(),draw_rect.bottomRight());
}

void SeqSecStrPainter::DrawArrow(QPainter* painter, const QRect& rect, int center, int std_diff, int diff, int max_diff)
{
  const QPoint points[4] = {
    QPoint(rect.left(), center+max_diff),
    QPoint(rect.right(), center+std_diff),
    QPoint(rect.right(), center-std_diff),
    QPoint(rect.left(), center-max_diff)
  };
  painter->setBrush(brush_);
  painter->setPen(brush_pen_);
  painter->drawPolygon(points, 4);
  painter->setPen(border_pen_);
  painter->drawLine(QPoint(rect.left(),center+diff),points[0]);
  painter->drawLine(QPoint(rect.left(),center-diff),points[3]);
  painter->drawLine(points[0],points[1]);
  painter->drawLine(points[2],points[3]);
}

}}
