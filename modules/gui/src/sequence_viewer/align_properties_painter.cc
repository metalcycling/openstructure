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



#include "align_properties_painter.hh"

#include <QtGui>
namespace ost { namespace gui {

namespace {

QMap<QString,QColor> GetColorMap(){
  QMap<QString,QColor> map;
  map["G"]=QColor(175,175,175,100);
  map["A"]=QColor(175,175,175,100);
  map["V"]=QColor(175,175,175,100);
  map["L"]=QColor(175,175,175,100);
  map["I"]=QColor(175,175,175,100);
  map["F"]=QColor(255,165,100);
  map["Y"]=QColor(255,165,100);
  map["W"]=QColor(255,165,100);
  map["C"]=QColor(255,255,0,100);
  map["M"]=QColor(255,255,0,100);
  map["S"]=QColor(0,255,0,100);
  map["T"]=QColor(0,255,0,100);
  map["K"]=QColor(255,0,0,100);
  map["R"]=QColor(255,0,0,100);
  map["H"]=QColor(255,0,0,100);
  map["D"]=QColor(0,0,255,100);
  map["E"]=QColor(0,0,255,100);
  map["N"]=QColor(0,0,255,100);
  map["Q"]=QColor(0,0,255,100);
  map["P"]=QColor(0,255,255,100);
  return map;
}

}

QMap<QString,QColor> AlignPropertiesPainter::color_map_ = GetColorMap();


AlignPropertiesPainter::AlignPropertiesPainter(QObject* parent)
    : Painter(parent)
{}

void AlignPropertiesPainter::Paint(QPainter* painter, const QStyleOptionViewItem& option, const QModelIndex& index){
  painter->save();
  if (index.column()>0){
    QString text = index.data(Qt::DisplayRole).toString();
    if(color_map_.contains(text)){
      painter->fillRect(option.rect, color_map_[text]);
    }
  }
  painter->restore();
}

}}
