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



#include "base_row.hh"

#include <QtGui>
namespace ost { namespace gui {

BaseRow::BaseRow(QObject *parent) : QObject(parent)
{ this->SetFont(font_); }

BaseRow::BaseRow(QFont font, QObject *parent) : QObject(parent)
{
  this->SetFont(font);
}

int BaseRow::GetColumnCount() const
{
  return 0;
}

void BaseRow::Init()
{
  QFontMetrics metrics = QFontMetrics(font_);
  default_font_size_=QSize(metrics.boundingRect('W').width(),metrics.boundingRect('|').height());
  default_cell_size_ = QSize(metrics.boundingRect('W').width()+2,metrics.boundingRect('|').height()*2);
}

void BaseRow::InsertPainter(Painter* painter, int pos)
{
  if(pos == -1 || pos == painter_.size()){
    painter_.append(painter);
  }
  else if(this->IsPainterPosValid(pos)){
    painter_.insert(pos, painter);
  }
}

void BaseRow::RemovePainter(Painter* painter)
{
  painter_.removeAll(painter);
}

Painter* BaseRow::GetPainter(int pos)
{
  if(this->IsPainterPosValid(pos)){
    return painter_[pos];
  }
  return NULL;
}

int BaseRow::GetPainterCount()
{
  return painter_.size();
}

bool BaseRow::IsPainterPosValid(int pos)
{
  if(pos >= 0 && pos < painter_.size()){
    return true;
  }
  return false;
}

const QFont& BaseRow::GetFont() const
{
  return font_;
}

void BaseRow::SetFont(const QFont& font)
{
  font_ = font;
  QFontMetrics metrics = QFontMetrics(font_);
  default_font_size_=QSize(metrics.boundingRect('W').width(),metrics.boundingRect('|').height());
  default_cell_size_ = QSize(metrics.boundingRect('W').width()+2,metrics.boundingRect('|').height()*2);
}

const QSize& BaseRow::GetFontSize() const
{
  return default_font_size_;
}

const QSize& BaseRow::GetCellSize() const
{
  return default_cell_size_;
}

const PainterList& BaseRow::GetPainters() const
{
  return painter_;
}

QVariant BaseRow::GetData(int column, int role) const
{
  return QVariant();
}

bool BaseRow::SetData(int column, const QVariant& value, int role)
{
  return false;
}

Qt::ItemFlags BaseRow::Flags(int column) const
{
  return Qt::NoItemFlags;
}

void BaseRow::DoubleClicked(int column)
{ }

void BaseRow::SetSelection(const QSet<int>& added, const QSet<int>& removed)
{ }

void BaseRow::ZoomIn()
{
  QFont font = this->GetFont();
  int pointsize = font.pointSize();
  if(pointsize < 30){
    font.setPointSize(++pointsize);
  }
  this->SetFont(font);
}

void BaseRow::ZoomOut()
{
  QFont font = this->GetFont();
  int pointsize = font.pointSize();
  if(pointsize > 5){
    font.setPointSize(--pointsize);
  }
  this->SetFont(font);
}

}}
