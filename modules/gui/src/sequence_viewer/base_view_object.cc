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



#include "base_view_object.hh"

#include <QtGui>
namespace ost { namespace gui {

BaseViewObject::BaseViewObject(QObject* parent): QObject(parent)
{ }

void BaseViewObject::InsertRow(int pos, BaseRow* row)
{
  if(pos >= 0 && pos <= rows_.size()){
    rows_.insert(pos,row);
  }
}

void BaseViewObject::RemoveRow(BaseRow* row)
{
  rows_.removeAll(row);
}

BaseRow* BaseViewObject::GetRow(int pos)
{
   if(pos >= 0 && pos < rows_.size()){
     return rows_[pos];
   }
   return NULL;
}

int BaseViewObject::GetRowCount()
{
  return rows_.size();
}

void BaseViewObject::SetSelection(int row, const QSet<int>& added, const QSet<int>& removed)
{
  if(row>=0 && row < rows_.size()){
    rows_[row]->SetSelection(added,removed);
  }
}

QVariant BaseViewObject::GetData(int row, int column, int role)
{
  if(row<0 || row >= rows_.size())return QVariant();

  return rows_[row]->GetData(column,role);
}

int BaseViewObject::GetMaxColumnCount() const
{
  int columns = 0;
  for(int i = 0; i < rows_.size(); i++){
    int col_length = rows_[i]->GetColumnCount();
    if(columns < col_length){
      columns = col_length;
    }
  }
  return columns;
}

bool BaseViewObject::SetData(int row, int column, const QVariant& value, int role)
{
  if(row<0 || row >= rows_.size())return false;

  return rows_[row]->SetData(column, value, role);
}

const QStringList& BaseViewObject::GetDisplayModes()
{
  return display_modes_;
}

const QString& BaseViewObject::GetCurrentDisplayMode()
{
  return current_display_mode_;
}

void BaseViewObject::SetDisplayMode(const QString& mode)
{
  if(display_modes_.contains(mode)){
    this->current_display_mode_ = mode;
  }
}

void BaseViewObject::AddDisplayMode(const QString& mode)
{
  if(!display_modes_.contains(mode)){
    this->display_modes_.append(mode);
  }
}

void BaseViewObject::DoubleClicked(int row, int column)
{
  if(row>=0 || row < rows_.size()){
    rows_[row]->DoubleClicked(column);
  }
}

void BaseViewObject::ZoomIn()
{
  for(int i=0; i< rows_.size(); i++){
    rows_[i]->ZoomIn();
  }
}

void BaseViewObject::ZoomOut()
{
  for(int i=0; i< rows_.size(); i++){
    rows_[i]->ZoomOut();
  }
}

Qt::ItemFlags BaseViewObject::Flags(int row, int column) const
{
  if(row<0 || row >= rows_.size())return Qt::NoItemFlags;

  return rows_[row]->Flags(column);
}

}}
