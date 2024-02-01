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

#include "label_node.hh"

#include <ost/gfx/entity.hh>
#include <ost/gfx/gfx_node.hh>
#include <ost/gui/scene_win/scene_win_model.hh>
#include <QFont>
namespace ost { namespace gui {

LabelNode::LabelNode(QString name, SceneNode* parent):SceneNode(parent),name_(name){
}

QVariant LabelNode::GetData(int column, int role){
  if(column<0 || column > 2)return QVariant();

  if(column==1) {
    if (role==Qt::DisplayRole) {
      return QVariant(name_);
    } else if(role==Qt::FontRole) {
      QFont f("Helvetica");
      f.setBold(true);
      return QVariant(f);
    }
  }
  return QVariant();
}

int LabelNode::GetColumnCount() const{
  return 2;
}

bool LabelNode::SetData(int column, const QVariant& value, int role){
  if (column==0 && role == Qt::CheckStateRole) {
    if (value.toBool()) {
      this->GetParent()->SetData(0,value,Qt::CheckStateRole);
    } else {
      for(int i = 0; i<this->children().size(); i++){
        this->GetChild(i)->SetData(0,value,Qt::CheckStateRole);
      }
    }
    return true;
  }
  return false;
}

Qt::ItemFlags LabelNode::Flags(int column) const{
  if(column==0 || column==1){
    return Qt::ItemIsSelectable|Qt::ItemIsEnabled;
  }
  return Qt::NoItemFlags;
}

void LabelNode::SetName(QString name){
  name_ = name;
}

QString LabelNode::GetName() const{
  return name_;
}

}}

