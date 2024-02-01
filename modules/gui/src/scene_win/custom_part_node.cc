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

#include "custom_part_node.hh"

#include "entity_part_node.hh"
#include <ost/gui/scene_win/scene_win_model.hh>
#include <ost/gfx/entity.hh>
#include <ost/gfx/gfx_node.hh>

#include <QFont>
namespace ost { namespace gui {

CustomPartNode::CustomPartNode(QString name, gfx::EntityP entity, mol::QueryViewWrapper part, SceneNode* parent):EntityPartNode(name,entity,part,parent){
}

bool CustomPartNode::SetData(int column, const QVariant& value, int role){
  if (column==1 && role == Qt::EditRole) {
    this->SetName(value.toString());
    return true;
  }
  return false;
}

Qt::ItemFlags CustomPartNode::Flags(int column) const{
  if(column==0){
    return Qt::ItemIsSelectable|Qt::ItemIsEnabled;
  }
  else if(column==1){
    return Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable;
  }
  return Qt::NoItemFlags;
}

}}

