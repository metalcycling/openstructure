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

#include "root_node.hh"

#include <ost/gfx/scene.hh>
#include <ost/gfx/gfx_node.hh>
#include <QFont>

namespace ost { namespace gui {

RootNode::RootNode(SceneNode* parent):SceneNode(parent){}

RootNode::RootNode(QObject* parent):SceneNode(parent){}

QVariant RootNode::GetData(int column, int role){
  return QVariant();
}

int RootNode::GetColumnCount() const{
  return 2;
}

bool RootNode::SetData(int column, const QVariant& value, int role){
  return false;
}

Qt::ItemFlags RootNode::Flags(int column) const{
  if(column==0 || column==1){
    return Qt::ItemIsEnabled;
  }
  return Qt::NoItemFlags;
}

}}

