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

#include "gfx_scene_node.hh"
#include <ost/gui/scene_win/scene_win_model.hh>
#include <ost/gfx/scene.hh>
#include <ost/gfx/gfx_node.hh>

#include <QFont>
namespace ost { namespace gui {

GfxSceneNode::GfxSceneNode(gfx::GfxNodeP gfx_node, SceneNode* parent):
  SceneNode(parent), gfx_node_(gfx_node) { }

GfxSceneNode::GfxSceneNode(gfx::GfxNodeP gfx_node, QObject* parent):
  SceneNode(parent), gfx_node_(gfx_node) { }

QVariant GfxSceneNode::GetData(int column, int role){
  if (column<0 || column > 2)
    return QVariant();

  if (role==Qt::CheckStateRole && column==0) {
    return QVariant(gfx_node_->IsVisible() ? Qt::Checked : Qt::Unchecked);
  } else if(column==1) {
    if (role==Qt::DisplayRole) {
      return QVariant(QString(gfx_node_->GetName().c_str()));
    } else if(role==Qt::FontRole) {
      QFont f("Helvetica");
      if(gfx_node_->GetType()==0){
        f.setItalic(true);
      }
      return QVariant(f);
    }
  }
  return QVariant();
}

int GfxSceneNode::GetColumnCount() const{
  return 2;
}

bool GfxSceneNode::SetData(int column, const QVariant& value, int role)
{
  if (column==0 && role == Qt::CheckStateRole) {
    if (value.toBool()) {
     this->GetParent()->SetData(0,value,Qt::CheckStateRole);
     gfx_node_->Show();
    } else {
     gfx_node_->Hide();
    }
    return true;
  }
  else if (column==1 && role == Qt::EditRole) {
    gfx_node_->Rename(value.toString().toStdString());
    return true;
  }
  return false;
}

Qt::ItemFlags GfxSceneNode::Flags(int column) const{
  if(column==0){
    return Qt::ItemIsSelectable|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled|Qt::ItemIsDragEnabled;
  }
  else if(column==1){
    return Qt::ItemIsSelectable|Qt::ItemIsEnabled|Qt::ItemIsEditable|Qt::ItemIsDragEnabled;
  }
  return Qt::NoItemFlags;
}

gfx::GfxNodeP GfxSceneNode::GetGfxNode() const{
  return gfx_node_;
}

}}

