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

#include <ost/mol/query_view_wrapper.hh>
#include <ost/gfx/entity.hh>

#include <ost/gui/scene_win/scene_win_model.hh>
#include "render_mode_node.hh"

#include <QFont>
namespace ost { namespace gui {

RenderModeNode::RenderModeNode(gfx::EntityP entity, gfx::RenderMode::Type render_mode, SceneNode* parent):EntityPartNode(entity->GetRenderModeName(render_mode).c_str(),entity,mol::QueryViewWrapper(entity->GetRenderView(render_mode)),parent),
    entity_(entity),
    render_mode_(render_mode){
}

QVariant RenderModeNode::GetData(int column, int role){
  if(column<0 || column > 2)return QVariant();

  if (role==Qt::CheckStateRole && column==0) {
    return QVariant(entity_->IsRenderModeEnabled(render_mode_) ? Qt::Checked : Qt::Unchecked);
  } else if(column==1) {
    if (role==Qt::DisplayRole) {
      String name = entity_->GetRenderModeName(render_mode_);
      return QVariant(name.c_str());
    } else if(role==Qt::FontRole) {
      QFont f("Helvetica");
      return QVariant(f);
    }
  }
  return QVariant();
}

int RenderModeNode::GetColumnCount() const{
  return 2;
}

bool RenderModeNode::SetData(int column, const QVariant& value, int role){
  if (column==0 && role == Qt::CheckStateRole) {
    if (value.toBool()) {
      this->GetParent()->SetData(0,value,Qt::CheckStateRole);
      entity_->SetEnableRenderMode(render_mode_, true);
    } else {
      entity_->SetEnableRenderMode(render_mode_, false);
    }
    return true;
  }
  return false;
}

Qt::ItemFlags RenderModeNode::Flags(int column) const{
  if(column==0){
    return Qt::ItemIsSelectable|Qt::ItemIsUserCheckable|Qt::ItemIsEnabled;
  }
  else if(column==1){
    return Qt::ItemIsSelectable|Qt::ItemIsEnabled;
  }
  return Qt::NoItemFlags;
}

gfx::RenderMode::Type RenderModeNode::GetRenderMode() const {
  return render_mode_;
}

void RenderModeNode::Update()
{
  this->SetQueryView(mol::QueryViewWrapper(entity_->GetRenderView(this->GetRenderMode())));
}

}}

