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
#include <ost/gui/gosty_app.hh>
#include <boost/pointer_cast.hpp>


#include <ost/gfx/scene.hh>
#include <ost/gfx/gfx_node.hh>
#include <ost/gfx/entity.hh>
#include <ost/gui/scene_win/scene_win.hh>
#include <ost/gui/scene_win/scene_win_model.hh>
#include "render_modes_node.hh"
#include "render_mode_node.hh"
#include <QFont>

namespace ost { namespace gui {

RenderModesNode::RenderModesNode(gfx::EntityP entity, SceneNode* parent):
  LabelNode("Render Modes",parent),node_(entity) {
  SceneWinModel* model = GostyApp::Instance()->GetSceneWin()->GetModel();
  model->AddNode(parent, this);

  this->Update();

  GostyApp::Instance()->GetSceneWin()->GetModel()->AttachRenderModeObserver(this);
}

void RenderModesNode::RenderModeChanged(){
  this->Update();
}

void RenderModesNode::Update(){
  SceneWinModel* model = GostyApp::Instance()->GetSceneWin()->GetModel();
  gfx::EntityP entity = boost::dynamic_pointer_cast<gfx::Entity>(this->GetGfxNode());
  gfx::RenderModeTypes render_modes =  entity->GetNotEmptyRenderModes();
  for(unsigned int i=0; i<render_modes.size();i++){
    if(!render_types_.contains(render_modes[i])){
      RenderModeNode* node = new RenderModeNode(entity, render_modes[i],this);
      model->AddNode(this, node);
      render_types_.insert(render_modes[i],node);
    }
    else{
      render_types_[render_modes[i]]->Update();
    }
  }
  QSet<gfx::RenderMode::Type> types_to_delete;
  QMap<gfx::RenderMode::Type,RenderModeNode*>::iterator type;
  for (type = render_types_.begin(); type != render_types_.end(); ++type){
    if(find(render_modes.begin(), render_modes.end(), type.key()) == render_modes.end()){
      model->RemoveNode(type.value());
      types_to_delete.insert(type.key());
    }
  }

  QSet<gfx::RenderMode::Type>::iterator to_delete;
  for(to_delete = types_to_delete.begin(); to_delete != types_to_delete.end(); ++to_delete){
    render_types_.remove(*to_delete);
  }

}

gfx::GfxNodeP RenderModesNode::GetGfxNode() const{
  return node_;
}

}}

