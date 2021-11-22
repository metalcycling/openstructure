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
  Authors: Stefan Scheuber
*/
#include <boost/pointer_cast.hpp>


#include <ost/dyn_cast.hh>

#include <ost/gfx/scene.hh>

#include <ost/gui/scene_selection.hh>

#include <ost/gfx/map_iso.hh>

#include "custom_part_node.hh"
#include "entity_node.hh"
#include "entity_part_node.hh"
#include "gfx_scene_node.hh"

#include "context_menu.hh"

#include <QAction>
#include <QMenu>
#include <QItemSelection>
#include <QItemSelectionModel>
namespace ost { namespace gui {

ContextMenu::ContextMenu(QTreeView* view, SceneWinModel* model):
  QObject(model),view_(view),model_(model)
{
  QAction* action = new QAction("Center on Object",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(CenterOnObjects()));
  this->AddAction(action, GFX_OBJECT);

  action = new QAction("Deselect",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(Deselect()));
  this->AddAction(action, ENTITY);

  action = new QAction("Delete",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(Delete()));
  this->AddAction(action, GFX_NODE|NOT_SCENE);

  action = new QAction("Show",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(Show()));
  this->AddAction(action, GFX_NODE|NOT_VISIBLE);
  action = new QAction("Hide",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(Hide()));
  this->AddAction(action, GFX_NODE|NOT_HIDDEN);

  action = new QAction("Show",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(MakeVisible()));
  this->AddAction(action, ENTITY_VIEW | VIEWS_SAME_OBJECT);
  action = new QAction("Hide",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(MakeHidden()));
  this->AddAction(action, ENTITY_VIEW | VIEWS_SAME_OBJECT);
  action = new QAction("Show Exclusive",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(ShowExclusive()));
  this->AddAction(action, ENTITY_VIEW | VIEWS_SAME_OBJECT);
  action = new QAction("Hide Exclusive",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(HideExclusive()));
  this->AddAction(action, ENTITY_VIEW | VIEWS_SAME_OBJECT);
  action = new QAction("Select All",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(SelectAllViews()));
  this->AddAction(action, ENTITY_VIEW | VIEWS_SAME_OBJECT);
  action = new QAction("Deselect All",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(DeselectAllViews()));
  this->AddAction(action, ENTITY_VIEW | VIEWS_SAME_OBJECT);
  action = new QAction("Create Custom View",this);
  connect(action, SIGNAL(triggered()), this, SLOT(AddView()));
  this->AddAction(action, ENTITY_VIEW);

  action = new QAction("Delete",this);
  connect(action, SIGNAL(triggered()), this, SLOT(DeleteView()));
  this->AddAction(action, ENTITY_VIEW | CUSTOM_VIEW);

  action = new QAction("Rename",this);
  connect(action, SIGNAL(triggered()), this, SLOT(Rename()));
  this->AddAction(action, GFX_OBJECT | SINGLE);

  action = new QAction("Rename",this);
  connect(action, SIGNAL(triggered()), this, SLOT(Rename()));
  this->AddAction(action, ENTITY_VIEW | SINGLE);

  action = new QAction("View Density Slices",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(ViewDensitySlices()));
  this->AddAction(action, MAP);

  action = new QAction("Show Original Map",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(ShowOriginalMap()));
  this->AddAction(action, MAP | SINGLE | MAP_DOWNSAMPLED);

  action = new QAction("Show Downsampled Map",this);
  connect(action, SIGNAL(triggered()), SceneSelection::Instance(), SLOT(ShowDownsampledMap()));
  this->AddAction(action, MAP | SINGLE | MAP_ORIGINAL | MAP_DSAMPLED_AVAIL);

}

void ContextMenu::ShowMenu(const QPoint& pos)
{
  QModelIndexList indexes = view_->selectionModel()->selection().indexes();
  ContextActionTypes flags;
  flags = ~flags;

  gfx::EntityP view_entity = gfx::EntityP();
  int cnt = 0;
  if(! indexes.empty()){
    for(int i = 0; i < indexes.size(); i++){
      if(indexes[i].column()==0){
        cnt++;
        GfxSceneNode* gfx_scene_node = qobject_cast<GfxSceneNode*>(model_->GetItem(indexes[i]));
        if(gfx_scene_node){
          gfx::GfxNodeP gfx_node = gfx_scene_node->GetGfxNode();
          if(gfx_node->IsVisible()){
            flags &= ~NOT_HIDDEN;
          }
          else{
            flags &= ~NOT_VISIBLE;
          }
          if(gfx::Scene::Instance().GetRootNode() == gfx_node){flags &= ~NOT_SCENE;}
          if(!dyn_cast<gfx::GfxObj> (gfx_node)){flags &= ~GFX_OBJECT;}
          if(!dyn_cast<gfx::Entity> (gfx_node)){flags &= ~ENTITY;}
          if(!dyn_cast<gfx::MapIso>(gfx_node))
          {
            flags &= ~MAP;
          } else {
            gfx::MapIsoP mapisop = dyn_cast<gfx::MapIso> (gfx_node);
            if (mapisop->GetShownMapType() == gfx::ORIGINAL_MAP){
              flags &= ~MAP_DOWNSAMPLED;
            } else {
              flags &= ~MAP_ORIGINAL;
            }
            if (mapisop->IsDownsampledMapAvailable() == false){
              flags &= ~MAP_DSAMPLED_AVAIL;
            }
          }
          if(!dyn_cast<gfx::MapIso> (gfx_node)){
          flags &= ~MAP;

          }
          if(!dyn_cast<gfx::MapIso> (gfx_node)){flags &= ~MAP;}
        }
        else{
          flags &= ~(GFX_NODE | GFX_OBJECT | ENTITY | MAP);
          flags &= ~(NOT_VISIBLE | NOT_HIDDEN);
        }
        EntityPartNode* entity_part_node = qobject_cast<EntityPartNode*>(model_->GetItem(indexes[i]));
        if(!entity_part_node){
          flags &= ~(ENTITY_VIEW | CUSTOM_VIEW);
        }
        else{
          if(view_entity){
            if(!(entity_part_node->GetEntity() == view_entity)){
              flags &= ~VIEWS_SAME_OBJECT;
            }
          }
          else{
            view_entity = entity_part_node->GetEntity();
          }
          CustomPartNode* custom_part_node = qobject_cast<CustomPartNode*>(entity_part_node);
          if(!custom_part_node){
            flags &= ~CUSTOM_VIEW;
          }
        }
      }
    }

    if(cnt==1){
      flags &= ~MULTI;
    }
    else{
      flags &= ~SINGLE;
    }
    QMenu* menu = new QMenu();
    QVectorIterator<std::pair<QAction*, ContextActionTypes> > i(actions_);
    flags ^= NOT_VISIBLE;
    flags ^= NOT_HIDDEN;
    while (i.hasNext()) {
      const std::pair<QAction*, ContextActionTypes>& value = i.next();
      if (value.second == (value.second & flags)){
        menu->addAction(value.first);
      }
    }
    if(! menu->actions().empty()){
      menu->popup(pos);
    }
  }
}

void ContextMenu::AddAction(QAction* action,ContextActionTypes type){
  actions_.push_back(std::make_pair(action, type));
}

void ContextMenu::AddView(){
  gfx::EntityP entity = SceneSelection::Instance()->GetViewEntity();

  EntityNode* ent_node = qobject_cast<EntityNode*>(model_->FindGfxNode(entity));
  if(ent_node)
  {
    mol::EntityView view = SceneSelection::Instance()->GetViewUnion();
    SceneNode* child_node = new CustomPartNode("New View", entity, mol::QueryViewWrapper(view),ent_node->GetCustomViewNode());
    model_->AddNode(ent_node->GetCustomViewNode(),child_node);
    this->Rename(model_->GetIndexOf(child_node,1));
  }
}

void ContextMenu::Rename(QModelIndex index){
  if(index.isValid()){
    view_->expand(index.parent());
    view_->edit(index);
  }
}

void ContextMenu::Rename(){
  QModelIndexList indexes = view_->selectionModel()->selection().indexes();
  if(indexes.size() == 2)
  {
    this->Rename(indexes[1]);
  }
}

void ContextMenu::DeleteView(){
  QModelIndexList indexes = view_->selectionModel()->selection().indexes();
  for(int i = 0; i < indexes.size(); i++){
    if(indexes[i].column()==0){
      model_->RemoveNode(model_->GetItem(indexes[i]));
    }
  }
}

}} // ns

