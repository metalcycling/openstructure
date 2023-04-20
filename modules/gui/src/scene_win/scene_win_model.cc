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
  Authors: Marco Biasini, Ansgar Philippsen, Stefan Scheuber
*/

#include "scene_win_model.hh"

#include <boost/pointer_cast.hpp>

#include <ost/gfx/scene.hh>
#include <ost/gfx/gfx_object.hh>
#include <ost/gfx/entity_fw.hh>
#include <ost/gfx/entity.hh>

#include <ost/gui/scene_win/root_node.hh>
#include <ost/gui/scene_win/gfx_scene_node.hh>
#include <ost/gui/scene_win/entity_node.hh>

#include <QSize>
#include <QFont>
#include <QDebug>
#include <QDataStream>

namespace ost { namespace gui {

SceneWinModel::SceneWinModel(QObject* parent):
  QAbstractItemModel(parent)
{
  root_node_ = new RootNode(this);
  this->beginInsertRows(QModelIndex(),0,0);
  scene_node_ = new GfxSceneNode(gfx::Scene::Instance().GetRootNode(),root_node_);
  root_node_->AddChild(scene_node_);
  this->endInsertRows();
  gfx::Scene::Instance().AttachObserver(this);
}

SceneWinModel::~SceneWinModel()
{
  gfx::Scene::Instance().DetachObserver(this);
}

gfx::NodePtrList SceneWinModel::GetGfxNodes(QModelIndexList indexes){
  gfx::NodePtrList selected_nodes;
  for(int i = 0; i < indexes.size(); i++){
    if(indexes[i].column()==0){
      GfxSceneNode* gfx_scene_node = qobject_cast<GfxSceneNode*>(GetItem(indexes[i]));
      if(gfx_scene_node){
        gfx::GfxNodeP node = gfx_scene_node->GetGfxNode();
        selected_nodes.push_back(node);
      }
    }
  }
  return selected_nodes;
}

mol::QueryViewWrapperList SceneWinModel::GetQueryViewsList(QModelIndexList indexes){
  mol::QueryViewWrapperList selected_views;
  for(int i = 0; i < indexes.size(); i++){
    if(indexes[i].column()==0){
      EntityPartNode* part_node = qobject_cast<EntityPartNode*>(GetItem(indexes[i]));
      if(part_node){
        mol::QueryViewWrapper node = part_node->GetQueryView();
        selected_views.push_back(node);
      }
    }
  }
  return selected_views;
}

gfx::EntityP SceneWinModel::GetEntityOfViews(QModelIndexList indexes){
  mol::QueryViewWrapperList selected_views;
  gfx::EntityP entity;
  bool entity_set = false;
  for(int i = 0; i < indexes.size(); i++){
    if(indexes[i].column()==0){
      EntityPartNode* part_node = qobject_cast<EntityPartNode*>(GetItem(indexes[i]));
      if(part_node){
        if(entity_set){
          if (entity != part_node->GetEntity()){
            return gfx::EntityP();
          }
        }
        else{
          entity = part_node->GetEntity();
          entity_set = true;
        }
      }
    }
  }
  return entity;
}

void SceneWinModel::Update()
{
  emit this->dataChanged(QModelIndex(),QModelIndex());
}

QModelIndex SceneWinModel::index(int row, int col, 
                                 const QModelIndex& parent) const
{
  if (!hasIndex(row, col, parent)){
    return QModelIndex();
  }
  SceneNode* parent_node = GetItem(parent);

  if (parent_node->children().size() <= row)
    return QModelIndex();
  return createIndex(row, col, parent_node->GetChild(row));
}

QModelIndex SceneWinModel::parent(const QModelIndex& index) const
{
  if(!index.isValid())
    return QModelIndex();
  SceneNode* node = GetItem(index);
  if (!node)
    return QModelIndex();
  assert(node);
  SceneNode* parent_node = node->GetParent();

  if(parent_node==root_node_ || parent_node==NULL) return QModelIndex();

  return createIndex(parent_node->GetRow(), 0, parent_node);
}

int SceneWinModel::rowCount(const QModelIndex& parent) const
{
  SceneNode* parent_node = GetItem(parent);

  return parent_node->children().size();
}

int SceneWinModel::columnCount(const QModelIndex& parent) const
{
  SceneNode* parent_node = GetItem(parent);

  return parent_node->GetColumnCount();
}

QVariant SceneWinModel::data(const QModelIndex& index, int role) const
{
  if (!index.isValid()) return QVariant();

  SceneNode* node = GetItem(index);
  if(!node) return QVariant();
  QVariant data = node->GetData(index.column(),role);
  return data;
}

QVariant SceneWinModel::headerData(int section, Qt::Orientation orientation,
                                    int role) const
{
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole) {
    if(section==0) return QVariant("Visible");
    else if(section==1) return QVariant("Node");
  }
  return QVariant();
}

Qt::ItemFlags SceneWinModel::flags(const QModelIndex& index) const
{
  if (!index.isValid())
    return Qt::NoItemFlags;
  SceneNode* node = reinterpret_cast<SceneNode*>(index.internalPointer());
  if(node){
    return node->Flags(index.column());
  }

  return QAbstractItemModel::flags(index);
}

bool SceneWinModel::setData(const QModelIndex& index,
                             const QVariant& value, int role)
{
  if (!index.isValid()) {
    return false;
  }

  SceneNode* node = reinterpret_cast<SceneNode*>(index.internalPointer());
  if(node){
    bool ret = node->SetData(index.column(),value,role);
    if(ret && role==Qt::CheckStateRole){
      this->Update();
    }
    return ret;

  }
  return false;
}

QStringList SceneWinModel::mimeTypes() const
{
    QStringList types;
    types << "text/plain";
    return types;
}

Qt::DropActions SceneWinModel::supportedDragActions() const
{
    return Qt::MoveAction;
}

QMimeData* SceneWinModel::mimeData(const QModelIndexList &indexes) const
{
  QMimeData *mimeData = new QMimeData();
  QByteArray encoded_data;

  QDataStream stream(&encoded_data, QIODevice::WriteOnly);


  foreach (QModelIndex index, indexes) {
    if (index.isValid() && index.column()==1) {
      QString text = "scene['"+data(index, Qt::DisplayRole).toString()+"']";
      encoded_data.append(text.toUtf8());
    }
  }
  mimeData->setData("text/plain", encoded_data);
  return mimeData;
}

void SceneWinModel::NodeAdded(const gfx::GfxNodeP& node)
{
  gfx::EntityP e=boost::dynamic_pointer_cast<gfx::Entity>(node);
  SceneNode* parent_node=this->FindGfxNode(node->GetParent());
  if (parent_node==NULL) {
    parent_node=scene_node_;
  }
  if (e) {
    EntityNode* ent_node=new EntityNode(e,parent_node);
    node_map_[node.get()]=ent_node;
  }
  else{
    SceneNode* scene_node = new GfxSceneNode(node, parent_node);
    node_map_[node.get()]=scene_node;
    this->AddNode(parent_node, scene_node);
    for (gfx::GfxNodeVector::const_iterator i=node->GetChildren().begin(),
         e=node->GetChildren().end(); i!=e; ++i) {
      this->NodeAdded(*i);
    }
  }
}

void SceneWinModel::NodeRemoved(const gfx::GfxNodeP& node)
{
  SceneNode* scene_node = this->FindGfxNode(node);
  if(scene_node && scene_node->GetParent()){
    node_map_.remove(node.get());
    this->RemoveNode(scene_node);
  }
}

void SceneWinModel::SelectionChanged(const gfx::GfxObjP& obj, 
                                     const mol::EntityView& sel)
{
}

void SceneWinModel::RenderModeChanged(const gfx::GfxNodeP& node)
{
  if(render_observers_.contains(node.get())){
    render_observers_[node.get()]->RenderModeChanged();
  }
}

SceneNode* SceneWinModel::GetItem(const QModelIndex &index) const
{
  if (index.isValid()) {
    SceneNode *node = reinterpret_cast<SceneNode*>(index.internalPointer());
    if (node) return node;
  }
  return root_node_;
}

bool SceneWinModel::AddNode(SceneNode* parent, SceneNode* child)
{
  QModelIndex parent_index=this->GetIndexOf(parent);
  if(parent_index.isValid()){
    int row = parent->children().size();
    this->beginInsertRows(parent_index,row,row);
    parent->AddChild(child);
    this->endInsertRows();
    return true;
  }
  return false;
}

bool SceneWinModel::RemoveNode(SceneNode* node)
{
  QModelIndex index=this->GetIndexOf(node);
  SceneNode* parent=node->GetParent();
  if(parent && index.isValid()){
    QModelIndex parent_index =this->GetIndexOf(parent);
    int row = node->GetRow();
    this->beginRemoveRows(parent_index,row,row);
    delete node;
    this->endRemoveRows();
    return true;
  }
  return false;
}


void SceneWinModel::AttachRenderModeObserver(RenderModesNode* node)
{
  render_observers_.insert(node->GetGfxNode().get(),node);
}

void SceneWinModel::DetachRenderModeObserver(RenderModesNode* node)
{
  if(render_observers_.contains(node->GetGfxNode().get())){
    render_observers_.remove(node->GetGfxNode().get());
  }
}

QModelIndex SceneWinModel::GetIndexOf(SceneNode* node, int column) 
{
  std::vector<SceneNode*> trace;
  SceneNode* parent=node->GetParent();
  while (parent) {
    trace.push_back(parent);
    parent=parent->GetParent();
    
  }
  if (!trace.empty()) {
    trace.pop_back();
  }
  QModelIndex parent_index=this->index(0, 0, QModelIndex());
  for (std::vector<SceneNode*>::reverse_iterator i=trace.rbegin(), 
       e=trace.rend(); i!=e; ++i) {
    SceneNode* cur=*i;
    parent_index=this->GetIndex(cur, parent_index, 0);
  }
  return this->GetIndex(node, parent_index, column);
}

QModelIndex SceneWinModel::GetIndex(SceneNode* node, QModelIndex parent, 
                                    int column) {
  if(parent.isValid()){
    SceneNode* parent_node =reinterpret_cast<SceneNode*>(parent.internalPointer());
    if (parent_node == node) return parent;
    for (int i=parent_node->children().size()-1; i>=0; --i) {
      SceneNode* child = parent_node->GetChild(i);
      if(child == node){
        if(column<child->GetColumnCount())
          return index(i,column,parent);
        else
          return QModelIndex();
      }
      QModelIndex found = GetIndex(node, index(i,0,parent),column);
      if (found.isValid()) {
        return found;
      }
    }
  }
  return QModelIndex();
}

SceneNode* SceneWinModel::FindGfxNode(gfx::GfxNodeP node) 
{
  QMap<gfx::GfxNode*, SceneNode*>::iterator i=node_map_.find(node.get());
  if (i!=node_map_.end()) {
    return i.value();
  }
  return NULL;
}

}} // ns

