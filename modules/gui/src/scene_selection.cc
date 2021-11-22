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
#include <ost/dyn_cast.hh>

#include <ost/geom/vec3.hh>
#include <ost/geom/vecmat3_op.hh>

#include <ost/mol/view_op.hh>

#include <ost/gfx/scene.hh>
#include <ost/gfx/gfx_object.hh>
#include <ost/gfx/entity.hh>

#include <ost/gfx/map_iso.hh>
#include <ost/gui/main_area.hh>
#include <ost/gui/perspective.hh>
#include <ost/gui/data_viewer/data_viewer.hh>

#include <ost/gui/scene_win/scene_win.hh>

#include "scene_selection.hh"

namespace ost { namespace gui {
SceneSelection* SceneSelection::scene_selection_ = NULL;

SceneSelection::SceneSelection():nodes_(),views_()
{

}

void SceneSelection::SetActiveNodes(gfx::NodePtrList nodes, gfx::EntityP entity, mol::QueryViewWrapperList views){
  nodes_ = nodes;
  view_entity_ = entity;
  views_ = views;
}

gfx::GfxNodeP SceneSelection::GetActiveNode(unsigned int pos) const{
  if(pos < nodes_.size()){
    return nodes_[pos];
  }
  else{
    throw Error("Index out of bounds: There is no active node at the given position");
  }
}

int SceneSelection::GetActiveNodeCount() const{
  return nodes_.size();
}

int SceneSelection::GetActiveViewCount() const{
  return views_.size();
}

mol::EntityView SceneSelection::GetActiveView(unsigned int pos) const{
  if(pos < views_.size()){
    return views_[pos].GetEntityView();
  }
  else{
    throw Error("Index out of bounds: There is no active view at the given position");
  }
}


SceneSelection* SceneSelection::Instance() {
  if (!SceneSelection::scene_selection_) {
    SceneSelection::scene_selection_=new SceneSelection;
  }
  return SceneSelection::scene_selection_;
}

void SceneSelection::CenterOnObjects() {
  geom::Vec3 min=geom::Vec3(std::numeric_limits<Real>::max(),
                  std::numeric_limits<Real>::max(),
                  std::numeric_limits<Real>::max());
  geom::Vec3 max=geom::Vec3(-std::numeric_limits<Real>::max(),
                  -std::numeric_limits<Real>::max(),
                  -std::numeric_limits<Real>::max());
  bool changed=false;
  for(unsigned int i = 0; i < nodes_.size(); i++){
    gfx::GfxNodeP node = nodes_[i];
    if (node) {
      gfx::GfxObj* obj = dynamic_cast<gfx::GfxObj*> (node.get());
      if (obj) {
        geom::AlignedCuboid pair = obj->GetBoundingBox();
        min = geom::Min(min, pair.GetMin());
        max = geom::Max(max, pair.GetMax());
        std::cerr << obj->GetName() << " " << min << " " << max << std::endl;
        changed=true;
      }
    }
  }
  if(changed){
    gfx::Scene::Instance().SetCenter(geom::Vec3(((max.x+min.x)/2),
                                                ((max.y+min.y)/2),
                                                ((max.z+min.z)/2)));
  }
}

void SceneSelection::Delete() {
  QList<gfx::GfxNodeP> selected_objects;
  for(unsigned int i = 0; i < nodes_.size(); i++){
    gfx::GfxNodeP node = nodes_[i];
    selected_objects.append(node);
  }
  for(int i=0; i < selected_objects.size(); i++){
    gfx::Scene::Instance().Remove(selected_objects[i]);
  }
  gfx::Scene::Instance().RequestRedraw();
}

void SceneSelection::ViewDensitySlices() {
  for(unsigned int i = 0; i < nodes_.size(); i++){
    gfx::GfxNodeP node = nodes_[i];
    if (node) {
      gfx::MapIsoP obj = dyn_cast<gfx::MapIso> (node);
      if (obj) {
        // The following is a hack. I need to pass a reference to an ImagHandle
        // that never goes out of scope, so I get a reference from the MapIso using
        // GetMap and pass it to the CreateDataViewer
        GostyApp::Instance()->CreateDataViewer(obj->GetOriginalMap(),QString(obj->GetName().c_str()));
      }
    }
  }
}

void SceneSelection::ShowDownsampledMap()
{
  gfx::MapIsoP obj = dyn_cast<gfx::MapIso> (nodes_[0]);
  if (obj) {
    obj->ShowDownsampledMap();
  }
}


void SceneSelection::ShowOriginalMap()
{
  gfx::MapIsoP obj = dyn_cast<gfx::MapIso> (nodes_[0]);
  if (obj) {
    obj->ShowOriginalMap();
  }
}

void SceneSelection::Deselect(){
  for(unsigned int i = 0; i < nodes_.size(); i++){
    gfx::GfxNodeP node = nodes_[i];
    if (node) {
      gfx::Entity* obj = dynamic_cast<gfx::Entity*> (node.get());
      if (obj) {
        obj->SetSelection(obj->GetView().CreateEmptyView());
      }
    }
  }
}

void SceneSelection::Show() {
  gfx::Scene::Instance().GetRootNode()->Show();
  for(unsigned int i = 0; i < nodes_.size(); i++){
    gfx::GfxNodeP node = nodes_[i];
    if (node) {
      node->Show();
      GostyApp::Instance()->GetSceneWin()->Update();
    }
  }
}

void SceneSelection::Hide() {
  for(unsigned int i = 0; i < nodes_.size(); i++){
    gfx::GfxNodeP node = nodes_[i];
    if (node) {
      node->Hide();            
      GostyApp::Instance()->GetSceneWin()->Update();      
    }
  }
}

mol::EntityView SceneSelection::GetViewUnion() {
  mol::EntityView view;
  if(! views_.empty()){
    view = views_[0].GetEntityView().Copy();
    for(unsigned int i = 1; i < views_.size(); i++){
      view = mol::Union(view,views_[i].GetEntityView());
    }
  }
  return view;
}

void SceneSelection::ShowExclusive(){
  view_entity_->SetVisible(view_entity_->GetView(),false);
  this->MakeVisible();
}

void SceneSelection::HideExclusive(){
  view_entity_->SetVisible(view_entity_->GetView(),true);
  this->MakeHidden();
}

void SceneSelection::MakeVisible(){
  for(unsigned int i= 0; i < views_.size(); i++){
    view_entity_->SetVisible(views_[i].GetEntityView(),true);
  }
}

void SceneSelection::MakeHidden(){
  for(unsigned int i= 0; i < views_.size(); i++){
    view_entity_->SetVisible(views_[i].GetEntityView(),false);
  }
}

void SceneSelection::SelectAllViews(){
  mol::EntityView union_view = this->GetViewUnion();
  if(union_view.IsValid()){
    view_entity_->SetSelection(union_view);
  }
}

void SceneSelection::DeselectAllViews(){
  mol::EntityView sel = view_entity_->GetSelection();
  if(sel.IsValid()){
    sel = mol::Difference(sel, this->GetViewUnion());
    view_entity_->SetSelection(sel);
  }
}



gfx::EntityP SceneSelection::GetViewEntity() const{
  return view_entity_;
}

}}

