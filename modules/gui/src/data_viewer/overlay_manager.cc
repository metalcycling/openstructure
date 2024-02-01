//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include <iostream>

#include "overlay_manager.hh"
#include "overlay_base.hh"
#include "data_viewer_panel.hh"

namespace ost { namespace img { namespace gui {

OverlayManager::OverlayManager(QWidget* ref):
  ref_(ref),ov_list_(),next_id_(0),active_ov_(),obs_list_()
{}

OverlayManager::~OverlayManager()
{
  for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it)
    (*it)->OnRelease(this);
}

void OverlayManager::Attach(OverlayManagerObserver* obs)
{
  obs_list_.push_back(obs);
}

void OverlayManager::Detach(OverlayManagerObserver* obs)
{
  obs_list_.remove(obs);
}

int OverlayManager::AddOverlay(const OverlayPtr& ov, bool make_active) 
{
  int this_id = next_id_++;
  ov_list_.push_back(OverlayInfo(ov,this_id,true));
  
  for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it)
    (*it)->OnAddOverlay(this, this_id);

  if(make_active) {
    ActivateOverlay(this_id);
  }

  return this_id;
}

void OverlayManager::DeleteOverlay(const String& name)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
    if(active_ov_==pos->ov){
      active_ov_=ov_list_.begin()->ov;
    }
    ov_list_.erase(pos);
    for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it) {
      (*it)->OnRemoveOverlay(this, pos->id);
    }
  }
}

void OverlayManager::DeleteOverlay(int id)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
    if(active_ov_==pos->ov){
      active_ov_=ov_list_.begin()->ov;
    }
    ov_list_.erase(pos);
    for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it) {
      (*it)->OnRemoveOverlay(this, pos->id);
    }
  }
}

void OverlayManager::DeleteAllOverlays()
{
  for(OverlayList::iterator oit=ov_list_.begin();oit!=ov_list_.end();++oit) {
    for(ObserverList::iterator xit=obs_list_.begin();xit!=obs_list_.end();++xit) {
      (*xit)->OnRemoveOverlay(this, oit->id);
    }
  }
  ov_list_.clear();
}

void OverlayManager::DeleteActiveOverlay()
{
  OverlayList::iterator pos=find_ov(active_ov_);
  if (pos!=ov_list_.end()) {
    for(ObserverList::iterator xit=obs_list_.begin();xit!=obs_list_.end();++xit) {
      (*xit)->OnRemoveOverlay(this, pos->id);
    }
  pos=ov_list_.erase(pos);
  active_ov_.reset();
  }
  
}

void OverlayManager::ActivateOverlay(const String& name)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
    active_ov_=pos->ov;
    for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it) {
      (*it)->OnActivate(this, pos->id);
    }
  }
}

void OverlayManager::ActivateOverlay(int id)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
    active_ov_=pos->ov;
    for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it) {
      (*it)->OnActivate(this, pos->id);
    }
  }
  ref_->update();
}

void OverlayManager::SetOverlayVisibility(const String& name, bool visible)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
    pos->visible=visible;
    for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it) {
      (*it)->OnVisibilityChange(this, pos->id, visible);
    }
  }
  ref_->update();
}

void OverlayManager::SetOverlayVisibility(int id, bool visible)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
    pos->visible=visible;
    for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it) {
      (*it)->OnVisibilityChange(this, pos->id, visible);
    }
  }
  ref_->update();
}

void OverlayManager::SetOverlayLock(const String& name, bool lock)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
    pos->locked=lock;
    for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it) {
      (*it)->OnLockChange(this, pos->id, lock);
    }
  }
  ref_->update();
}

void OverlayManager::SetOverlayLock(int id, bool lock)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
    pos->locked=lock;
    for(ObserverList::iterator it=obs_list_.begin();it!=obs_list_.end();++it) {
      (*it)->OnLockChange(this, pos->id, lock);
    }
  }
  ref_->update();
}

bool OverlayManager::IsVisible(const String& name)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
    return pos->visible;
  }
  return false;
}

bool OverlayManager::IsVisible(int id)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
    return pos->visible;
  }
  return false;
}

bool OverlayManager::IsLocked(const String& name)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
    return pos->locked;
  }
  return false;
}

bool OverlayManager::IsLocked(int id)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
    return pos->locked;
  }
  return false;
}

bool OverlayManager::IsActive(const String& name)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
   return pos->ov==active_ov_;
  }
  return false;
}

bool OverlayManager::IsActive(int id)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
   return pos->ov==active_ov_;
  }
  return false;
}

String OverlayManager::GetOverlayName(int id)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
    return pos->ov->GetName();
  }
  return String();;
}

int OverlayManager::GetOverlayId(const String& name)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
    return pos->id;
  }
  return 0;
}

OverlayPtr OverlayManager::GetActiveOverlay() const 
{
  return active_ov_;
}

OverlayPtr OverlayManager::RetrieveOverlay(const String& name)
{
  OverlayList::iterator pos=find_ov(name);
  if(pos!=ov_list_.end()) {
    return pos->ov;
  }
  return OverlayPtr();
}

OverlayPtr OverlayManager::RetrieveOverlay(int id)
{
  OverlayList::iterator pos=find_ov(id);
  if(pos!=ov_list_.end()) {
    return pos->ov;
  }
  return OverlayPtr();
}

void OverlayManager::OnDraw(QPainter& pnt,  DataViewerPanel* dvp) const
{
  bool active_visible=false;
  for(OverlayList::const_iterator it=ov_list_.begin();it!=ov_list_.end();++it) {
    pnt.save();
    if(it->ov && it->visible) {
      if(it->ov!=active_ov_) {
        it->ov->OnDraw(pnt, dvp, false );
      } else {
        active_visible=it->visible;
      }
    }
    pnt.restore();
  }
  if(active_visible) {
    pnt.save();
      active_ov_->OnDraw(pnt,dvp,true);
    pnt.restore();
  }
}

bool OverlayManager::OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, const QPoint& lastmouse) const
{
  bool ret=false;
  if(active_ov_) {
    ret = active_ov_->OnMouseEvent(e,dvp, lastmouse);
    dvp->update();
  }
  return ret;
}

bool OverlayManager::OnKeyEvent(QKeyEvent* e,  DataViewerPanel* dvp) const
{
  bool ret=false;
  if(active_ov_) {
    ret=active_ov_->OnKeyEvent(e,dvp);
    dvp->update();
  }
  return ret;
}

std::vector<int> OverlayManager::GetIDList() const
{
  std::vector<int> l;
  for(OverlayList::const_iterator it=ov_list_.begin();it!=ov_list_.end();++it) {
    l.push_back(it->id);
  }
  return l;
}

void OverlayManager::RequestRedraw()
{
  if(ref_) {
    ref_->update();
  }
}

// private methods:

OverlayManager::OverlayList::iterator OverlayManager::find_ov(int id) 
{
  for(OverlayList::iterator it=ov_list_.begin();it!=ov_list_.end();++it) {
    if(id==it->id) {
      return it;
      break;
    }
  }
  return ov_list_.end();
}

OverlayManager::OverlayList::iterator OverlayManager::find_ov(const String& name) 
{
  for(OverlayList::iterator it=ov_list_.begin();it!=ov_list_.end();++it) {
    if(name==it->ov->GetName()) {
      return it;
      break;
    }
  }
  return ov_list_.end();
}

OverlayManager::OverlayList::iterator OverlayManager::find_ov(OverlayPtr anoverlay)
{
  for(OverlayList::iterator it=ov_list_.begin();it!=ov_list_.end();++it) {
    if(anoverlay==it->ov) {
      return it;
      break;
    }
  }
  return ov_list_.end();
}

}}}  //ns
