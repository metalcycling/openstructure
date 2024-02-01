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

#ifndef IMG_GUI_OVERLAY_MANAGER_HH
#define IMG_GUI_OVERLAY_MANAGER_HH

#include <vector>
#include <list>



#include <ost/base.hh>
#include <ost/gui/module_config.hh>


#include "overlay_base_fw.hh"
#include "overlay_manager_fw.hh"
#include "overlay_manager_observer.hh"

#include <QWidget>
#include <QPainter>
#include <QEvent>
namespace ost { namespace img { namespace gui {

class DataViewerPanel;

class DLLEXPORT_OST_GUI OverlayManager
{
  struct OverlayInfo {
    OverlayInfo():
      ov(),id(0),visible(true),locked(false) {}

    OverlayInfo(const OverlayPtr& o, int i, bool v):
      ov(o), id(i), visible(v),locked(false) {}

    OverlayPtr ov;
    int id;
    bool visible;
    bool locked;
  };

  typedef std::vector<OverlayInfo> OverlayList;

  typedef std::list<OverlayManagerObserver*> ObserverList;

public:
  OverlayManager(QWidget* ref=0);
  ~OverlayManager();

  //! attach an observer to this manager
  void Attach(OverlayManagerObserver* obs);

  //! detach an observer from this manager
  void Detach(OverlayManagerObserver* obs);

  //! add overlay, returning a unique id
  int AddOverlay(const OverlayPtr& ov, bool make_active=true);

  //! activate overlay by name
  void ActivateOverlay(const String& name);
  //! activate overlay by id
  void ActivateOverlay(int id);

  //! show/hide overlay by name
  void SetOverlayVisibility(const String& name, bool visible);
  //! show/hide overlay by id
  void SetOverlayVisibility(int id, bool visible);

  //! lock/unlock overlay by name
  void SetOverlayLock(const String& name, bool visible);
  //! lock/unlock overlay by id
  void SetOverlayLock(int id, bool visible);

  //! return currently active overlay
  OverlayPtr GetActiveOverlay() const;

  //! retrieve overlay by name
  OverlayPtr RetrieveOverlay(const String& name);
  //! retrieve overlay by id
  OverlayPtr RetrieveOverlay(int id);

  //! return list with all current overlay ids
  std::vector<int> GetIDList() const;
  
  //! delete overlay by name
  void DeleteOverlay(const String& name);
  //! delete overlay by id
  void DeleteOverlay(int id);

  void DeleteAllOverlays();

  void DeleteActiveOverlay();

  void OnDraw(QPainter& pnt,  DataViewerPanel* dvp) const;

  bool OnMouseEvent(QMouseEvent* e,  DataViewerPanel* dvp, const QPoint& lastmouse) const;

  bool OnKeyEvent(QKeyEvent* e,  DataViewerPanel* dvp) const;
  
  bool IsVisible(const String& name);
  bool IsVisible(int id);
  bool IsLocked(const String& name);
  bool IsLocked(int id);
  bool IsActive(const String& name);
  bool IsActive(int id);
  String GetOverlayName(int id);
  int GetOverlayId(const String& name);
  
  void RequestRedraw();

private:
  QWidget* ref_;
  OverlayList ov_list_;
  int next_id_;

  OverlayPtr active_ov_;

  ObserverList obs_list_;

  OverlayList::iterator find_ov(int id);
  OverlayList::iterator find_ov(const String& name);
  OverlayList::iterator find_ov(OverlayPtr anoverlay);
};


}}}  //ns

#endif
