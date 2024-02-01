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

#ifndef IMG_OV_GUI_OVDIALOG_H
#define IMG_OV_GUI_OVDIALOG_H


#include <vector>
#include <map>
#include "overlay_manager_observer.hh"

#include <QWidget>
#include <QRadioButton>
#include <QCheckBox>
#include <QMenu>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QButtonGroup>
namespace ost { namespace img { namespace gui {

class OverlayCustomActCheckBox: public QRadioButton
{
  Q_OBJECT;
public:
  OverlayCustomActCheckBox(int i, OverlayManager* m);

public slots:
  void OnToggle(bool c);

private:
  int id_;
  OverlayManager* man_;
};

////////////////////////////

class OverlayCustomVisCheckBox: public QCheckBox
{
  Q_OBJECT;
public:
  OverlayCustomVisCheckBox(int i, OverlayManager* m);

public slots:
  void OnToggle(bool c);

private:
  int id_;
  OverlayManager* man_;
};

////////////////////////////

class OverlayCustomLockCheckBox: public QCheckBox
{
  Q_OBJECT;
public:
  OverlayCustomLockCheckBox(int i, OverlayManager* m);

public slots:
  void OnToggle(bool c);

private:
  int id_;
  OverlayManager* man_;
};

////////////////////////////

class OverlayManagerGUI: public QWidget, public OverlayManagerObserver {
  Q_OBJECT;

  struct OverlayEntry {
    int row;
    QRadioButton* a;
    QCheckBox* v;
    QCheckBox* l;
    QPushButton* b;
  };
  
  typedef std::map<int, OverlayEntry> EntryMap;

public:
  OverlayManagerGUI(QWidget* p, OverlayManagerPtr ov_manager);

  // observer interface
  virtual void OnLockChange(OverlayManager* m, int id, bool lock);
  virtual void OnVisibilityChange(OverlayManager* m, int id, bool is_visible);
  virtual void OnActivate(OverlayManager* m, int id) ;
  virtual void OnAddOverlay(OverlayManager* m, int id) ;
  virtual void OnRemoveOverlay(OverlayManager* m, int id) ;
  virtual void OnRelease(OverlayManager* m) ;

  // 
  void OnEventActivate(QAction* e);
  void OnEventShowHide(QAction* e);
  void OnEventLock(QAction* e);
  void OnEventMenuPopup(QAction* e);
  void OnEventMenu(QAction* e);

  QButtonGroup* GetButtonGroup(void);

public slots:
  void SetInfoText(const QString& t);

signals:
  void SettingsChanged();

private:
  OverlayManagerPtr ov_manager_;
  EntryMap entry_map_;
  QButtonGroup* active_group_;
  QGridLayout* layout_;
  QLabel* infotext_;
  
  void UpdateView();
};

}}}  //ns

#endif
