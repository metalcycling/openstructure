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


#include <ost/message.hh>

#include "overlay_manager_gui.hh"
#include "overlay_manager.hh"
#include "overlay.hh"

#include <QScrollArea>
namespace ost { namespace img { namespace gui {

////////////////////////////

OverlayCustomActCheckBox::OverlayCustomActCheckBox(int i, OverlayManager* m):
  QRadioButton(),
  id_(i),
  man_(m)
{
}

void OverlayCustomActCheckBox::OnToggle(bool c)
{
  if(c) man_->ActivateOverlay(id_);
}

////////////////////////////

OverlayCustomVisCheckBox::OverlayCustomVisCheckBox(int i, OverlayManager* m):
  QCheckBox(),
  id_(i),
  man_(m)
{
}

void OverlayCustomVisCheckBox::OnToggle(bool c)
{
  man_->SetOverlayVisibility(id_, c);
}

////////////////////////////

OverlayCustomLockCheckBox::OverlayCustomLockCheckBox(int i, OverlayManager* m):
  QCheckBox(),
  id_(i),
  man_(m)
{
}

void OverlayCustomLockCheckBox::OnToggle(bool c)
{
  man_->SetOverlayLock(id_, c);
}

////////////////////////////

OverlayManagerGUI::OverlayManagerGUI(QWidget* p,OverlayManagerPtr ov_manager):
  QWidget(p),
  ov_manager_(ov_manager)
{
  ov_manager_->Attach(this);

  active_group_ = new QButtonGroup(this);
  active_group_->setExclusive(true);

  QVBoxLayout *vbox = new QVBoxLayout;

  layout_ = new QGridLayout;
  layout_->setHorizontalSpacing(0);
  layout_->setVerticalSpacing(3);
  layout_->addWidget(new QLabel("A"),0,0,Qt::AlignCenter);
  layout_->addWidget(new QLabel("V"),0,1,Qt::AlignCenter);
  layout_->addWidget(new QLabel("L"),0,2,Qt::AlignCenter);
  layout_->addWidget(new QLabel("Overlay"),0,3,Qt::AlignLeft);
  layout_->setColumnStretch(3,1);

  vbox->addLayout(layout_);

  infotext_=new QLabel("");
  QScrollArea* infoarea = new QScrollArea(this);
  infoarea->setWidget(infotext_);
  vbox->addWidget(infoarea);

  setLayout(vbox);
}

void OverlayManagerGUI::OnLockChange(OverlayManager* m, int id, bool lock)
{
  EntryMap::iterator eit = entry_map_.find(id);
  if(eit!=entry_map_.end()) {
    eit->second.l->setChecked(lock);
  }
}

void OverlayManagerGUI::OnVisibilityChange(OverlayManager* m, int id, bool is_visible)
{
  EntryMap::iterator eit = entry_map_.find(id);
  if(eit!=entry_map_.end()) {
    eit->second.v->setChecked(is_visible);
  }
}

void OverlayManagerGUI::OnActivate(OverlayManager* m, int id)
{
  EntryMap::iterator eit = entry_map_.find(id);
  if(eit!=entry_map_.end()) {
    eit->second.a->setChecked(true);
  }
}

void OverlayManagerGUI::OnAddOverlay(OverlayManager* m, int id) 
{
  OverlayPtr ov = m->RetrieveOverlay(id);
  connect(ov.get(),SIGNAL(InfoTextChanged(const QString&)),this,SLOT(SetInfoText(const QString&)));
  OverlayEntry oe;
  oe.row=0.0;
  oe.a=new OverlayCustomActCheckBox(id,m);
  connect(oe.a,SIGNAL(toggled(bool)),oe.a,SLOT(OnToggle(bool)));
  active_group_->addButton(oe.a);
  oe.v=new OverlayCustomVisCheckBox(id,m);
  connect(oe.v,SIGNAL(toggled(bool)),oe.v,SLOT(OnToggle(bool)));
  oe.l=new OverlayCustomLockCheckBox(id,m);
  connect(oe.l,SIGNAL(toggled(bool)),oe.l,SLOT(OnToggle(bool)));
  oe.b=new QPushButton(ov->GetName().c_str(),this);
  oe.b->setFlat(true);
  QMenu* x = ov->GetMenu();
  if(m && x) {
    oe.b->setMenu(x);
    connect(x, SIGNAL(triggered(QAction*)), ov.get(), SLOT(OnMenuEventSlot(QAction*)));
    connect(ov.get(),SIGNAL(SettingsChanged()),this,SIGNAL(SettingsChanged()));
  }

  entry_map_[id]=oe;

  int row = entry_map_.size();
  layout_->addWidget(oe.a,row,0,Qt::AlignCenter);
  layout_->addWidget(oe.v,row,1,Qt::AlignCenter);
  layout_->addWidget(oe.l,row,2,Qt::AlignCenter);
  layout_->addWidget(oe.b,row,3,Qt::AlignLeft);

  if(row==1) {
    oe.a->setChecked(true);
  }
  
  oe.v->setChecked(true);
  oe.l->setChecked(false);
}

void OverlayManagerGUI::OnRemoveOverlay(OverlayManager* m, int id) 
{
  EntryMap::iterator eit = entry_map_.find(id);
  if(eit!=entry_map_.end()) {
    layout_->removeWidget(eit->second.a);
    layout_->removeWidget(eit->second.v);
    layout_->removeWidget(eit->second.l);
    layout_->removeWidget(eit->second.b);
    active_group_->removeButton(eit->second.a);
    eit->second.a->deleteLater();
    eit->second.v->deleteLater();
    eit->second.l->deleteLater();
    eit->second.b->deleteLater();
    entry_map_.erase(eit);
  }
}

void OverlayManagerGUI::OnRelease(OverlayManager* m) 
{
}

void  OverlayManagerGUI::OnEventActivate(QAction* e)
{
}

void  OverlayManagerGUI::OnEventShowHide(QAction* e)
{
}

void  OverlayManagerGUI::OnEventLock(QAction* e)
{
}

void OverlayManagerGUI::OnEventMenuPopup(QAction* e)
{
}

void OverlayManagerGUI::OnEventMenu(QAction* e)
{
}

QButtonGroup* OverlayManagerGUI::GetButtonGroup(void)
{
  return active_group_;
}

void OverlayManagerGUI::UpdateView()
{
}

void OverlayManagerGUI::SetInfoText(const QString& t)
{
  infotext_->setText(t);
  infotext_->adjustSize();
}

}}}  //ns
