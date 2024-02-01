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
#include <cassert>

#include <ost/gui/widget_geom_handler.hh>
#include "main_area.hh"

#include <QMdiSubWindow>
#include <QSettings>
/*
  Author: Marco Biasini
 */

namespace ost { namespace gui {

namespace {

class PersistentMdiSubWindow : public QMdiSubWindow {
public:
  PersistentMdiSubWindow(const QString& name, QWidget* widget, QWidget* parent=NULL):
    QMdiSubWindow(parent),geom_handler_(name,this,parent)
  {
	this->init(name, widget);
  }

  PersistentMdiSubWindow(const QString& name, QWidget* widget, int width, int height, int x, int y, QWidget* parent=NULL):
    QMdiSubWindow(NULL),geom_handler_(name,this,parent)
  {
    this->init(name, widget);
    this->resize(width, height);
    this->move(x, y);
  }

  void init(const QString& name, QWidget* widget){
    this->setWidget(widget);
    widget->setParent(this);
  }

  void LoadGeom()
  {
    geom_handler_.LoadGeom("ui/perspective/main_area/");
  }

  void SaveGeom()
  {
    geom_handler_.SaveGeom("ui/perspective/main_area/");
  }

  ~PersistentMdiSubWindow()
  {
    this->SaveGeom();
  }

private:
  WidgetGeomHandler geom_handler_;
};

}

MainArea::MainArea(QWidget* parent):
  QMdiArea(parent)
{
  
}

void MainArea::AddPersistentWidget(const QString& title, const QString& name, 
                                   QWidget* widget, int width, int height, int x, int y)
{
  PersistentMdiSubWindow* pmdi=new PersistentMdiSubWindow(name, widget, width, height, x, y);
  this->addSubWindow(pmdi);
  pmdi->setWindowTitle(title);
  if (this->viewMode()==QMdiArea::SubWindowView) {
    pmdi->LoadGeom();
  }
  pmdi->show();
}

void MainArea::AddPersistentWidget(const QString& title, const QString& name,
                                   QWidget* widget, Qt::WindowStates window_state)
{
  PersistentMdiSubWindow* pmdi=new PersistentMdiSubWindow(name, widget, this);
  this->addSubWindow(pmdi);
  pmdi->setWindowTitle(title);
  if (this->viewMode()==QMdiArea::SubWindowView) {
    pmdi->LoadGeom();
  }
  pmdi->setWindowState(window_state);
  pmdi->show();
}

void MainArea::EnableTabbedMode(bool tm)
{
  
  QMdiArea::ViewMode mode=tm ? QMdiArea::TabbedView : QMdiArea::SubWindowView;
  this->setViewMode(mode);
  //this->setDocumentMode(tm);
}

void MainArea::AddWidget(const QString& title, QWidget* widget)
{
  QMdiSubWindow* mdi=new QMdiSubWindow(this);
  mdi->setAttribute(Qt::WA_DeleteOnClose);
  mdi->setWindowTitle(title);
  mdi->setWidget(widget);
  widget->setParent(mdi);
  this->addSubWindow(mdi);
  connect(widget, SIGNAL(destroyed()), mdi, SLOT(close()));
  mdi->show();
}

void MainArea::ShowSubWindow(QWidget* widget)
{
  assert(widget->parent());
  dynamic_cast<QWidget*>(widget->parent())->show();
}

void MainArea::HideSubWindow(QWidget* widget)
{
  assert(widget->parent());
  dynamic_cast<QWidget*>(widget->parent())->hide();
}

}}
