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

#include <ost/platform.hh>

#include <ost/log.hh>

#include <ost/gui/widget_registry.hh>
#include <ost/gui/perspective.hh>
#include <ost/gui/file_browser.hh>
#include <ost/gui/main_area.hh>
#include <ost/gui/messages/message_box_widget.hh>

#include <ost/gui/panels/panel_manager.hh>
#include <ost/gui/panels/button_bar.hh>
#include <ost/gui/panels/panel_bar.hh>

#include <QTextEdit>
#include <QSizeGrip>
#include <QMainWindow>
#include <QSettings>
#include <QDebug>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QToolButton>
#include <QSizeGrip>
#include <QKeySequence>
#include <QStatusBar>
#include <QPushButton>
#include <QMenuBar>
#include <QMap>
#include <QString>
#include <QStatusBar>
#include <QSlider>

/*
  Author: Marco Biasini
 */
namespace ost { namespace gui {

Perspective::Perspective(QMainWindow* parent):
  QObject(parent), central_(new QWidget(parent)), 
  menu_bar_(new QMenuBar(parent)),
  main_area_(new MainArea),
  quick_access_bar_(new QWidget),
  verbosity_slider_(new QSlider(Qt::Horizontal,quick_access_bar_))
{
#if !defined(__APPLE__)
  parent->setMenuBar(menu_bar_);
#endif

  QVBoxLayout* l=new QVBoxLayout(central_);
  l->setMargin(0);
  l->setSpacing(0);
  parent->setCentralWidget(central_);

  status_bar_ = new QStatusBar(main_area_);
  status_bar_->setSizeGripEnabled (false);

  panels_ = new PanelManager(main_area_);
  l->addWidget(panels_,1);
  this->SetupQuickAccessBar();
  l->addWidget(quick_access_bar_, 0);
}

void Perspective::SetupQuickAccessBar()
{
  QDir icon_path(GetSharedDataPath().c_str());
  icon_path.cd("gui");
  icon_path.cd("icons");

  QHBoxLayout* l2=new QHBoxLayout(quick_access_bar_);
  l2->setMargin(0);
  l2->setSpacing(0);
  l2->setContentsMargins(0, 0, 0, 0);
  QPushButton* toggle_side_bar=new QPushButton(quick_access_bar_);
  toggle_side_bar->setAttribute(Qt::WA_MacSmallSize);
  toggle_side_bar->setFlat(true);
  toggle_side_bar->setToolTip("Glory A-Button");
  toggle_side_bar->setIcon(QIcon(icon_path.absolutePath()+QDir::separator()+
                                 QString("show_sidebar_icon.png")));
  toggle_side_bar->setIconSize(QSize(10,10));  
  toggle_side_bar->setCheckable(true);
  toggle_side_bar->setFixedSize(QSize(20, 20));  
  toggle_side_bar->setChecked(false);
  connect(toggle_side_bar, SIGNAL(clicked()), panels_,
          SLOT(ToggleHide()));
  l2->addWidget(toggle_side_bar, 0);

  QPushButton* add_side_bar_widget=new QPushButton(quick_access_bar_);
  add_side_bar_widget->setAttribute(Qt::WA_MacSmallSize);
  add_side_bar_widget->setFlat(true);  
  add_side_bar_widget->setFixedSize(QSize(20, 20));
  add_side_bar_widget->setToolTip("Add widget to left sidebar");
  add_side_bar_widget->setIcon(QIcon(icon_path.absolutePath()+QDir::separator()+
                                     QString("add_icon.png")));
  add_side_bar_widget->setIconSize(QSize(10,10));  
  connect(add_side_bar_widget, SIGNAL(clicked()), this, 
          SLOT(AddSideBarWidget()));          

  l2->addWidget(add_side_bar_widget, 0);

  verbosity_slider_->setMinimum(0);
  verbosity_slider_->setMaximum(5);
  verbosity_slider_->setValue(Logger::Instance().GetVerbosityLevel());
  verbosity_slider_->setToolTip("Verbosity Level");
  verbosity_slider_->setTickPosition(QSlider::TicksBothSides);
  verbosity_slider_->setFixedWidth(100);
  verbosity_slider_->setTracking(false);
  connect(verbosity_slider_, SIGNAL(valueChanged(int)),
          this,SLOT(VerbosityChanged(int)));
  l2->addWidget(verbosity_slider_, 0);

  l2->addWidget(status_bar_);
  l2->addWidget(new MessageBoxWidget(quick_access_bar_));
  l2->addWidget(new QSizeGrip(quick_access_bar_));
}

void Perspective::StatusMessage(const String& m)
{
  status_bar_->showMessage(QString(m.c_str()));
}

QMenuBar* Perspective::GetMenuBar()
{
  return menu_bar_;
}

QMenu* Perspective::GetMenu(const QString& name)
{
  QList<QAction *> actions=menu_bar_->actions();
  for (QList<QAction*>::iterator i=actions.begin(),
       e=actions.end(); i!=e; ++i) {
    if ((*i)->text()==name) {
      return (*i)->menu();
    }
  }
  QMenu* new_menu=new QMenu(name, menu_bar_);
  menu_bar_->addMenu(new_menu);
  return new_menu;
}

PanelManager* Perspective::GetPanels()
{
  return panels_;
}

void Perspective::AddSideBarWidget()
{
  QMenu* sb=panels_->GetAvailableWidgetsMenu();
  if (!sb) {
    return;
  }
  QPoint pos=dynamic_cast<QWidget*>(sender())->mapToGlobal(QPoint());
  if (QAction* a=sb->exec(pos)) {
    QString str_val=a->data().toString();
    WidgetRegistry* wf=WidgetRegistry::Instance();    
    panels_->AddWidget(LEFT_PANEL,wf->Create(str_val, NULL));
  }
  delete sb;
}

void Perspective::ChangeVerbositySlider(int value)
{
  verbosity_slider_->setValue(value);
}


void Perspective::VerbosityChanged(int value)
{
  Logger::Instance().PushVerbosityLevel(value);
}

bool Perspective::Restore()
{
  return true;
}

MainArea* Perspective::GetMainArea()
{
  return main_area_;
}

void Perspective::HideAllBars()
{
  panels_->ToggleHide();
  menu_bar_->hide();
}

bool Perspective::Save()
{
  return panels_->Save("ui/perspective/panels");
}

Perspective::~Perspective()
{
  this->Save();
}

}}
