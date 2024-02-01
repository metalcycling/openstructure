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
#include <ost/platform.hh>

#include <ost/gui/widget.hh>
#include <ost/gui/widget_registry.hh>
#include <ost/gui/perspective.hh>

#include "panel_manager.hh"
#include "panel_bar.hh"
#include "panel_bar_widget_holder.hh"

#include <cassert>
#include <QApplication>
#include <QMenu>
#include <QHBoxLayout>
#include <QToolButton>
#include <QDir>
#include <QIcon>
#include <QAction>
#include <QToolBar>
#include <QVBoxLayout>
#include <QAction>
#include <QComboBox>
#include <QDebug>
#include <QMimeData>
#include <QDrag>

/*
  Author: Marco Biasini, Stefan Scheuber
 */

namespace ost { namespace gui {


PanelBarWidgetHolder::PanelBarWidgetHolder(Widget* widget, Qt::Orientation orientation,
                                           QWidget* parent):
  QWidget(parent), toolbar_(new QToolBar(parent)),
  widget_(widget), label_(NULL)
{
  this->SetupToolBar();
  QVBoxLayout* v_layout=new QVBoxLayout(this);
  this->setLayout(v_layout);
  widget->setVisible(true);  
  v_layout->setMargin(0);
  v_layout->setSpacing(0);    
  v_layout->addWidget(toolbar_, 0);
  v_layout->addWidget(widget, 1);

  this->setSizePolicy(widget->sizePolicy());  
  if(orientation == Qt::Horizontal){
    if(this->sizePolicy().horizontalPolicy()==QSizePolicy::Fixed){
      this->setFixedWidth(widget->minimumWidth());
    } else if (this->sizePolicy().verticalPolicy()==QSizePolicy::Fixed) {
      v_layout->addStretch(1);
    }
  }
  if(orientation == Qt::Vertical){
    if(this->sizePolicy().verticalPolicy()==QSizePolicy::Fixed) {
      this->setFixedHeight(this->sizeHint().height());
    }
  }

  this->setAcceptDrops(true);
}

void PanelBarWidgetHolder::SetupToolBar()
{
  QDir icon_path(GetSharedDataPath().c_str());
  icon_path.cd("gui");
  icon_path.cd("icons");

  QString class_name=widget_->GetUniqueID();
  WidgetRegistry* wf=WidgetRegistry::Instance();
  toolbar_->setAttribute(Qt::WA_MacSmallSize);
  QAction* label_ = toolbar_->addAction(wf->GetFullName(class_name));
  connect(label_, SIGNAL(triggered(bool)), this, SLOT(LabelClick(bool)));
  //toolbar_->addSeparator();
  QWidget* spacer = new QWidget(this);
  QHBoxLayout* hbox = new QHBoxLayout(spacer);
  spacer->setLayout(hbox);
  hbox->addStretch();
  ActionList actions = widget_->GetActions();
  toolbar_->addWidget(spacer);
  toolbar_->setFloatable(false);
  toolbar_->setContentsMargins(0, 0, 0, 0);
  toolbar_->setToolButtonStyle(Qt::ToolButtonIconOnly);
  toolbar_->setIconSize(QSize(16,16));
  for (int i=0; i<actions.size();i++) {
    toolbar_->addAction(actions[i]);
  }  
  QAction* delete_me=toolbar_->addAction(" x ");
  delete_me->setToolTip("Close "+ label_->text());
  delete_me->setIcon(QIcon(icon_path.absolutePath()+QDir::separator()+
                           QString("close_icon.png")));
  connect(delete_me, SIGNAL(triggered()), this, SLOT(Remove()));
}

void PanelBarWidgetHolder::Remove()
{
  GostyApp::Instance()->GetPerspective()->GetPanels()->RemoveWidget(widget_);
}

void PanelBarWidgetHolder::SetWidget(Widget* widget)
{
  QBoxLayout* l=dynamic_cast<QBoxLayout*>(layout());
  assert(l);
  l->removeWidget(widget_);
  l->addWidget(widget, 1);
  QString class_name=widget->GetUniqueID();
  WidgetRegistry* wf=WidgetRegistry::Instance();
  label_->setText(wf->GetFullName(class_name));
  widget_=widget;
}

void PanelBarWidgetHolder::LabelClick(bool checked)
{
  PanelManager* panels = GostyApp::Instance()->GetPerspective()->GetPanels();
  QMenu* m = panels->GetAvailableWidgetsMenu();
  if (QAction* a=m->exec(toolbar_->mapToGlobal(QPoint(0, 0)))) {
    QString new_class_name=a->data().toString();
    if (new_class_name!=widget_->GetUniqueID()) {
      panels->ReplaceWidget(widget_,new_class_name);
    }
  }
}

PanelBarWidgetHolder::~PanelBarWidgetHolder(){
  widget_->setParent(NULL);
  widget_=NULL;
}

void PanelBarWidgetHolder::mousePressEvent(QMouseEvent *event)
{
  if (event->button() == Qt::LeftButton){
    drag_start_pos_ = event->pos();
  }
}

void PanelBarWidgetHolder::mouseMoveEvent(QMouseEvent *event)
{
  if (!(event->buttons() & Qt::LeftButton)) return;

  if ((event->pos() - drag_start_pos_).manhattanLength()
        < QApplication::startDragDistance()) return;

  QDrag *drag = new QDrag(this);
  QMimeData *mime_data = new QMimeData;
  QVariant self;
  self.setValue(widget_);
  mime_data->setData("OpenStructure/Widget",QByteArray());
  mime_data->setProperty("OpenStructure/Widget",self);

  drag->setMimeData(mime_data);
  drag->exec();
  PanelManager* panels = GostyApp::Instance()->GetPerspective()->GetPanels();
  panels->EndDrag();
}

void PanelBarWidgetHolder::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("OpenStructure/Widget")){
    PanelManager* panels = GostyApp::Instance()->GetPerspective()->GetPanels();
    panels->StartDrag();
    event->acceptProposedAction();
  }
}

void PanelBarWidgetHolder::dropEvent(QDropEvent *event)
{
  const QMimeData* mime_data = event->mimeData();
  QVariant variant =  mime_data->property("OpenStructure/Widget");
  Widget* widget = variant.value<Widget*>();
  if(widget && widget != widget_){
    GostyApp::Instance()->GetPerspective()->GetPanels()->MoveNextTo(widget_,widget);
  }

  event->acceptProposedAction();
}

}}
