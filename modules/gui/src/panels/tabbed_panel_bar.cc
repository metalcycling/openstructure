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
#include <ost/gui/perspective.hh>
#include <ost/gui/widget_registry.hh>
#include <ost/gui/panels/panel_manager.hh>
#include "tabbed_panel_bar.hh"

#include <QApplication>
#include <QSettings>
#include <QDir>
#include <QMimeData>
#include <QDrag>

namespace ost { namespace gui {

TabbedPanelBar::TabbedPanelBar(PanelBar* parent):
  PanelWidgetContainer(parent)
{
  layout_ = new QHBoxLayout(this);
  layout_->setMargin(0);
  layout_->setSpacing(0);
  tab_widget_ = new TabbedDragWidget(this);
  #if QT_VERSION>=0x40503
    tab_widget_->setDocumentMode(true);
  #endif
  layout_->addWidget(tab_widget_);

  toolbar_ = new QToolBar(this);
  toolbar_->setFloatable(false);
  toolbar_->setContentsMargins(0, 0, 0, 0);
  toolbar_->setToolButtonStyle(Qt::ToolButtonIconOnly);
  toolbar_->setIconSize(QSize(16,16));
  tab_widget_->setCornerWidget(toolbar_, Qt::TopRightCorner);
  connect(tab_widget_,SIGNAL(currentChanged(int)),this,SLOT(CurrentChanged(int)));
  this->setAcceptDrops(true);
}

bool TabbedPanelBar::Save(const QString& prefix)
{
  return true;
}

bool TabbedPanelBar::Restore(const QString& prefix)
{
  return true;
}

void TabbedPanelBar::WidgetMoved(Widget* widget, int index){
  QString class_name = widget->GetUniqueID();
  if(tab_widget_->indexOf(widget)>=0){
    tab_widget_->removeTab(tab_widget_->indexOf(widget));
  }
  tab_widget_->insertTab(index,widget,WidgetRegistry::Instance()->GetFullName(class_name));
}

void TabbedPanelBar::ShowWidget(Widget* widget, int pos, bool show){
  if(tab_widget_->indexOf(widget)>=0){
    if(!show){
      tab_widget_->removeTab(tab_widget_->indexOf(widget));
    }
  }
  else{
    if(show){
      WidgetMoved(widget,pos);
    }
  }
}

QString TabbedPanelBar::GetName(){
  return "Tabbed Mode";
}

void TabbedPanelBar::CurrentChanged(int index){
  if(Widget* widget = qobject_cast<Widget*>(tab_widget_->currentWidget())){
    ActionList actions = widget->GetActions();
    toolbar_->clear();
    if(!actions.isEmpty()){
      for(int i = 0;i<actions.size();i++){
        toolbar_->addAction(actions[i]);
      }
    }
  }
}

void TabbedPanelBar::dragEnterEvent(QDragEnterEvent *event)
{
  if (event->mimeData()->hasFormat("OpenStructure/Widget")){
    PanelManager* panels = GostyApp::Instance()->GetPerspective()->GetPanels();
    panels->StartDrag();
    event->acceptProposedAction();
  }
}

void TabbedPanelBar::dropEvent(QDropEvent *event)
{
  const QMimeData* mime_data = event->mimeData();
  QVariant variant =  mime_data->property("OpenStructure/Widget");
  Widget* widget = variant.value<Widget*>();
  Widget* current_widget = qobject_cast<Widget*>(tab_widget_->currentWidget());
  if(widget && current_widget && widget != current_widget){
    GostyApp::Instance()->GetPerspective()->GetPanels()->MoveNextTo(current_widget,widget);
  }
  event->acceptProposedAction();
}

TabbedPanelBar::~TabbedPanelBar(){
  toolbar_->clear();
}

//TabbedDragWidget
TabbedDragWidget::TabbedDragWidget(QWidget* parent):
  QTabWidget(parent) {
  setTabBar(new DragTabBar(this));
}

//DragTabBar
DragTabBar::DragTabBar(QTabWidget* parent):
  QTabBar(parent),tab_widget_(parent) {}

void DragTabBar::mousePressEvent(QMouseEvent *event)
{
  if (event->button() == Qt::LeftButton){
    drag_start_pos_ = event->pos();
  }
  QTabBar::mousePressEvent(event);
}

void DragTabBar::mouseMoveEvent(QMouseEvent *event)
{
  if (!(event->buttons() & Qt::LeftButton)) return;

  if ((event->pos() - drag_start_pos_).manhattanLength()
        < QApplication::startDragDistance()) return;

  Widget* widget = qobject_cast<Widget*>(tab_widget_->currentWidget());
  if(!widget)return;

  QDrag *drag = new QDrag(this);
  QMimeData *mime_data = new QMimeData;
  QVariant self;
  self.setValue(widget);
  mime_data->setData("OpenStructure/Widget",QByteArray());
  mime_data->setProperty("OpenStructure/Widget",self);

  drag->setMimeData(mime_data);
  drag->exec();
  PanelManager* panels = GostyApp::Instance()->GetPerspective()->GetPanels();
  panels->EndDrag();
}



}}

