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
 Author: Stefan Scheuber
 */

#include <ost/gui/gosty_app.hh>
#include <ost/platform.hh>

#include <ost/gui/perspective.hh>
#include <ost/gui/widget_registry.hh>
#include <ost/gui/messages/message_widget.hh>
#include <ost/gui/panels/panel_manager.hh>

#include "log_reader.hh"
#include "message_box_widget.hh"

#include <QDir>
#include <QLabel>
#include <QMenu>
#include <QVBoxLayout>


namespace ost {namespace gui {

MessageLevel::MessageLevel(QPixmap icon, QWidget* parent): QWidget(parent), message_count_(0), count_label_(new QLabel(QString::number(message_count_)+"  ",this)){
  QHBoxLayout* layout = new QHBoxLayout(this);
  layout->setMargin(0);
  layout->setSpacing(0);
  QLabel* icon_label = new QLabel(this);
  icon_label->setFixedSize(20,20);
  icon_label->setPixmap(icon.scaled(15,15, Qt::IgnoreAspectRatio, 
                                    Qt::SmoothTransformation));
  layout->addWidget(icon_label);
  layout->addWidget(count_label_);
}

void MessageLevel::SetMessageCount(int message_count){
  message_count_ = message_count;
  count_label_->setText(QString::number(message_count_)+"  ");
}

int MessageLevel::GetMessageCount(){
  return message_count_;
}



MessageBoxWidget::MessageBoxWidget(QWidget* parent) :
  QWidget(parent) {
  QHBoxLayout* layout = new QHBoxLayout(this);
  layout->setMargin(0);
  layout->setSpacing(0);
  MessageWidget* message_widget = GostyApp::Instance()->GetMessageWidget();
  this->connect(message_widget,SIGNAL(MessageCountChanged(QMessageBox::Icon)),this,SLOT(Update(QMessageBox::Icon)));
  this->connect(message_widget,SIGNAL(AllCleared()),this,SLOT(UpdateAll()));
}

void MessageBoxWidget::Update(QMessageBox::Icon icon) {
  MessageWidget* message_widget = GostyApp::Instance()->GetMessageWidget();
  int count = message_widget->GetMessagesCount(icon);
  if(!level_map_.contains(icon)){
    level_map_[icon] = new MessageLevel(message_widget->GetIcon(icon,this),this);
    this->layout()->addWidget(level_map_[icon]);
  }
  level_map_[icon]->SetMessageCount(count);
  if(count > 0){
    level_map_[icon]->show();
  }
  else{
    level_map_[icon]->hide();
  }
}

void MessageBoxWidget::UpdateAll() {
  QMapIterator<QMessageBox::Icon,MessageLevel*> i(level_map_);
  while (i.hasNext()) {
      i.next();
      this->Update(i.key());
  }
}

void MessageBoxWidget::mouseDoubleClickEvent(QMouseEvent* event){
  MessageWidget* message_widget = GostyApp::Instance()->GetMessageWidget();
  PanelManager* manager = GostyApp::Instance()->GetPerspective()->GetPanels();
  if(message_widget->isHidden()){
    manager->MoveWidget(message_widget,RIGHT_PANEL,0);
  }
}

MessageBoxWidget::~MessageBoxWidget() {
}

}
} // ns
