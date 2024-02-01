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

#include "button_box.hh"

namespace ost { namespace gui {

ButtonBox::ButtonBox(const QString& label, QWidget* parent):
 Widget(NULL, parent)
{
  layout_ = new QHBoxLayout();
  layout_->setMargin(2);
  layout_->setSpacing(2);
  buttons_layout_ = new QHBoxLayout();
  label_ = new QLabel(label,this);
  this->setLayout(layout_);
  layout_->addWidget(label_);
  layout_->addLayout(buttons_layout_);
  buttons_ = new QList<EventButton *>();
  this->setAcceptDrops(true);
}

bool ButtonBox::Save(const QString& prefix)
{
  return true;
}

bool ButtonBox::Restore(const QString& prefix)
{
  return true;
}

void ButtonBox::AddButton(const QString& name, Widget* widget, bool pressed){
  EventButton* button = new EventButton(name,widget,pressed);
  buttons_layout_->addWidget(button);
  connect(button, SIGNAL(ButtonPressed(Widget*)), this,
        SLOT(ButtonPressed(Widget*)));
  connect(button, SIGNAL(ButtonDragged(EventButton*)), this,
        SIGNAL(ButtonDragged(EventButton*)));
  connect(button, SIGNAL(ButtonDropped(EventButton*)), this,
        SIGNAL(ButtonDropped(EventButton*)));
  buttons_->append(button);
  //button->setShortcut("Ctrl+Shift+"+QString::number(buttons_->size()));
}

void ButtonBox::ButtonPressed(Widget* widget){
  emit this->ToggleWidget(widget);
}

void ButtonBox::RemoveButton(EventButton* button){
  buttons_->removeAt(buttons_->indexOf(button));
  buttons_layout_->removeWidget(button);
  emit this->Removed(button->GetWidget());
}

void ButtonBox::InsertButton(EventButton* position, EventButton* new_button){
  int posIndex = buttons_->indexOf(position);
  if(posIndex == -1)posIndex = buttons_->size();
    buttons_->insert(posIndex,new_button);
    buttons_layout_->insertWidget(posIndex, new_button);
    new_button->disconnect(SIGNAL(ButtonPressed(Widget*)));
    new_button->disconnect(SIGNAL(ButtonDragged(EventButton*)));
    new_button->disconnect(SIGNAL(ButtonDropped(EventButton*)));
    connect(new_button, SIGNAL(ButtonPressed(Widget*)), this,
          SLOT(ButtonPressed(Widget*)));
    connect(new_button, SIGNAL(ButtonDragged(EventButton*)), this,
          SIGNAL(ButtonDragged(EventButton*)));
    connect(new_button, SIGNAL(ButtonDropped(EventButton*)), this,
          SIGNAL(ButtonDropped(EventButton*)));
    emit this->Inserted(new_button->GetWidget(),posIndex,new_button->isChecked());
}

void ButtonBox::MoveButton(EventButton* position, EventButton* mov_button){
  int position_index = buttons_->indexOf(position);
  int move_index = buttons_->indexOf(mov_button);
  if(move_index >= 0 && position_index >= 0){
    buttons_->removeAt(move_index);
    buttons_layout_->removeWidget(mov_button);
    buttons_->insert(position_index,mov_button);
    buttons_layout_->insertWidget(position_index, mov_button);
    emit this->Inserted(mov_button->GetWidget(),position_index,mov_button->isChecked());
  }
}

void ButtonBox::dropEvent(QDropEvent* event){
  event->setDropAction(Qt::MoveAction);
  event->accept();
  emit this->ButtonDropped(NULL);
}

void ButtonBox::dragEnterEvent(QDragEnterEvent* event){
  event->accept();
}


}}

