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

#include "button_bar.hh"

namespace ost { namespace gui {

ButtonBar::ButtonBar(QWidget* parent):
 Widget(NULL, parent)
{
  layout_ = new QHBoxLayout();
  layout_->setMargin(0);
  layout_->setSpacing(0);
  this->setLayout(layout_);
  drag_box_=0;
  drag_button_=0;
}

bool ButtonBar::Save(const QString& prefix)
{
  return true;
}

bool ButtonBar::Restore(const QString& prefix)
{
  return true;
}

void ButtonBar::AddButtonBox(ButtonBox* button_box){
  layout_->addWidget(button_box);
  connect(button_box, SIGNAL(ButtonDragged(EventButton*)), this,
        SLOT(ButtonDragged(EventButton*)));
  connect(button_box, SIGNAL(ButtonDropped(EventButton*)), this,
        SLOT(ButtonDropped(EventButton*)));
}

void ButtonBar::ButtonDragged(EventButton* button){
  QObject* sender = const_cast<QObject*>(QObject::sender());
  drag_box_ = qobject_cast<ButtonBox*>(sender);
  drag_button_ = button;
}

void ButtonBar::ButtonDropped(EventButton* drop_button){
  QObject* sender = const_cast<QObject*>(QObject::sender());
  ButtonBox* drop_box = qobject_cast<ButtonBox*>(sender);
  if(drag_box_ && drop_box){
    if(drag_box_ == drop_box){
      drop_box->MoveButton(drop_button, drag_button_);
    }
    else{
      drag_box_->RemoveButton(drag_button_);
      drop_box->InsertButton(drop_button, drag_button_);
    }
  }
}

}}

