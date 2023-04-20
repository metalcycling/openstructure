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

#include "event_button.hh"
#include <QMimeData>
#include <QDrag>

namespace ost { namespace gui {

EventButton::EventButton(const QString& name, Widget* widget, bool pressed, QWidget* parent):
  QToolButton(parent)
{
  this->setText(name);
  connect(this, SIGNAL(clicked()), this,
  SLOT(clicked()));
  widget_ = widget;
  this->setCheckable(true);
  this->setChecked(pressed);
  this->setAcceptDrops(true);
  this->setToolTip("Show/Hide");
}

void EventButton::clicked(){
  emit this->ButtonPressed(widget_);
}

void EventButton::mouseMoveEvent(QMouseEvent* event){
   emit this->ButtonDragged(this);

   QMimeData* mimeData = new QMimeData();
   QDrag* drag = new QDrag(this);
   drag->setMimeData(mimeData);
   drag->setHotSpot(event->pos() - this->rect().topLeft());
   drag->exec(Qt::MoveAction);
}

void EventButton::dropEvent(QDropEvent* event){
  event->setDropAction(Qt::MoveAction);
  event->accept();
  emit this->ButtonDropped(this);
}

void EventButton::dragEnterEvent(QDragEnterEvent* event){
  event->accept();
}

Widget* EventButton::GetWidget(){
  return widget_;
}

}}
