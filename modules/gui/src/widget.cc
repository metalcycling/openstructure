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
#include <QVBoxLayout>

#include "widget.hh"

/*
  Author: Marco Biasini
 */

namespace ost { namespace gui {

Widget::Widget(QWidget* widget, QWidget* parent):
  QWidget(parent), internal_(widget), destroy_on_close_(true)
{
  if (widget) {
    widget->setParent(this);
  }
}

void Widget::SetDestroyOnClose(bool flag)
{
  destroy_on_close_=flag;
}

ActionList Widget::GetActions(){
  return QList<QAction *>();
}

QWidget* Widget::GetInternalWidget()
{
  return internal_;
}

void Widget::SetInternalWidget(QWidget* widget)
{
  assert(internal_==NULL);
  internal_=widget;
  QBoxLayout* l=new QVBoxLayout(this);
  l->setMargin(0);
  l->setSpacing(0);
  l->addWidget(widget);
  widget->setParent(this);
}

bool Widget::DestroyOnClose() const
{
  return destroy_on_close_;
}

}}
