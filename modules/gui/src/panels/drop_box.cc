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


#include "drop_box.hh"

#include <QHBoxLayout>
#include <QLabel>
namespace ost { namespace gui {

DropBox::DropBox(QWidget* parent):
 Widget(NULL, parent)
{
  QHBoxLayout* layout = new QHBoxLayout();
  this->setLayout(layout);
  QLabel* label = new QLabel("drop it here");
  layout->addWidget(label);
  label->setBackgroundRole(QPalette::Text);
  QPalette pal;
  QColor bgColor("red");
  QColor fgColor("green");
  pal.setColor(QPalette::Base, bgColor);

  pal.setColor(QPalette::Text, fgColor);
  label->setPalette(pal);
}

bool DropBox::Save(const QString& prefix)
{
  return true;
}

bool DropBox::Restore(const QString& prefix)
{
  return true;
}

}}

