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
#include <iostream>
#include "widget_geom_handler.hh"

#include <QSettings>
namespace ost { namespace gui {

WidgetGeomHandler::WidgetGeomHandler(const QString& name, QWidget* widget, QWidget* parent)
{
  this->setObjectName(name);
  widget_ = widget;
}

void WidgetGeomHandler::LoadGeom(const QString& prefix)
{
  QSettings settings;
  settings.beginGroup(prefix+this->objectName());
  if (settings.contains("pos") && settings.contains("size")) {
    QPoint pos = settings.value("pos").toPoint();
    Qt::WindowStates state = Qt::WindowStates(settings.value("state").toInt());
    widget_->move(pos);
    widget_->setWindowState(state);
  }
  settings.endGroup();
}

void WidgetGeomHandler::SaveGeom(const QString& prefix)
{
  QSettings settings;
  settings.beginGroup(prefix+this->objectName());
  settings.setValue("pos", widget_->pos());
  settings.setValue("size", widget_->size());
  settings.setValue("state",int(widget_->windowState()));
  settings.endGroup();
}

}}
