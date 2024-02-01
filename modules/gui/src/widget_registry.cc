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

#include "widget_registry.hh"

#include <QDebug>
namespace ost { namespace gui {
  

bool WidgetRegistry::RegisterWidgetFactory(WidgetFactory* factory)
{
  factory->setParent(this);
  id_to_factory_.insert(factory->GetUniqueID(), factory);
  return true;
}

bool WidgetRegistry::UnregisterWidgetFactory(const QString& id)
{
  if(id_to_factory_.contains(id)){
    id_to_factory_[id]->setParent(NULL);
    id_to_factory_.remove(id);
    return true;
  }
  return false;
}

/// \brief create widget by id
Widget* WidgetRegistry::Create(const QString& id, QWidget* parent)
{
  if (id_to_factory_.contains(id)) {
    return id_to_factory_.value(id)->Create(parent);
  }
  return NULL;
}

QString WidgetRegistry::GetFullName(const QString& id)
{
  if (id_to_factory_.contains(id)) {
    return id_to_factory_.value(id)->GetFullName();
  }
  return QString("");
}

WidgetRegistry* WidgetRegistry::Instance()
{
  static WidgetRegistry* registry=new WidgetRegistry;
  return registry;
}

}}
