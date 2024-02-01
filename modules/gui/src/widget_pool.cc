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

#include "widget_pool.hh"
#include "widget_registry.hh"

namespace ost { namespace gui {

  WidgetPool::WidgetPool(QObject* parent):
    QObject(parent)
  { }

  void WidgetPool::Add(const QString& name, int limit)
  {
    QString n=name;
    n.replace(".", "::");
    if (pool_.contains(n)) {
      return;
    }
    pool_[n]=PooledWidget(n, limit);
  }

  bool WidgetPool::IsAvailable(const QString& name)
  {
    QString n=name;
    n.replace(".", "::");
    if (pool_.contains(n)) {
      int av=pool_[n].available;
      return av==-1 || av>0;
    }
    return false;
  }

  Widget* WidgetPool::Create(const QString& name)
  {
    QString n=name;
    n.replace(".", "::");
    if (this->IsAvailable(n)) {
      return WidgetRegistry::Instance()->Create(n, NULL);
    }
    return NULL;
  }

  bool WidgetPool::Take(const QString& name)
  {
    QString n=name;
    n.replace(".", "::");
    if (this->IsAvailable(n)) {
      int av=pool_[n].available;
      pool_[n].available=av==-1 ? av : av-1;
      return true;
    }
    return false;
  }

  bool WidgetPool::Give(const QString& name)
  {
    QString n=name;
    n.replace(".", "::");
    if (pool_.contains(n)) {
      int av=pool_[n].available;
      pool_[n].available=av==-1 ? av : av+1;
      return true;
    }
    return false;
  }

  QList<QString> WidgetPool::GetAvailableWidgets()
  {
    QList<QString> av;
    for (QMap<QString, PooledWidget>::iterator i=pool_.begin(),
         e=pool_.end(); i!=e; ++i) {

      if (i.value().available==-1 || i.value().available>0) {
        av.push_back(i.key());
      }
    }
    return av;
  }

}}
