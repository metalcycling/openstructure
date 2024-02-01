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
  Authors: Marco Biasini, Andreas Schenk
 */


#include "python_completer.hh"
#include "python_namespace_tree_model.hh"
#include <iostream>
#include <QStringList>


namespace ost { namespace gui {



PythonCompleter::PythonCompleter(QObject * parent):
 OstCompleterBase(parent),
 model_(new PythonNamespaceTreeModel())
{
  setModel(model_);
}

PythonCompleter::~PythonCompleter()
{
  delete model_;
}

QStringList PythonCompleter::splitPath(const QString & path) const
{
  return path.split(".");
}

QString PythonCompleter::pathFromIndex(const QModelIndex & index) const
{
  if(!index.isValid()){
    return QString();
  }
  QString path=model()->data(index,Qt::DisplayRole).toString();
  QModelIndex parent_index=model()->parent(index);
  while(parent_index.isValid()){
    path=model()->data(parent_index,Qt::DisplayRole).toString()+"."+path;
    parent_index=model()->parent(parent_index);
  }
  return path;
}

void PythonCompleter::SetCompletionPrefix(const QString& prefix)
{
  dynamic_cast<PythonNamespaceTreeModel*>(this->model())->Reset();  
  this->setCompletionPrefix(prefix);
}

}}//ns

