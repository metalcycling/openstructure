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
#ifndef PYTHON_NAMESPACE_TREE_MODEL_HH_
#define PYTHON_NAMESPACE_TREE_MODEL_HH_

/*
  Author: Andreas Schenk
 */

#include <QAbstractItemModel>


namespace ost{namespace gui{


//fw decl
class PythonNamespaceTreeItem;


class PythonNamespaceTreeModel: public QAbstractItemModel
{
  Q_OBJECT
public:
  PythonNamespaceTreeModel();
  virtual ~PythonNamespaceTreeModel();
  QVariant data(const QModelIndex &index, int role) const;
  Qt::ItemFlags flags(const QModelIndex &index) const;
  QVariant headerData(int section, Qt::Orientation orientation,
                      int role=Qt::DisplayRole) const;
  virtual bool hasChildren(const QModelIndex& parent=QModelIndex()) const;
  virtual bool canFetchMore(const QModelIndex& parent) const;
  virtual void fetchMore(const QModelIndex& parent);
  int rowCount(const QModelIndex& parent=QModelIndex()) const;
  int columnCount(const QModelIndex&parent=QModelIndex()) const;
  QModelIndex index(int row, int column, const QModelIndex& parent=QModelIndex()) const;
  QModelIndex parent(const QModelIndex &index) const;
  void Reset() {this->beginResetModel();
                this->endResetModel(); }
public slots:
  void NamespaceChanged();

protected:
  PythonNamespaceTreeItem* root_item_;
};

}}//ns

#endif /*PYTHON_NAMESPACE_TREE_MODEL_HH_*/
