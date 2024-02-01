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
#include "python_interpreter.hh"
#include "python_namespace_tree_model.hh"
#include "python_namespace_tree_item.hh"

namespace ost{namespace gui{

PythonNamespaceTreeModel::PythonNamespaceTreeModel():
  QAbstractItemModel(),
  root_item_(new PythonNamespaceTreeItem(PythonInterpreter::Instance().GetMainModule(),"main"))
{
  connect(&PythonInterpreter::Instance(), SIGNAL(Finished(unsigned int,bool)),
          this, SLOT(NamespaceChanged(void)));
}

PythonNamespaceTreeModel::~PythonNamespaceTreeModel()
{
  delete root_item_;
}

void PythonNamespaceTreeModel::NamespaceChanged()
{
  // todo should only be called after all commands are executed
  delete root_item_;
  root_item_=new PythonNamespaceTreeItem(PythonInterpreter::Instance().GetMainModule(),"main");
  dataChanged(QModelIndex(),QModelIndex());
}

QModelIndex PythonNamespaceTreeModel::index(int row, int column, const QModelIndex &parent) const
{
  if (!hasIndex(row, column, parent)){
    return QModelIndex();
  }
  PythonNamespaceTreeItem* parent_item;
  if (!parent.isValid()){
    parent_item = root_item_;
  }else{
    parent_item = static_cast<PythonNamespaceTreeItem*>(parent.internalPointer());
  }
  PythonNamespaceTreeItem *child_item = parent_item->GetChild(row);
  if (child_item){
    return createIndex(row, column, child_item);
  }else{
    return QModelIndex();
  }
}

QModelIndex PythonNamespaceTreeModel::parent(const QModelIndex &index) const
{
  if (!index.isValid()){
    return QModelIndex();
  }
  PythonNamespaceTreeItem* child_item=static_cast<PythonNamespaceTreeItem*>(index.internalPointer());
  PythonNamespaceTreeItem* parent_item=child_item->GetParent();
  if (parent_item == root_item_){
    return QModelIndex();
  }
  return createIndex(parent_item->GetRow(), 0, parent_item);
}

int PythonNamespaceTreeModel::rowCount(const QModelIndex &parent) const
 {
  PythonNamespaceTreeItem *parent_item;
  if (parent.column() > 0){
    return 0;
  }

  if (!parent.isValid()){
    parent_item = root_item_;
  }else{
    parent_item = static_cast<PythonNamespaceTreeItem*>(parent.internalPointer());
  }
  return parent_item->ChildCount();
 }

int PythonNamespaceTreeModel::columnCount(const QModelIndex &parent) const
{
  return 1;
}

QVariant PythonNamespaceTreeModel::data(const QModelIndex &index, int role) const
{
  if (!index.isValid()|| !(role == Qt::DisplayRole || role == Qt::EditRole)){
    return QVariant();
  }
  PythonNamespaceTreeItem* item = static_cast<PythonNamespaceTreeItem*>(index.internalPointer());
  return QVariant(item->GetName());
}

Qt::ItemFlags PythonNamespaceTreeModel::flags(const QModelIndex &index) const
{
  if (!index.isValid()){
    return Qt::ItemFlags();
  }
  return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}


QVariant PythonNamespaceTreeModel::headerData(int section, 
                                              Qt::Orientation orientation,
                                              int role) const
 {
  if (orientation == Qt::Horizontal && role == Qt::DisplayRole){
    return QVariant(root_item_->GetName());
  }
  return QVariant();
 }

bool PythonNamespaceTreeModel::canFetchMore( const QModelIndex & parent ) const
{
  PythonNamespaceTreeItem *parent_item;
  if (!parent.isValid()){
    parent_item = root_item_;
  }else{
    parent_item = static_cast<PythonNamespaceTreeItem*>(parent.internalPointer());
  }
  return parent_item->CanFetchMore();
}

void PythonNamespaceTreeModel::fetchMore( const QModelIndex & parent )
{
  PythonNamespaceTreeItem *parent_item;
  if (!parent.isValid()){
    parent_item = root_item_;
  }else{
    parent_item = static_cast<PythonNamespaceTreeItem*>(parent.internalPointer());
  }
  int first=parent_item->ChildCount();
  parent_item->FetchMore();  
  this->beginInsertRows(parent, first, parent_item->ChildCount());  
  this->endInsertRows();
}

bool PythonNamespaceTreeModel::hasChildren(const QModelIndex & parent) const
{
  PythonNamespaceTreeItem *parent_item;
  if (!parent.isValid()){
    parent_item = root_item_;
  }else{
    parent_item = static_cast<PythonNamespaceTreeItem*>(parent.internalPointer());
  }
  return parent_item->HasChildren();
}

}}//ns
