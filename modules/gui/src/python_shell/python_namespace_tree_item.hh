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
#ifndef PYTHON_NAMESPACE_TREE_ITEM_HH_
#define PYTHON_NAMESPACE_TREE_ITEM_HH_

/*
   Author: Andreas Schenk
 */
#include <boost/python.hpp>
#include <QList>
#include <QString>

namespace ost{namespace gui{
namespace bp = boost::python;

class PythonNamespaceTreeItem
{
public:
  PythonNamespaceTreeItem(const bp::object& ns, const QString& name, PythonNamespaceTreeItem* parent=0);
  virtual ~PythonNamespaceTreeItem();
  PythonNamespaceTreeItem* GetParent() const;
  PythonNamespaceTreeItem* GetChild(unsigned int index);
  unsigned int ChildCount();
  void DeleteChildren();
  unsigned int GetRow() const;
  QString GetName()  const;
  bool HasChildren() const;
  bool CanFetchMore() const;
  void FetchMore();
protected:
  PythonNamespaceTreeItem* parent_;
  QList<PythonNamespaceTreeItem*> children_;
  bp::object namespace_;
  QString name_;
  bool initialized_;
};

}}//ns

#endif /*PYTHON_NAMESPACE_TREE_ITEM_HH_*/
