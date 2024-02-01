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

#include "python_namespace_tree_item.hh"
#include <ost/base.hh>




namespace ost{namespace gui{

PythonNamespaceTreeItem::PythonNamespaceTreeItem(const bp::object& ns, 
                                                 const QString& name, 
                                                 PythonNamespaceTreeItem* parent):
  parent_(parent),
  children_(),
  namespace_(ns),
  name_(name),
  initialized_(false)
{
}

PythonNamespaceTreeItem::~PythonNamespaceTreeItem()
{
  DeleteChildren();
}

void PythonNamespaceTreeItem::DeleteChildren()
{
  qDeleteAll(children_.begin(), children_.end());
  children_.clear();
  initialized_=false;
}
unsigned int PythonNamespaceTreeItem::ChildCount()
{
  if(CanFetchMore()){
    FetchMore();
  }
  return children_.size();
}

PythonNamespaceTreeItem* PythonNamespaceTreeItem::GetParent() const
{
  return parent_;
}
PythonNamespaceTreeItem* PythonNamespaceTreeItem::GetChild(unsigned int index) 
{
  if(CanFetchMore()){
    FetchMore();
  }
  return children_.value(index);
}
unsigned int PythonNamespaceTreeItem::GetRow() const
{
  if(parent_){
    return parent_->children_.indexOf(const_cast<PythonNamespaceTreeItem*>(this));
  }
  return 0;
}

QString PythonNamespaceTreeItem::GetName()  const
{
  return name_;
}

bool PythonNamespaceTreeItem::HasChildren() const 
{
  if (initialized_) {
    return ! children_.empty();
  } else {
    return true;
  }
}

bool PythonNamespaceTreeItem::CanFetchMore() const
{
  return !initialized_;
}

void PythonNamespaceTreeItem::FetchMore()
{
  // todo should imediately return if worker thread is working
  // todo fix completion for builtins
  initialized_=true;  
  bp::object dir=bp::import("__main__").attr("__builtins__").attr("dir");  
  bp::list keys=bp::extract<bp::list>(dir(namespace_)); 
  unsigned int dict_length=bp::len(keys);
  for (unsigned int i=0;i<dict_length;++i) {
    QString child_name=QString::fromStdString(bp::extract<String>(keys[i]));
    if(child_name.startsWith("__")){
      continue;
    }
    bp::object child_namespace;
    try{
      String keystring=bp::extract<String>(keys[i]);
      child_namespace=namespace_.attr(keystring.c_str());
    } catch(bp::error_already_set&) {
      PyErr_Clear();
    }
    children_.append(new PythonNamespaceTreeItem(child_namespace,
                                                 child_name, this));
  }
}

}}//ns
