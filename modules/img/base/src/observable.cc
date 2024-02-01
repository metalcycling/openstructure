//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  abstract observable concept

  Author: Ansgar Philippsen
*/

#include <ost/img/observable.hh>
#include <ost/img/data_observer.hh>

namespace ost { namespace img {

Observable::Observable() {
  list_.clear();
}

Observable::Observable(const Observable& o) {
  list_.clear();
}

Observable::~Observable() {
  for(ObserverIter it=list_.begin();it!=list_.end();++it) {
    (*it)->ObserverInvalidate();
    (*it)->ObserverRelease();
  }
}

Observable& Observable::operator=(const Observable& o) {
  list_.clear();
  return *this;
}

void Observable::Attach(DataObserver* d) {
  list_.push_back(d);
}

void Observable::Detach(DataObserver* o) {
  list_.remove(o);
}

void Observable::Notify() const {
  for(ObserverConstIter it=list_.begin();it!=list_.end();++it)
    (*it)->ObserverUpdate();
}
void Observable::Notify(const Extent& e) const {
  for(ObserverConstIter it=list_.begin();it!=list_.end();++it)
    (*it)->ObserverUpdate(e);
}
void Observable::Notify(const Point& p) const {
  for(ObserverConstIter it=list_.begin();it!=list_.end();++it)
    (*it)->ObserverUpdate(p);
}

int Observable::GetListSize() const {
  return list_.size();
}

long Observable::MemSize() const {
  return sizeof(list_) + list_.size()*sizeof(DataObserver*);
}


}} // namespace
