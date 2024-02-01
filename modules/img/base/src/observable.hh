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
  observable concept

  Author: Ansgar Philippsen
*/

#ifndef IMG_OBSERVABLE_H
#define IMG_OBSERVABLE_H

#include <list>
#include "extent.hh" 

namespace ost { namespace img {

// fw decl
class DataObserver;

//! templated observable class
/*
  manages a list of observers, which must
  offer the methods ObserverRelease, ObserverInvalidate
  and ObserverUpdate
*/
class DLLEXPORT Observable {
  typedef std::list<DataObserver *> ObserverList;
  typedef  ObserverList::iterator ObserverIter;
  typedef  ObserverList::const_iterator ObserverConstIter;
public:
  Observable() ;
  /*
    copy logic: the observers are not copied
  */
  Observable(const Observable& o) ;

  ~Observable() ;

  /*
    assignement logic: the observers are not copied
  */
  Observable& operator=(const Observable& o) ;

  void Attach(DataObserver* d) ;

  void Detach(DataObserver* o) ;

  void Notify() const ;
  void Notify(const Extent& e) const ;
  void Notify(const Point& p) const ;
  int GetListSize() const ;

  long MemSize() const ;

  private:
    ObserverList list_;
};


}} // namespace

#endif
