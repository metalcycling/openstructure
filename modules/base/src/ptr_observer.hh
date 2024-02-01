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

#ifndef OST_PTR_OBSERVER_HH
#define OST_PTR_OBSERVER_HH

#include <ost/message.hh>
#include <ost/module_config.hh>
/*
  ptr observer concept

  Author: Ansgar Philippsen
*/

class DLLEXPORT InvalidatedPointer: public ost::Error {
public:
  InvalidatedPointer():
    ost::Error("access attempt on invalidated Pointer")
  {}
};

template <class T>
class TEMPLATE_EXPORT PtrObserver {
public:
  PtrObserver(T* t):
    observed_(t),
    valid_(true)
  {}

  ~PtrObserver() {}

  void Invalidate() {valid_=false;}

  bool IsValid() const {return valid_;}

  T& Ref() {
    if(!valid_) throw InvalidatedPointer();
    return *observed_;
  }

  T* Ptr() {
    if(!valid_) throw InvalidatedPointer();
    return observed_;
  }

  bool operator==(const PtrObserver<T>& ref) const {
    return observed_ == ref.observed_;
  }

private:
  PtrObserver(const PtrObserver<T>& o);
  PtrObserver<T>& operator=(const PtrObserver<T>& u);

  T* observed_;
  bool valid_;
};

#endif
