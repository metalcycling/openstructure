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
#ifndef OST_GUI_QPTR_HH
#define OST_GUI_QPTR_HH

#include <iostream>
/*
   Author: Marco Biasini
 */
namespace ost { namespace gui {

/// \brief storage class for exporting objects to Python that use the
///    hierarchical memory management mechanism of Qt.
template <typename T>
class QPtr
{
public:
  typedef T  element_type;
  typedef T  value_type;
  typedef T* pointer_type;
  
public:
  template <typename P> friend class QPtr;
  QPtr(): ptr_(NULL), master_(true), released_(new bool(false))
  { }
  
  QPtr(T* ptr): ptr_(ptr), master_(true), released_(new bool(false))
  { 
  }
  
  void swap(QPtr<T>& rhs)
  {
    std::swap(ptr_, rhs.ptr_);
    std::swap(master_, rhs.master_);
    std::swap(released_, rhs.released_);
  }
  template <typename P>
  QPtr(QPtr<P>& rhs): ptr_(rhs.get()), master_(false), 
                      released_(rhs.released_)
  { }
  
  void Release()
  {
    *released_=true;
  }
  
  bool IsReleased() const 
  {
    return (*released_)==true; 
  }
  pointer_type get() 
  { 
    return ptr_; 
  }
  
  const pointer_type get() const
  {
    return ptr_;
  }
  
  ~QPtr()
  {
    if (!master_) {
      return;
    }
    if ((*released_)==false && ptr_) {
      delete ptr_;
      ptr_=NULL;
    }
    delete released_;
  }
  
  operator bool() const
  {
    return ptr_!=NULL;
  }
private:
  pointer_type   ptr_;
  bool           master_;
  bool*          released_;
}; 


template <typename T>
T* get_pointer(QPtr<T>& rhs)
{
  return rhs.get();
}

template <typename T>
const T* get_pointer(const QPtr<T>& rhs)
{
  return rhs.get();
}

}}

#endif
