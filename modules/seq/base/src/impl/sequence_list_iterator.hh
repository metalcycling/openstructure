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
#ifndef OST_SEQ_IMPL_SEQUENCE_LIST_ITERATOR_HH
#define OST_SEQ_IMPL_SEQUENCE_LIST_ITERATOR_HH

/*
  Author: Marco Biasini
 */

#include <ost/seq/module_config.hh>
namespace ost { namespace seq { namespace impl {

/// \internal
template <typename T, typename I>
class TEMPLATE_EXPORT SequenceListIterator 
  : public std::iterator<std::random_access_iterator_tag, T>{
public:
  typedef T ValueType;
  typedef I IteratorType;
  typedef SequenceListIterator<T, I> ClassType;
protected:
  void UpdateVal()
  {
    if (current_it_<end_) {
      current_val_=ValueType(*current_it_);
    }
  }  
public:
  SequenceListIterator(IteratorType it, IteratorType end):
    current_it_(it), end_(end)
  {
    this->UpdateVal();
  }
  
  ValueType& operator*()
  {
    return current_val_;
  }
  
  ValueType* operator->()
  {
    return &current_val_;
  }
  
  ClassType& operator++() {
    ++current_it_;
    this->UpdateVal();
    return *this;
  }
  
  ClassType& operator+=(int n) {
    current_it_+=n;
    this->UpdateVal();
    return *this;
  }
  ClassType operator+(int n) {
    ClassType ans(current_it_, end_);
    ans+=n;
    return *this;
  }
  ClassType& operator-=(int n) {
    current_it_-=n;
    this->UpdateVal();
    return *this;
  }
  
  ClassType operator++(int) {
    ClassType ans(*this);
    ++current_it_;
    this->UpdateVal();
    return ans;
  }  
  
  bool operator==(const ClassType& rhs) const 
  { 
    return current_it_==rhs.current_it_; 
  }
  
  bool operator!=(const ClassType& rhs) const 
  { 
    return !(*this==rhs); 
  }  
private:
  IteratorType current_it_;
  IteratorType end_;  
  ValueType    current_val_;
};

}}}

#endif
