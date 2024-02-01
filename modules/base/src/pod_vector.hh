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
#ifndef OST_POD_VECTOR_HH
#define OST_POD_VECTOR_HH

/*
   Author: Marco Biasini
 */

#include <ost/module_config.hh>

#include <cassert>
#include <memory>
#include <iterator>
#include <string>
#include <iostream>
namespace ost {

/// \brief vector container that treats its data as POD - even if it isn't in 
///      the strict sense.
template <typename T>
class TEMPLATE_DEF_EXPORT PodVector {
public:
  PodVector(): begin_(NULL), end_(NULL), capacity_(NULL) {}
  PodVector(const PodVector<T>& rhs)
  {
    this->do_reserve(rhs.capacity());
    this->raw_set(rhs.data(), rhs.size());
  }
  ~PodVector()
  {
    if (begin_) {
      free(begin_);
    }
  }
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;
  typedef T* iterator;
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef const T* const_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  
  iterator begin() { return begin_; }
  iterator end() { return end_; }
  const_iterator begin() const { return begin_; }
  const_iterator end() const { return end_; }
  
  reverse_iterator rbegin() { return reverse_iterator(begin_); }
  reverse_iterator rend() { return reverse_iterator(end_); }
  const_reverse_iterator rbegin() const 
  { 
    return const_reverse_iterator(begin_); 
  }
  const_reverse_iterator rend() const 
  { 
    return const_reverse_iterator(end_); 
  }
  
  size_t capacity() const { return capacity_-begin_; }
  size_t size() const { return end_-begin_; }
  
  reference front() 
  { 
    return *begin_; 
  }
  const_reference front() const 
  { 
    return *begin_; 
  }
  
  void push_back(const_reference val)
  {
    if (end_>=capacity_) {
      this->do_reserve(this->capacity()*2);
    }
    memcpy(end_, &val, sizeof(value_type));
    ++end_;
  }
  
  void pop_back()
  {
    --end_;
  }
  
  void clear()
  {
    this->do_reserve(0);
    end_=begin_;
  }
  
  void reserve(size_t n)
  {
    return do_reserve(n);
  }
  
  reference operator[](size_t index)
  {
    return begin_[index];
  }
  const_reference operator[](size_t index) const
  {
    return begin_[index];
  }
  T* data() { return begin_; }

  const T* data() const { return begin_; }

  void raw_set(T* data, size_t n)
  {
    this->do_reserve(n);
    memcpy(begin_, data, sizeof(T)*n);
    end_=begin_+n;
  }
private:
  void do_reserve(size_t n);
  value_type*    begin_;
  value_type*    end_;
  value_type*    capacity_;
};

template <typename T>
void PodVector<T>::do_reserve(size_t n)
{
  size_t old_size=this->size();
  n=(n==0 ? 1 : n);
  if (n==this->capacity()) {
    return;
  }
  begin_=reinterpret_cast<T*>(realloc(begin_, n*sizeof(T)));
  end_=begin_+old_size;
  capacity_=begin_+n;  
}

}

#endif
