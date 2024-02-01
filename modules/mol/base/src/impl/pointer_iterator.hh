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
#ifndef OST_MOL_IMPL_POINTER_ITERATOR_HH
#define OST_MOL_IMPL_POINTER_ITERATOR_HH

#include <vector>

namespace ost { namespace mol { namespace impl {

template <typename T>
class pointer_it :  public std::iterator<std::forward_iterator_tag, T>{
public:
  pointer_it(T* s): s_(s) { }
  
  pointer_it<T>& operator++()
  {
    ++s_;
    return *this;
  }
  pointer_it<T>& operator+(int rhs)
  {
    s_+=rhs;
    return *this;
  }
  
  bool operator==(const pointer_it<T>& rhs) const
  {
    return rhs.s_==s_;
  }
  
  bool operator!=(const pointer_it<T>&rhs) const
  {
    return !(*this==rhs);
  }

  T& operator*()
  {
    return *s_;
  }

  T* operator->()
  {
    return s_;
  }
private:
  T* s_;
};

template <typename T>
inline pointer_it<T> begin(const std::vector<T>& values)
{
  return pointer_it<T>(values.empty() ? NULL : const_cast<T*>(&values.front()));
}

template <typename T>
inline pointer_it<T> end(const std::vector<T>& values)
{
  return pointer_it<T>(values.empty() ? NULL : const_cast<T*>(&values.back()+1));
}

}}}

#endif
