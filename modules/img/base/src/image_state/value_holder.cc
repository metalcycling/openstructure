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
  value holder for image states

  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_IMAGE_STATE_INST_H
#error "value_holder.cc must be included from image_state_inst.hh"
#endif

#include <cassert>
#include <algorithm>

#include <ost/img/value_util.hh>
#include <ost/img/size.hh>

#include "value_holder.hh"

namespace ost { namespace img { namespace image_state {

namespace {

template <typename V>
std::size_t calc_volume(const Size& s)
{
  unsigned int tmp;
  unsigned int numbits = 0;
  //find and add highest bit set for each dimension
  tmp=s[0];
  while (tmp >>= 1)
  {
    ++numbits;
  }
  tmp=s[1];
  while (tmp >>= 1)
  {
    ++numbits;
  }
  tmp=s[2];
  while (tmp >>= 1)
  {
    ++numbits;
  }
  tmp=sizeof(V);
  while (tmp >>= 1)
  {
    ++numbits;
  }
  // check if number of pixels (2**(numbits+1)) larger than what std::size_t can handle
  if(numbits>=sizeof(std::size_t)*8){
      throw std::bad_alloc();
  }
  return static_cast<std::size_t>(s[0])*static_cast<std::size_t>(s[1])*static_cast<std::size_t>(s[2]);
}

} //ns

template <typename V>
ValueHolder<V>::ValueHolder(const Size& s):
#ifdef USE_ROW_ORDER
  depth_(s[2]),
  height_depth_(s[1]*depth_),
  data_(s[0]*height_depth_)
#else
  height_(s[1]),
  width_height_(s[0]*height_),
    data_(calc_volume<V>(s))
#endif
{
}

template <typename V>
ValueHolder<V>::ValueHolder(const Size& s, const Size& ps):
#ifdef USE_ROW_ORDER
  depth_(s[2]),
  height_depth_(s[1]*depth_),
#else
  height_(s[1]),
  width_height_(s[0]*height_),
#endif
    data_(calc_volume<V>(ps))
{
}

template <typename V>
ValueHolder<V>::ValueHolder(const ValueHolder<V>& h):
#ifdef USE_ROW_ORDER
  depth_(h.depth_),
  height_depth_(h.height_depth_),
#else
  height_(h.height_),
  width_height_(h.width_height_),
#endif
  data_(h.data_.begin(),h.data_.end())
{
}

template <typename V>
ValueHolder<V>& ValueHolder<V>::operator=(const ValueHolder<V>& h) 
{
  if(this!=&h) {
    data_.assign(h.data_.begin(),h.data_.end());
    #ifdef USE_ROW_ORDER
      depth_=h.depth_;
      height_depth_=h.height_depth_;
    #else
      height_=h.height_;
      width_height_=h.width_height_;
    #endif
  }
  return *this;
}

template <typename V>
ValueHolder<V>::~ValueHolder()
{
}

template <typename V>
void ValueHolder<V>::Swap(ValueHolder& vh)
{
  #ifdef USE_ROW_ORDER
    std::swap<>(depth_,vh.depth_);
    std::swap<>(height_depth_,vh.height_depth_);
  #else
    std::swap<>(height_,vh.height_);
    std::swap<>(width_height_,vh.width_height_);
  #endif
  std::swap<>(data_,vh.data_);
}




template <typename V>
DataType ValueHolder<V>::GetDataType() 
{
  return Val2Type<V>();
}

template <typename V>
std::size_t ValueHolder<V>::MemSize() const
{
  return sizeof(V)*data_.size();
}



}}} // namespaces
