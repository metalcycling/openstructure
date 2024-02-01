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
  value holder for image state

  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef VALUE_HOLDER_H
#define VALUE_HOLDER_H

#include <ost/base.hh>
#include <ost/img/module_config.hh>
#include <ost/img/data_types.hh>
#include <ost/img/size.hh>

#include "index.hh"
#include "type_fw.hh"

namespace value_holder_test {

template <typename V>
void ReleaseAndReconstruct();

}

#define USE_ROW_ORDER 1



namespace ost { namespace img { namespace image_state {

/*
  Value Holder

  Provide a dynamically allocated 3D array with
  optimized access, based on a zero based Index, which
  in turn encapsulates an unsigned integer triplet

  Should provide reasonable exception safety, ie an exception
  during construction will not cause a memory leak.
*/
template <typename V>
class TEMPLATE_EXPORT ValueHolder {
public:
  typedef V* VPtr;
  typedef const V* ConstVPtr;
  typedef VPtr* VPtrPtr;

  /*! @name Construction, Release, etc
   */
  //@{
  //! initialization with size
  ValueHolder(const Size& s);

  //! initialization with size and physical size
  ValueHolder(const Size& s, const Size& ps);

  //! copy ctor provides full copy!
  ValueHolder(const ValueHolder<V>& h);


  //! assignement provides full copy!
  /*!
    Not exception safe, ie the state of value holder is
    undefined if something goes wrong here. This is on
    purpose to allow memory efficient assignement!
  */
  ValueHolder& operator=(const ValueHolder<V>& h);

  //! free allocated memory upon destruction
  ~ValueHolder();


  //! swap data with another value holder
  void Swap(ValueHolder& vh);
  //@}

  static DataType GetDataType();

  std::size_t  MemSize() const;



  //! return direct r/w access to the value without boundary check
  /*!
    The lookup is based on an integer triplet encapsulated within
    Index.
  */
  V& Value(const Index& i)
  {
  #ifdef USE_ROW_ORDER
    assert(i.w<depth_);
    return data_[i.u*height_depth_ + i.v*depth_ +i.w];
  #else
    assert(i.v<height_);
    return data_[i.w*width_height_ + i.v*height_ +i.u];
  #endif
  }

  //! return direct ro access to the value without boundary check
  /*!
    The lookup is based on an integer triplet encapsulated within
    Index.
  */
  const V& Value(const Index& i) const
  {
  #ifdef USE_ROW_ORDER
    assert(i.w<depth_);
    return data_[i.u*height_depth_ + i.v*depth_ +i.w];
  #else
    assert(i.v<height_);
    return data_[i.w*width_height_ + i.v*height_ +i.u];
  #endif
  }

  //! return pointer to raw data
  VPtr GetData() {return &data_[0];}
  //! const version of GetData()
  ConstVPtr GetData() const {return &data_[0];}


  ConstVPtr GetEnd() const {return &data_[0]+data_.size();}


private:
#ifdef USE_ROW_ORDER
  std::size_t depth_;
  std::size_t height_depth_;
#else
  std::size_t height_;
  std::size_t width_height_;
#endif

  // actual data storage
  std::vector<V> data_;

};

}}} // namespaces

#endif
