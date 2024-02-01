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
#ifndef OST_IO_BINARY_DATA_SOURCE_HH
#define OST_IO_BINARY_DATA_SOURCE_HH
/*
  Author: Marco Biasini
  
  minimalistic binary serialization interface
 */

#include <iostream>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>

#include <ost/io/module_config.hh>
namespace ost { namespace io {

class BinaryDataSource {
public:
    BinaryDataSource(std::istream& stream) 
      : stream_(stream) {
    }
    template <typename SERIALIZABLE>
    BinaryDataSource& operator >> (SERIALIZABLE& object) {
      Serialize(*this, object);
      return *this;
    }
    template <typename SERIALIZABLE>
    BinaryDataSource& operator & (SERIALIZABLE& object) {
      return this->operator>>(object);
    }
    std::istream& Stream() {
      return stream_;
    }
    
    bool IsSource() { return true; }
private:
  std::istream& stream_;
};

namespace detail {

template <bool B, typename T>
struct SerializeHelper;

template <typename T>
struct SerializeHelper<false, T> {
  SerializeHelper(BinaryDataSource& source, T& value)
  {
    value.Serialize(source);
  }
};

template <typename T>
struct SerializeHelper<true, T> {
  SerializeHelper(BinaryDataSource& source, T& value)
  {
    source.Stream().read(reinterpret_cast<char*>(&value), sizeof(T));
  }
};
}

template <typename T>
void Serialize(BinaryDataSource& source, T& value)
{  
  detail::SerializeHelper<boost::is_integral<T>::value || 
                          boost::is_floating_point<T>::value,
                          T>(source, value);
}

inline void RawSerialize(BinaryDataSource& sink, 
                         char* value, size_t size)
{
    sink.Stream().read(value, size);
}

inline void Serialize(BinaryDataSource& source, String& str)
{
  static char tmp_buf[1024];
  size_t str_len;
  source.Stream().read(reinterpret_cast<char*>(&str_len), sizeof(size_t));
  source.Stream().read(tmp_buf, str_len);
  str.reserve(str_len+1);
  str.assign(tmp_buf, str_len); 
}
}}

#endif
