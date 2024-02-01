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
#ifndef OST_IO_BINARY_DATA_SINK_HH
#define OST_IO_BINARY_DATA_SINK_HH
/*
  Author: Marco Biasini
  
  minimalistic binary serialization interface
 */

#include <iostream>

#include <ost/io/module_config.hh>

namespace ost { namespace io {

class BinaryDataSink {
public:
    BinaryDataSink(std::ostream& stream) 
      : stream_(stream) {
    }
    template <typename SERIALIZABLE>
    BinaryDataSink& operator << (SERIALIZABLE object) {
      Serialize(*this, object);
      return *this;
    }
    template <typename SERIALIZABLE>
    BinaryDataSink& operator & (SERIALIZABLE object) {
      return this->operator<<(object);
    }
    std::ostream& Stream() {
      return stream_;
    }
    
    bool IsSource() { return false; }
private:
  std::ostream& stream_;
};

namespace detail {

template <bool B, typename T>
struct SerializeWriteHelper;

template <typename T>
struct SerializeWriteHelper<false, T> {
  SerializeWriteHelper(BinaryDataSink& sink, T& value)
  {
    value.Serialize(sink);
  }
};

template <typename T>
struct SerializeWriteHelper<true, T> {
  SerializeWriteHelper(BinaryDataSink& sink, T& value)
  {
    sink.Stream().write(reinterpret_cast<char*>(&value), sizeof(T));
  }
};
}

template <typename T>
void Serialize(BinaryDataSink& sink, const T& value)
{  
  detail::SerializeWriteHelper<boost::is_integral<T>::value || 
                               boost::is_floating_point<T>::value,
                               T>(sink, const_cast<T&>(value));
}


inline void RawSerialize(BinaryDataSink& sink, char* value, size_t size)
{
    sink.Stream().write(reinterpret_cast<char*>(value), size);
}

inline void Serialize(BinaryDataSink& sink, const String& str) {
  size_t str_len=str.length();
  sink.Stream().write(reinterpret_cast<char*>(&str_len), sizeof(size_t));
  sink.Stream().write(str.c_str(), str.length());
}
}}

#endif
