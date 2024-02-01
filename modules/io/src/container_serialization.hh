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
#ifndef OST_CONTAINER_SERIALIZATION_HH
#define OST_CONTAINER_SERIALIZATION_HH
/*
  Author: Marco Biasini
 */
#include <ost/stdint.hh>
#include <vector>
#include <map>

#include <boost/type_traits/is_pod.hpp>
namespace ost {  namespace io {

template <typename T>
void Serialize(BinaryDataSink& sink, const std::vector<T>& value) {
  sink << uint32_t(value.size());
  if (boost::is_pod<T>::value) {
    sink.Stream().write(reinterpret_cast<const char*>(&value.front()), 
                        sizeof(T)*value.size());
  } else {
    typename std::vector<T>::const_iterator i=value.begin();
    for (; i!=value.end(); ++i) {
      sink << *i;
    }
  }
}

template <typename T>
void Serialize(BinaryDataSource& source, std::vector<T>& value) {
  uint32_t vector_size;
  source >> vector_size;
  value.resize(vector_size);
  if (boost::is_pod<T>::value) {
    source.Stream().read(reinterpret_cast<char*>(&value.front()), 
                         sizeof(T)*vector_size);
  } else {
    typename std::vector<T>::iterator i=value.begin();
    for (; i!=value.end(); ++i) {
      source >> *i;
    }    
  }
}

template <typename K, typename V>
void Serialize(BinaryDataSink& sink, const std::map<K, V>& value) {
  sink << uint32_t(value.size());
  typename std::map<K, V>::const_iterator i=value.begin();
  for (; i!=value.end(); ++i) {
    sink << i->first;
    sink << i->second;
  }
}

template <typename K, typename V>
void Serialize(BinaryDataSource& source, std::map<K, V>& value) {
  uint32_t vector_size;
  source >> vector_size;
  for (uint32_t i=0; i<vector_size; ++i) {
    K k;
    source >> k;
    std::pair<typename std::map<K, V>::iterator, bool> p;
    p=value.insert(std::make_pair(k, V()));
    source >> p.first->second;
  }
}

}}

#endif
