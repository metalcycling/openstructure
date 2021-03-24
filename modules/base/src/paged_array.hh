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


#ifndef OST_PAGED_ARRAY_HH
#define OST_PAGED_ARRAY_HH

/*
  Author: Marco Biasini, Gabriel Studer
 */

#include <vector>
#include <fstream>
#include <ost/stdint.hh>
#include <ost/message.hh>

namespace ost{ namespace db{


/// \brief Vector style container that splits content in pages, suited for large 
///        amounts of data. Comes with serialization functionality.
template <typename T, uint64_t P>
class PagedArray  {
private:
  typedef std::vector<T> Page;
public:
  typedef T value_type;
  
  PagedArray() {}
  
  T& operator[](uint64_t i) 
  { 
    return pages_[this->GetPage(i)][this->GetOffset(i)]; 
  }
  
  const T& operator[](uint64_t i) const
  { 
    return pages_[this->GetPage(i)][this->GetOffset(i)]; 
  }  
  
  void push_back(const T& t)
  { 
    if (pages_.empty()) {
      pages_.push_back(Page());
      pages_.back().reserve(P);
    }
    if (pages_.back().size()==P) {
      pages_.push_back(Page());
      pages_.back().reserve(P);
    }
    pages_.back().push_back(t);
  }
  uint64_t size() const 
  {
    if (pages_.empty()) {
      return 0;
    } else {
      return static_cast<uint64_t>((pages_.size()-1)*P+pages_.back().size());  
    }
  }

  bool empty() const 
  {
    return pages_.empty();
  }

  void clear() 
  {
    pages_.clear();
  }
  
  // clears everything in [from, to[ and shifts every datapoint after the specified
  // range at the location defined by from
  void ClearRange(uint64_t from, uint64_t to) 
  {
    if(to <= from) return; // invalid range
    uint64_t current_size = this->size();
    if(from >= current_size) return; // nothing to delete
    if(to > current_size) to = current_size; // from is in the valid range, but
                                             // to is too large. Let's just kill
                                             // [from, this->size()[

    uint64_t num_elements_to_shift = current_size - to;
    for(uint64_t i = 0; i < num_elements_to_shift; ++i) {
      (*this)[from + i] = (*this)[to + i];
    }

    uint64_t num_elements_deleted = to - from;
    uint64_t new_size = current_size - num_elements_deleted;
    uint64_t new_last_element_idx = new_size - 1;
    uint64_t new_last_page_idx = this->GetPage(new_last_element_idx);
    pages_.resize(new_last_page_idx + 1);
    uint64_t offset = this->GetOffset(new_last_element_idx);
    pages_.back().resize(offset + 1);
  }
  

  void Write(std::ofstream& out_stream, bool treat_as_pod=true) const
  {
    if (treat_as_pod) {
      uint64_t s=this->size();
      out_stream.write(reinterpret_cast<char*>(&s), sizeof(uint64_t));
      for (typename std::vector<Page>::const_iterator i=pages_.begin(), 
           e=pages_.end(); i!=e; ++i) {
        out_stream.write(reinterpret_cast<const char*>(&i->front()), 
                         i->size()*sizeof(T));        
      }
    } else {
      throw ost::Error("Cannot write non POD paged array!");
    }
  }
  void Read(std::ifstream& in_stream, bool treat_as_pod=true)
  {
    if (treat_as_pod) {
      uint64_t s;
      in_stream.read(reinterpret_cast<char*>(&s), sizeof(uint64_t));
      uint64_t num_pages=s/P+1;
      pages_.resize(num_pages);
      for (uint64_t i=0; i<num_pages; ++i) {
        uint64_t left=std::min(P, s-i*P);
        Page& p=pages_[i];
        p.reserve(P);
        p.resize(left);
        in_stream.read(reinterpret_cast<char*>(&p.front()), 
                       p.size()*sizeof(T));
        if(!in_stream.good()) {
          throw ost::Error("Failed to read data!");
        }
      }
    } else {
      throw ost::Error("Cannot load non POD paged array!");
    }
  }

  // portable serialization
  // (cleanly element by element with fixed-width base-types)
  template <typename DS>
  void Serialize(DS& ds) {
    // assumes that T is not a fundamental type!
    ds & pages_;
  }

  bool operator==(const PagedArray& rhs) const {
    return pages_ == rhs.pages_;
  }
  bool operator!=(const PagedArray& rhs) const {
    return !this->operator==(rhs);
  }
  
private:
  uint64_t GetPage(uint64_t i) const { return i/P; }
  uint64_t GetOffset(uint64_t i) const { return i % P; }
  std::vector<Page> pages_;
};

}} // ns

#endif
