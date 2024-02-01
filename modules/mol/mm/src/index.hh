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

#ifndef OST_MM_INDEX_HH
#define OST_MM_INDEX_HH

/*
  Author: Marco Biasini
 */ 

#include <algorithm>

#include <cstring>

namespace ost { namespace mol{ namespace mm{

namespace impl {
template <uint D>
class IndexBase {
public:
  enum { Dimension = D };
   IndexBase(const IndexBase& rhs) {
    memcpy(data_, rhs.data_, sizeof(uint[D]));
  }
  IndexBase() {
    memset(data_, 0, sizeof(uint[D]));
  }

  IndexBase& operator=(const IndexBase& rhs) {
    memcpy(data_, rhs.data_, sizeof(uint[D]));
    return *this;
  }
  uint operator[](uint idx) const {
    assert(idx < D);
    return data_[idx];
  }
  uint& operator[](uint idx) {
    assert(idx < D);
    return data_[idx];
  } 
  //allows index to be inserted into set
  inline bool operator < (const IndexBase<D>& rhs) const{
    return std::lexicographical_compare(data_, data_+D, rhs.data_, rhs.data_+D);
  } 

  inline bool operator==(const IndexBase<D>& rhs) const{
    return std::equal(data_,data_+D,rhs.data_);
  }

  inline bool operator!=(const IndexBase<D>& rhs) const{
    return !(*this == rhs);
  }
  
private:
  uint data_[D];   
};

} // namespace impl

template <uint D>
class Index;

template <>
class Index<1> : public impl::IndexBase<1> {
public:
  Index() : impl::IndexBase<1>() {}
  Index(uint a) {
    (*this)[0]=a;
  }
  template <typename DS>
  void Serialize(DS& ds){
    ds & (*this)[0];    
  }
};
template <>
class Index<2> : public impl::IndexBase<2> {
public:
  Index() : impl::IndexBase<2>() {}
  Index(uint a, uint b) {
    (*this)[0]=a;
    (*this)[1]=b;    
  }
  template <typename DS>
  void Serialize(DS& ds){
    ds & (*this)[0];    
    ds & (*this)[1];    
  }
};
template <>
class Index<3> : public impl::IndexBase<3> {
public:
  Index() : impl::IndexBase<3>() {}
  Index(uint a, uint b, uint c) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
  }
  template <typename DS>
  void Serialize(DS& ds){
    ds & (*this)[0];    
    ds & (*this)[1];    
    ds & (*this)[2];    
  }
};
template <>
class Index<4> : public impl::IndexBase<4> {
public:
  Index() : impl::IndexBase<4>() {}
  Index(uint a, uint b, uint c, uint d) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
    (*this)[3]=d;            
  }
  template <typename DS>
  void Serialize(DS& ds){
    ds & (*this)[0];    
    ds & (*this)[1];    
    ds & (*this)[2];    
    ds & (*this)[3];    
  }
};
template <>
class Index<5> : public impl::IndexBase<5> {
public:
  Index() : impl::IndexBase<5>() {}
  Index(uint a, uint b, uint c, uint d, uint e) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
    (*this)[3]=d;            
    (*this)[4]=e;                
  }
  template <typename DS>
  void Serialize(DS& ds){
    ds & (*this)[0];    
    ds & (*this)[1];    
    ds & (*this)[2];    
    ds & (*this)[3];    
    ds & (*this)[4];    
  }
};
template <>
class Index<6> : public impl::IndexBase<6> {
public:
  Index() : impl::IndexBase<6>() {}
  Index(uint a, uint b, uint c, uint d, uint e, uint f) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
    (*this)[3]=d;            
    (*this)[4]=e;
    (*this)[5]=f;    
  }
  template <typename DS>
  void Serialize(DS& ds){
    ds & (*this)[0];    
    ds & (*this)[1];    
    ds & (*this)[2];    
    ds & (*this)[3];    
    ds & (*this)[4];    
    ds & (*this)[5];    
  }
};
template <>
class Index<7> : public impl::IndexBase<7> {
public:
  Index() : impl::IndexBase<7>() {}
  Index(uint a, uint b, uint c, uint d, uint e, uint f, uint g) {
    (*this)[0]=a;
    (*this)[1]=b;    
    (*this)[2]=c;        
    (*this)[3]=d;            
    (*this)[4]=e;
    (*this)[5]=f;    
    (*this)[6]=g;        
  }
  template <typename DS>
  void Serialize(DS& ds){
    ds & (*this)[0];    
    ds & (*this)[1];    
    ds & (*this)[2];    
    ds & (*this)[3];    
    ds & (*this)[4];    
    ds & (*this)[5];    
    ds & (*this)[6];    
  }
};
template<uint D>
class IndexIterator {
public:
  typedef Index<D>  IndexType;
  IndexIterator(const IndexType& s, const IndexType& e)
    : start_(s), end_(e), current_(s) {
    
  }  
  
  IndexIterator<D>& operator++() {
    uint current_it=0;
    while (++current_[current_it] > end_[current_it]) {
      current_it++;
      if (current_it < D) {
        current_[current_it-1] = start_[current_it-1];
      } else {
        break;
      }
    }
    return *this;
  }
  const IndexType& operator *() const {
    return current_;
  }
  bool AtEnd() {
    return current_[D-1] > end_[D-1];
  }
private:
  IndexType  start_;
  IndexType  end_;
  IndexType  current_;

};

}}}

#endif
