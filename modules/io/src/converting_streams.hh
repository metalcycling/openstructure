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

#ifndef OST_IO_CONVERTING_STREAMS_H
#define OST_IO_CONVERTING_STREAMS_H

#include <iostream>
#include "convert.hh"

#include <ost/stdint.hh>
#include <ost/config.hh>


namespace ost { namespace io {


template<int CONVERSIONTYPE>
class BinaryOStream : public std::basic_ostream<char>
{
public:
  BinaryOStream(std::basic_ostream<char>& ostr):
    std::basic_ostream<char>(ostr.rdbuf())
  {
  }

  BinaryOStream& operator<<(const char& value)
  {
    std::basic_ostream<char>::write(&value,1);
    return *this;
  }
  BinaryOStream& operator<<(const int8_t& value)
  {
    return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const uint8_t& value)
  {
    return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const int16_t& value)
  {
    return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const uint16_t& value)
  {
    return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const int32_t& value)
  {
    return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const uint32_t& value)
  {
    return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const int64_t& value)
  {
    return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const uint64_t& value)
  {
    return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const float& value)
  {
   return write_helper(&value,1);
  }
  BinaryOStream& operator<<(const double& value)
  {
    return write_helper(&value,1);
  }

  BinaryOStream& write(const char* value,std::streamsize n)
  {
    std::basic_ostream<char>::write(value,n);
    return *this;
  }
  BinaryOStream& write(const int8_t* value,std::streamsize n)
  {
    return write_helper(value,n);
  }
  BinaryOStream& write(const uint8_t* value,std::streamsize n)
  {
    return write_helper(value,n);
  }
  BinaryOStream& write(const int16_t* value,std::streamsize n)
  {
    return write_helper(value,n);
  }
  BinaryOStream& write(const uint16_t* value,std::streamsize n)
  {
    return write_helper(value,n);
  }
  BinaryOStream& write(const int32_t* value,std::streamsize n)
  {
    return write_helper(value,n);
  }
  BinaryOStream& write(const uint32_t* value,std::streamsize n)
  {
    return write_helper(value,n);
  }
  BinaryOStream& write(const int64_t* value,std::streamsize n)
  {
    return write_helper(value,n);
  }
  BinaryOStream& write(const uint64_t* value,std::streamsize n)
  {
    return write_helper(value,n);
  }
  BinaryOStream& write(const float* value,std::streamsize n)
  {
   return write_helper(value,n);
  }

  BinaryOStream& write(const double* value,std::streamsize n)
  {
    return write_helper(value,n);
  }

protected:
  template<typename T>
  BinaryOStream& write_helper(const T* value,std::streamsize n)
  {
    for(std::streamsize i=0;i<n;++i){
      T tmp=Convert<CONVERSIONTYPE,T>::To(value[i]);
      std::basic_ostream<char>::write(reinterpret_cast<char*>(&tmp),sizeof(T));
    }
    return *this;
  }

};

template<int CONVERSIONTYPE>
class BinaryIStream : public std::basic_istream<char>
{
public:
  BinaryIStream( std::basic_istream<char>& istr):
    std::basic_istream<char>(istr.rdbuf())
  {
  }

  BinaryIStream& operator>>(int8_t& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(uint8_t& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(int16_t& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(uint16_t& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(int32_t& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(uint32_t& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(int64_t& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(uint64_t& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(float& value)
  {
    return read_helper(&value,1);
  }
  BinaryIStream& operator>>(double& value)
  {
    return read_helper(&value,1);
  }

  BinaryIStream& read(char* value,std::streamsize n)
  {
    std::basic_istream<char>::read(value,n);
    return *this;
  }
  BinaryIStream& read(uint8_t* value,std::streamsize n)
  {
    return read_helper(value,n);
  }
  BinaryIStream& read(int8_t* value,std::streamsize n)
  {
    return read_helper(value,n);
  }
  BinaryIStream& read(int16_t* value,std::streamsize n)
  {
    return read_helper(value,n);
  }
  BinaryIStream& read(uint16_t* value,std::streamsize n)
  {
    return read_helper(value,n);
  }
  BinaryIStream& read(int32_t* value,std::streamsize n)
  {
    return read_helper(value,n);
  }
  BinaryIStream& read(uint32_t* value,std::streamsize n)
  {
    return read_helper(value,n);
  }
  BinaryIStream& read(int64_t* value,std::streamsize n)
  {
    return read_helper(value,n);
  }
  BinaryIStream& read(uint64_t* value,std::streamsize n)
  {
    return read_helper(value,n);
  }
  BinaryIStream& read(float* value,std::streamsize n)
  {
   return read_helper(value,n);
  }
  BinaryIStream& read(double* value,std::streamsize n)
  {
    return read_helper(value,n);
  }

protected:
  template<typename T>
  BinaryIStream& read_helper(T* value,std::streamsize n)
  {
    std::basic_istream<char>::read(reinterpret_cast<char*>(value),n*sizeof(T));
    for(std::streamsize i=0;i<n;++i){
      Convert<CONVERSIONTYPE,T>::FromIP(&value[i]);
    }
    return *this;
  }
};



}  // io ns



} // ost ns


#endif


