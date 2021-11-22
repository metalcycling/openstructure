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


#ifndef OST_DB_LINEAR_BINARY_CONTAINER_HH
#define OST_DB_LINEAR_BINARY_CONTAINER_HH

#include <vector>
#include <boost/shared_ptr.hpp>
#include <ost/base.hh>
#include <ost/stdint.hh>
#include <ost/geom/vec3.hh>
#include <ost/paged_array.hh>


namespace ost { namespace db {


class LinearCharacterContainer;
class LinearPositionContainer;
typedef boost::shared_ptr<LinearCharacterContainer> LinearCharacterContainerPtr;
typedef boost::shared_ptr<LinearPositionContainer> LinearPositionContainerPtr;


class LinearCharacterContainer {

public:

  LinearCharacterContainer(): num_elements_(0) { }

  static LinearCharacterContainerPtr Load(const std::string& filename);

  void Save(const std::string& filename) const;

  void AddCharacters(const std::string& characters);

  void ClearRange(std::pair<uint64_t, uint64_t> range);

  char GetCharacter(uint64_t idx) const;

  void GetCharacters(std::pair<uint64_t, uint64_t> range, String& s) const;

  uint64_t GetNumElements() const { return num_elements_; }

  bool operator==(const LinearCharacterContainer& other) const {
    return data_ == other.data_; // no check for num_elements_
                                 // they are coupled anyway
  }

  bool operator!=(const LinearCharacterContainer& other) const {
    return !(*this == other);
  }

private:

  PagedArray<char, 1048576> data_;
  // we track number of entries manually for fast checking of indices
  // the .size() function of data_ is too slow
  uint64_t num_elements_;
};




class LinearPositionContainer {

public:

  LinearPositionContainer(): num_elements_(0) { }

  static LinearPositionContainerPtr Load(const std::string& filename);

  void Save(const std::string& filename) const;

  void AddPositions(const geom::Vec3List& positions);

  void ClearRange(std::pair<uint64_t, uint64_t> range);

  void GetPosition(uint64_t idx, geom::Vec3& pos) const;

  void GetPositions(std::pair<uint64_t, uint64_t> range, 
                    geom::Vec3List& positions) const;

  uint64_t GetNumElements() const { return num_elements_; }

  bool operator==(const LinearPositionContainer& other) const {
    return data_ == other.data_; // no check for num_elements_
                                 // they are coupled anyway
  }

  bool operator!=(const LinearPositionContainer& other) const {
    return !(*this == other);
  } 

private:

  // 21 bits starting at index 1 are set
  static const uint64_t X_POS_MASK = 9223367638808264704; 
  static const int X_POS_SHIFT = 42;

  // 21 bits starting at index 22 are set
  static const uint64_t Y_POS_MASK = 4398044413952;
  static const int Y_POS_SHIFT = 21;

  // 21 bits starting at index 43 (the last 21 bits) are set
  static const uint64_t Z_POS_MASK = 2097151;
  static const int Z_POS_SHIFT = 0;

  // range that can be represented as 21 bit integer
  static const int CUSTOM_INT_MIN = -1048576;
  static const int CUSTOM_INT_MAX = 1048576;

  PagedArray<uint64_t, 1048576> data_;
  // we track number of entries manually for fast checking of indices
  // the .size() function of data_ is too slow
  uint64_t num_elements_;  
};


}} //ns

#endif
