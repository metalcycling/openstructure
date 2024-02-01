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


#include <ost/db/binary_container.hh>
#include <ost/message.hh>


namespace ost{ namespace db{

LinearCharacterContainerPtr LinearCharacterContainer::Load(const String& fn) {

  std::ifstream in_stream(fn.c_str(), std::ios::binary);
  if (!in_stream) {
    throw ost::Error("Could not open " + fn);
  }
  LinearCharacterContainerPtr loaded_container(new LinearCharacterContainer);
  loaded_container->data_.Read(in_stream);
  loaded_container->num_elements_ = loaded_container->data_.size();
  return loaded_container;
}

void LinearCharacterContainer::Save(const String& filename) const {

  std::ofstream out_stream(filename.c_str(), std::ios::binary);
  if (!out_stream) {
    throw ost::Error("Could not open " + filename);
  }

  data_.Write(out_stream);
}

void LinearCharacterContainer::AddCharacters(const String& characters) {
  for(uint i = 0; i < characters.size(); ++i) {
    data_.push_back(characters[i]);
  }
  num_elements_ = data_.size();
}

void LinearCharacterContainer::ClearRange(std::pair<uint64_t, uint64_t> range) {

  if(range.second <= range.first) {
    throw ost::Error("Cannot clear invalid range in "
                     "LinearCharacterContainer!");
  } 

  if(range.first >= num_elements_  || range.second > num_elements_) {
    throw ost::Error("Cannot clear invalid range in "
                     "LinearCharacterContainer!");
  }

  data_.ClearRange(range.first, range.second);
  num_elements_ = data_.size();
}

char LinearCharacterContainer::GetCharacter(uint64_t idx) const {

  if(idx >= num_elements_) {
    throw ost::Error("Cannot access invalid position in "
                     "LinearCharacterContainer!");
  }

  return data_[idx];
}

void LinearCharacterContainer::GetCharacters(std::pair<uint64_t, uint64_t> range,
                                             String& s) const {

  if(range.second <= range.first) {
    throw ost::Error("Cannot access invalid position in "
                     "LinearCharacterContainer!");
  } 

  if(range.first >= num_elements_ || range.second > num_elements_) {
    throw ost::Error("Cannot access invalid position in "
                     "LinearCharacterContainer!");
  }

  uint64_t size = range.second - range.first;
  s = String(size, '-');
  for(uint64_t i = 0; i < size; ++i) {
    s[i] = data_[range.first + i];
  }
}


LinearPositionContainerPtr LinearPositionContainer::Load(const String& fn) {

  std::ifstream in_stream(fn.c_str(), std::ios::binary);
  if (!in_stream) {
    throw ost::Error("Could not open " + fn);
  }
  LinearPositionContainerPtr loaded_container(new LinearPositionContainer);
  loaded_container->data_.Read(in_stream);
  loaded_container->num_elements_ = loaded_container->data_.size();
  return loaded_container;
}

void LinearPositionContainer::Save(const String& filename) const {

  std::ofstream out_stream(filename.c_str(), std::ios::binary);
  if (!out_stream) {
    throw ost::Error("Could not open " + filename);
  }

  data_.Write(out_stream);
}


void LinearPositionContainer::AddPositions(const geom::Vec3List& positions) {

  for(uint i = 0; i < positions.size(); ++i) {

    int64_t data = 0;
    bool low_accuracy = false;

    int64_t int_x_pos = std::floor(positions[i][0] * Real(100) + Real(0.5));
    int64_t int_y_pos = std::floor(positions[i][1] * Real(100) + Real(0.5));
    int64_t int_z_pos = std::floor(positions[i][2] * Real(100) + Real(0.5));

    if(int_x_pos < CUSTOM_INT_MIN || int_x_pos > CUSTOM_INT_MAX ||
       int_y_pos < CUSTOM_INT_MIN || int_y_pos > CUSTOM_INT_MAX ||
       int_z_pos < CUSTOM_INT_MIN || int_z_pos > CUSTOM_INT_MAX) {

      // let's try to bring it into the allowed range using a factor of
      // ten (results in only one digit of accuracy)
      int_x_pos = std::floor(positions[i][0] * Real(10) + Real(0.5));
      int_y_pos = std::floor(positions[i][1] * Real(10) + Real(0.5));
      int_z_pos = std::floor(positions[i][2] * Real(10) + Real(0.5));

      if(int_x_pos < CUSTOM_INT_MIN || int_x_pos > CUSTOM_INT_MAX ||
         int_y_pos < CUSTOM_INT_MIN || int_y_pos > CUSTOM_INT_MAX ||
         int_z_pos < CUSTOM_INT_MIN || int_z_pos > CUSTOM_INT_MAX) {
         // we don't go lower in accuracy, let's just throw an error
         throw ost::Error("Out of bound position observed!");
      }

      low_accuracy = true;
    }

    int_x_pos += CUSTOM_INT_MAX;
    int_y_pos += CUSTOM_INT_MAX;
    int_z_pos += CUSTOM_INT_MAX;

    data += (int_x_pos << X_POS_SHIFT);
    data += (int_y_pos << Y_POS_SHIFT);
    data += (int_z_pos << Z_POS_SHIFT);

    // as it is constructed, data is guaranteed to be positive and we only use
    // 63 bits.
    // Let's misuse the sign bit as a flag, whether we have a low accuracy
    // representation
    if(low_accuracy) {
      data += (static_cast<int64_t>(1) << 63);
    }

    data_.push_back(data);
  }

  num_elements_ = data_.size();
}


void LinearPositionContainer::ClearRange(std::pair<uint64_t, uint64_t> range) {

  if(range.second <= range.first) {
    throw ost::Error("Cannot clear invalid range in "
                     "LinearPositionContainer!");
  } 

  if(range.first >= num_elements_  || range.second > num_elements_) {
    throw ost::Error("Cannot clear invalid range in "
                     "LinearPositionContainer!");
  }

  data_.ClearRange(range.first, range.second);
  num_elements_ = data_.size();
}


void LinearPositionContainer::GetPosition(uint64_t idx, geom::Vec3& pos) const {

  if(idx >= num_elements_) {
    throw ost::Error("Cannot access invalid position in "
                     "LinearPositionContainer!");
  }

  int64_t data = data_[idx];
  Real factor = (data < 0) ? 0.1 : 0.01;
  pos[0] = factor * Real(static_cast<int>(((data & X_POS_MASK) >> X_POS_SHIFT)) - 
                         CUSTOM_INT_MAX);
  pos[1] = factor * Real(static_cast<int>(((data & Y_POS_MASK) >> Y_POS_SHIFT)) - 
                         CUSTOM_INT_MAX);
  pos[2] = factor * Real(static_cast<int>(((data & Z_POS_MASK) >> Z_POS_SHIFT)) - 
                         CUSTOM_INT_MAX);
}

void LinearPositionContainer::GetPositions(std::pair<uint64_t, uint64_t> range, 
                                           geom::Vec3List& positions) const {

  if(range.second <= range.first) {
    throw ost::Error("Cannot access invalid position in "
                     "LinearPositionContainer!");
  } 

  if(range.first >= num_elements_ || range.second > num_elements_) {
    throw ost::Error("Cannot access invalid position in "
                     "LinearPositionContainer!");
  }

  int num_positions = range.second - range.first;
  positions.resize(num_positions);

  for(int i = 0; i < num_positions; ++i) {
    this->GetPosition(range.first + i, positions[i]);
  }
}

}} //ns
