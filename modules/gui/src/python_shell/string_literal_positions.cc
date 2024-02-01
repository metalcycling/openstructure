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
#include "string_literal_positions.hh"

namespace ost{namespace gui{

StringLiteralPositions::StringLiteralPositions():
  QTextBlockUserData(),
  start_positions_(),
  end_positions_()
{
}

StringLiteralPositions::~StringLiteralPositions()
{
}

int StringLiteralPositions::GetStringLiteralIndex(int position)
{
  unsigned int i=0;
  for(;i<start_positions_.size();++i){
    if(position<start_positions_[i]){
      break;
    }
  }
  if(i==0 ){
    return -1;
  }else{
    if(position<=end_positions_[i-1]){
      return i-1;
    }else{
      return -1;
    }
  }
}

void StringLiteralPositions::Append(int start, int end)
{
  start_positions_.push_back(start);
  end_positions_.push_back(end);
}
int StringLiteralPositions::GetStart(int index)
{
  return start_positions_[index];
}
int StringLiteralPositions::GetEnd(int index)
{
  return end_positions_[index];
}

}}//ns
