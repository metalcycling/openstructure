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
#include "tool_option.hh"

namespace ost { namespace gui {

template class TEMPLATE_DEF_EXPORT ToolOptionNum<int,ToolOption::INT>;
template class TEMPLATE_DEF_EXPORT ToolOptionNum<float,ToolOption::FLOAT>;

ToolOption::ToolOption(const String& key, const String& verbose_name, 
                       Type type):
  key_(key), verbose_name_(verbose_name), type_(type)
{
  
}

const String& ToolOption::GetKey() const
{
  return key_;
}

const String& ToolOption::GetVerboseName() const 
{
  return verbose_name_;
}

ToolOption::Type ToolOption::GetType() const 
{
  return type_;
}


ToolOptionEnum::ToolOptionEnum(const String& key, const String& verbose_name):
  ToolOption(key, verbose_name, ToolOption::ENUM), index_(-1)
{
  
}

void ToolOptionEnum::Add(const String& text, int tag)
{
  tuples_.push_back(Tuple());
  tuples_.back().text=text;
  tuples_.back().tag=tag;
}

ToolOptionButton::ToolOptionButton(const String& key,
                                   const String& verbose_name,
                                   QObject* receiver,
                                   const char *slot_method):
    ToolOption(key, verbose_name, ToolOption::BUTTON), slot_method_(slot_method), receiver_(receiver)
{

}

}}
