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
/*
  Author: Marco Biasini
 */

#include "tool_options.hh"


namespace ost { namespace gui {

ToolOptions::ToolOptions()
{
  groups_.push_back(new ToolOptionGroup("", ""));
}

void ToolOptions::AddOption(ToolOption* option)
{
  this->FindGroupByName("")->AddOption(option);
}

ToolOptionGroup* ToolOptions::FindGroupByName(const String& group_name) const
{
  for (ToolOptionGroups::const_iterator i=groups_.begin(), 
       e=groups_.end(); i!=e; ++i) {
    if ((*i)->GetKey()==group_name) {
      return *i;
    }
  }
  return NULL;
}

const ToolOptionGroups& ToolOptions::GetGroups() const
{
  return groups_;
}

ToolOption* ToolOptions::GetOption(const String& group_name, 
                                     const String& key) const
{
  ToolOptionGroup* group=this->FindGroupByName(group_name);
  return group ? group->FindOptionByKey(key) : NULL;
}

ToolOption* ToolOptions::GetOption(const String& key) const
{
  return this->GetOption("", key);
}

void ToolOptions::AddOption(const String& group_name, 
                            ToolOption* option)
{
  ToolOptionGroup* group=this->FindGroupByName(group_name);
  if (!group) {
    return;
  }
  group->AddOption(option);
}

}}
