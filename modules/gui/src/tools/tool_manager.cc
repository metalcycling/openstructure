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
#include "tool_manager.hh"

/*
  Author: Marco Biasini
 */

namespace ost { namespace gui {

ToolManager& ToolManager::Instance()
{
  static ToolManager manager;
  return manager;
}

void ToolManager::AddTool(Tool* tool)
{
  tools_.push_back(tool);
  tool->setParent(this);
  emit this->ToolAdded(tool);
}

void ToolManager::SetActiveTool(Tool* tool)
{
  if (active_tool_==tool) {
    return;
  }
  if (active_tool_) {
    active_tool_->Deactivate();
  }
  active_tool_=tool;
  if (active_tool_) {
    active_tool_->Activate();
  }
  emit this->ActiveToolChanged(tool);
}

Tool* ToolManager::FindToolByName(const String& name) const
{
  for (std::vector<Tool*>::const_iterator i=tools_.begin(), 
       e=tools_.end(); i!=e; ++i) {
    if ((*i)->GetName()==name) {
      return *i;
    }
  }
  return NULL;
}

Tool* ToolManager::GetActiveTool() const
{
  return active_tool_;
}

ToolManager::ToolManager(): active_tool_(NULL) {
  
}

const ToolList& ToolManager::GetTools() const
{
  return tools_;
}

}}
