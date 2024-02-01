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
#ifndef OST_GUI_TOOL_MANAGER_HH
#define OST_GUI_TOOL_MANAGER_HH

#include <ost/gui/module_config.hh>

#include <ost/gui/tools/tool.hh>

#include <QObject>
/*
  Author: Marco Biasini
 */

namespace ost { namespace gui {

typedef std::vector<Tool*> ToolList;

class DLLEXPORT_OST_GUI ToolManager : public QObject {
  Q_OBJECT
public:
  static ToolManager& Instance();
  /// \brief add tool to tool manager
  /// 
  /// The ToolManager takes ownership of the tool.
  void AddTool(Tool* tool);
  
  void SetActiveTool(Tool* tool);
  Tool* FindToolByName(const String& name) const;
  Tool* GetActiveTool() const;
  
  const ToolList& GetTools() const;
signals:
  void ToolAdded(Tool*);
  void ActiveToolChanged(Tool*);
private:
  ToolManager();
  ToolList   tools_;
  Tool*      active_tool_;
};

}}

#endif
