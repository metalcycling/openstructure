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
#ifndef OST_GUI_TOOL_OPTIONS_HH
#define OST_GUI_TOOL_OPTIONS_HH

/*
  Author: Marco Biasini
 */
 

#include <ost/gui/module_config.hh>
#include <ost/gui/tools/tool_option_group.hh>

#include <QObject>
namespace ost { namespace gui {


/// \brief  tool options
/// 
/// Most \ref Tool "tools" are configurable. To simplify and standardise the
/// process of providing widgets to configure these tools, these options are
/// defined using an  instance of ToolOptions. A widget to modify these options
/// is then  automaticallly created when the tool is added to the ToolManager.
class DLLEXPORT_OST_GUI ToolOptions : public QObject {
  Q_OBJECT
public:
  ToolOptions();
  /// \brief Add option to default option group
  void AddOption(ToolOption* option);
  
  /// \brief Add option to group
  void AddOption(const String& group, ToolOption* option);
  
  ToolOptionGroup* FindGroupByName(const String& group_name) const;
  
  ToolOption* GetOption(const String& group_name, const String& key) const;
  ToolOption* GetOption(const String& key) const;
  const ToolOptionGroups& GetGroups() const;
private:
  ToolOptionGroups groups_;
};

}}

#endif
