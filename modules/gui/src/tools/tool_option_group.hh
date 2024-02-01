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
#ifndef OST_GUI_TOOL_OPTION_GROUP_HH
#define OST_GUI_TOOL_OPTION_GROUP_HH

/*
  Author: Marco Biasini
 */


#include <ost/gui/module_config.hh>
#include <ost/gui/tools/tool_option.hh>

#include <QObject>
namespace ost { namespace gui {


class DLLEXPORT_OST_GUI ToolOptionGroup : public QObject {
  Q_OBJECT
public: 
  ToolOptionGroup(const String& group_key, const String& verbose_name);
  
  /// \brief Add int option to default option group
  void AddOption(ToolOption* tool_option);
  
  const String& GetKey() const;
  const String& GetVerboseName() const;
  
  const ToolOptionList& GetOptions() const;
  
  ToolOption* FindOptionByKey(const String& key) const;
private:
  String         key_;
  String         verbose_name_;
  ToolOptionList options_;
};

typedef std::vector<ToolOptionGroup*>  ToolOptionGroups;

}}

#endif
