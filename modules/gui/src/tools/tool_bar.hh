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
#ifndef OST_GUI_TOOLS_TOOL_BAR_HH
#define OST_GUI_TOOLS_TOOL_BAR_HH

/*
  Author: Marco Biasini
 */


#include <ost/gfx/gfx_node.hh>
#include <ost/gui/module_config.hh>
#include <ost/gui/tools/tool.hh>

#include <QToolBar>
#include <QList>

class QActionGroup;

namespace ost { namespace gui {

class DLLEXPORT_OST_GUI ToolBar : public QToolBar {
  Q_OBJECT
public:
  ToolBar(QWidget* parent=NULL);
  
  QAction* AddAction(const QString& text, QIcon& icon);
public slots:
  /// \brief update status of all tools based on the active gfx::GfxNodes
  void OnActiveGfxNodeChanged(gfx::NodePtrList nodes);
  
  void OnToolAdded(const Tool*& tool);
  
  void OnAction(QAction* action);
  
  void ActivateTool(QAction* action);

private:
  QList<QAction*> GetActions();  
  QActionGroup*   group_;
};


}}

#endif
