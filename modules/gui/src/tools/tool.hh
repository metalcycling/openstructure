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
#ifndef OST_GUI_TOOL_HH
#define OST_GUI_TOOL_HH

/*
  Author: Marco Biasini
 */
#include <vector>


#include <ost/gui/module_config.hh>
#include <ost/gfx/gfx_node.hh>
#include <ost/gui/tools/mouse_event.hh>
#include <ost/gui/tools/tool_options.hh>

#include <QIcon>

namespace ost { namespace gui {

class ToolOptions;
class Tool;
class GLCanvas;

/// \brief  base class for tools
/// 
/// Parameters of the tool that are configurable should be implemented using the
/// ToolOptions facility. This automatically takes care of remembering the
/// used options over different sessions (in the future) and creates a widget 
/// to modify these options.
/// 
/// \section input_events Input Events
/// 
/// Tools receive input events when the Control key (Command on macOS) is
/// pressed. This includes both keyboard events and mouse events such as mouse
/// move events, clicks and Real clicks.
/// 
/// Tools only receive mouse move events when one of the mouse buttons is 
/// pressed. 
/// 
/// #CanOperateOn() determines if the tool can manipulate a certain 
/// gfx::GfxNodeP. This method is used to enabled/disable certain actions (in
/// the toolbar, for example).
/// 
/// \sa \ref color_tool.py "Color Tool Example"
class DLLEXPORT_OST_GUI Tool : public QObject {
  Q_OBJECT
protected:
  Tool(const String& name);
public:
  virtual ~Tool();
  /// \brief left mouse click
  virtual void Click(const MouseEvent& event) { }
  
  /// \brief Real click
  virtual void DoubleClick(const MouseEvent& event) { }

  virtual void MousePress(const MouseEvent& event) { }
  
  virtual void MouseMove(const MouseEvent& event) { }
  
  virtual void RenderGL() { }

  
  /// \brief called when the tool becomes active
  virtual void Activate() { }
  
  /// \brief called when the tool becomes inactive
  virtual void Deactivate() { }
  
  /// \brief whether the tool is capable of manipulating the given node
  virtual bool CanOperateOn(gfx::NodePtrList nodes)=0;
  ToolOptions* GetToolOptions() const;
  
  const String& GetName() const;

  virtual QIcon GetIcon() const;

private:
  String       name_;
  ToolOptions* tool_options_;
}; 


/// \example color_tool.py
/// 
/// Shows how to write a tool in Python, register it and interact with the 
/// graphical objects in the scene.

}}

#endif
