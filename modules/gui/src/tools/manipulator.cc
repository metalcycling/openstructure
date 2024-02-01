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
#include <ost/gui/gosty_app.hh>

#include <ost/platform.hh>

#include <ost/gfx/scene.hh>
#include <ost/gfx/gfx_object.hh>
#include <ost/gfx/gfx_node.hh>

#include <ost/gui/scene_selection.hh>

#include "manipulator.hh"

#include <QDir>
namespace ost { namespace gui {

Manipulator::Manipulator():
  Tool("Rigid Body Manipulator")
{

}

void Manipulator::MouseMove(const MouseEvent& event)
{
  static bool warn=true;
  // reference is GLCanvas::HandleMouseMoveEvent
  if (event.GetButtons() & MouseEvent::LeftButton) {
    if (event.GetButtons() & MouseEvent::MiddleButton) {
      if (event.IsShiftPressed()) {
        // slab in view, ignore here
      } else {
        Real d=static_cast<Real>(event.GetDelta().x()+event.GetDelta().y());
        SendCommand(gfx::INPUT_COMMAND_ROTZ, -d);
      }
    } else {
      if (event.IsShiftPressed()) {
        SendCommand(gfx::INPUT_COMMAND_TRANSX,
                    static_cast<Real>(event.GetDelta().x()));
        SendCommand(gfx::INPUT_COMMAND_TRANSY,
                    static_cast<Real>(-event.GetDelta().y()));
      } else {
        SendCommand(gfx::INPUT_COMMAND_ROTY,
                    static_cast<Real>(event.GetDelta().x())); // rotation around y
        SendCommand(gfx::INPUT_COMMAND_ROTX,
                    static_cast<Real>(event.GetDelta().y())); // rotation around x
      }
    }
  } else if (event.GetButtons() & MouseEvent::MiddleButton) {
    if (event.IsShiftPressed()) {
      if(warn) {
        LOG_WARNING("Manipulation Tool: MMB + Shift for z rotation is deprecated; use LMB+MMB instead");
        warn=false;
      }
      // rotation around z
      // WHY? this does not match the scene view transformation, which has a slab tf here
      // will be removed in a next update...
      Real d=static_cast<Real>(event.GetDelta().x()+event.GetDelta().y());
      SendCommand(gfx::INPUT_COMMAND_ROTZ, d);        
    } else {
      SendCommand(gfx::INPUT_COMMAND_TRANSZ,
                  static_cast<Real>(event.GetDelta().y())); // translate along z
    }
  }  
}

void Manipulator::DoubleClick(const MouseEvent& event)
{
  // todo pick object and set it as input grabber...
}

void Manipulator::SendCommand(gfx::InputCommand command, Real delta)
{
  gfx::InputEvent event(gfx::INPUT_DEVICE_MOUSE, command, delta);

  int active_node_count = SceneSelection::Instance()->GetActiveNodeCount();
  if(active_node_count > 0){
	for(int i = 0; i<active_node_count;i++){
      gfx::GfxNodeP node=SceneSelection::Instance()->GetActiveNode(i);
      if(node){
        gfx::GfxObj* obj = dynamic_cast<gfx::GfxObj*> (node.get());
          if (obj) {
            obj->OnInput(event);
          }
      }
	}
  }
}

bool Manipulator::CanOperateOn(gfx::NodePtrList nodes)
{
  for(unsigned int i=0; i < nodes.size(); i++){
    if(!dynamic_cast<gfx::GfxObj*>(nodes[i].get())){
      return false;
    }
  }
  return true;
}

QIcon Manipulator::GetIcon() const
{
  QDir icon_path(GetSharedDataPath().c_str());
  icon_path.cd("gui");
  icon_path.cd("icons");
  QIcon icon = QIcon(icon_path.absolutePath()+QDir::separator()+QString("rigid_body_manipulator_icon.png"));
  return icon;
}

}}
