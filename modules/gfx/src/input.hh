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
#ifndef OST_INPUT_HH
#define OST_INPUT_HH

/*
  constants and other things for input events;
  providing communication between the GUI and the
  GFX module

  Author: Ansgar Philippsen
*/

#include <ost/gfx/module_config.hh>

namespace ost { namespace gfx {

enum InputDevice {
  INPUT_DEVICE_DUMMY=0,
  INPUT_DEVICE_MOUSE
};

enum InputCommand {
  INPUT_COMMAND_NONE=0,
  INPUT_COMMAND_ROTX,
  INPUT_COMMAND_ROTY,
  INPUT_COMMAND_ROTZ,
  INPUT_COMMAND_TRANSX,
  INPUT_COMMAND_TRANSY,
  INPUT_COMMAND_TRANSZ,
  INPUT_COMMAND_SLABN,
  INPUT_COMMAND_SLABF,
  INPUT_COMMAND_AUTOSLAB,
  INPUT_COMMAND_TOGGLE_FOG,
  INPUT_COMMAND_CUSTOM1,
  INPUT_COMMAND_CUSTOM2,
  INPUT_COMMAND_CUSTOM3,
  INPUT_COMMAND_REBUILD
};

enum TransformTarget {
  TRANSFORM_VIEW=0,
  TRANSFORM_RIGID,
  TRANSFORM_TORSION,
  TRANSFORM_ANGLE,
  TRANSFORM_ROTAMER
};

class DLLEXPORT_OST_GFX InputEvent {
public:
  InputEvent(InputDevice dev, InputCommand com, float delta):
    dev_(dev), com_(com), index_(0), target_(0), delta_(delta)
  {}

  InputEvent(InputDevice dev, InputCommand com, int index, 
             int target, float delta):
    dev_(dev), com_(com), index_(index), target_(target), delta_(delta)
  {}

  InputDevice GetDevice() const {return dev_;}
  InputCommand GetCommand() const {return com_;}
  int GetIndex() const {return index_;}
  int GetTarget() const {return target_;}
  float GetDelta() const {return delta_;}
private:
  InputDevice dev_;
  InputCommand com_;
  int index_;
  int target_;
  float delta_;
};

}} // ns

#endif
