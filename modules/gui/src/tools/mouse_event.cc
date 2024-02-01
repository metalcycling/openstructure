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
#include "mouse_event.hh"

namespace ost { namespace gui {

MouseEvent::MouseEvent(Buttons buttons, const QPoint& from,
                       const QPoint& to, Qt::KeyboardModifiers mod)
  : buttons_(buttons), from_(from), to_(to), modifiers_(mod)
{  }

MouseEvent::Buttons MouseEvent::GetButtons() const
{
  return buttons_;
}

QPoint MouseEvent::GetDelta() const
{
  return to_-from_;
}

QPoint MouseEvent::GetLastPos() const
{
  return from_;
}

QPoint MouseEvent::GetPos() const
{
  return to_;
}

bool MouseEvent::IsShiftPressed() const
{
  return modifiers_ & Qt::ShiftModifier;
}

bool MouseEvent::IsMetaPressed() const
{
  return modifiers_ & Qt::MetaModifier;
}

bool MouseEvent::IsControlPressed() const
{
  return modifiers_ & Qt::ControlModifier;
}

}}

