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
#ifndef OST_GUI_MOUSE_EVENT_HH
#define OST_GUI_MOUSE_EVENT_HH

#include <ost/gui/module_config.hh>

#include <QPoint>
#include <QMouseEvent>

namespace ost { namespace gui {


class GLCanvas;

class DLLEXPORT_OST_GUI MouseEvent {
public:
  typedef enum {
    LeftButton=Qt::LeftButton,
    RightButton=Qt::RightButton,
    MiddleButton=Qt::MiddleButton
  } Button;
  typedef int Buttons;
public:
  MouseEvent(Buttons buttons, const QPoint& from,
             const QPoint& to, Qt::KeyboardModifiers modifiers);
  
  /// \brief direction of mouse movement
  QPoint GetDelta() const;
  QPoint GetLastPos() const;
  QPoint GetPos() const;
  Buttons GetButtons() const;
  bool IsShiftPressed() const;
  bool IsMetaPressed() const;
  bool IsControlPressed() const;
private:
  Buttons buttons_;  
  QPoint  from_;
  QPoint  to_;
  Qt::KeyboardModifiers modifiers_;
};

}}

#endif
