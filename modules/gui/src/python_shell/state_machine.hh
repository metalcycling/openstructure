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

#ifndef PYTHON_SHELL_STATE_MACHINE_HH
#define PYTHON_SHELL_STATE_MACHINE_HH

#include <QObject>

//fw decl
class QEvent;

namespace ost { namespace gui {

//fw decl
class State;

class StateMachine: public QObject{
Q_OBJECT
public:
  StateMachine(QObject* parent);
  void addState(State* state);
  void setInitialState(State* state);
  void start();
  bool isActive(State* state);
  void setActive(State* state);
  bool eventFilter (QObject * watched, QEvent * event);
protected:
  State* active_state_;
};

}}//ns

#endif // PYTHON_SHELL_STATE_MACHINE_HH
