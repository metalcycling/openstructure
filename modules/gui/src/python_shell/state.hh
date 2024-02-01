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

#ifndef PYTHON_SHELL_STATE_HH
#define PYTHON_SHELL_STATE_HH

#include <QObject>
#include <QList>

//fw decl
class QKeyEvent;
class QMouseEvent;

namespace ost { namespace gui {

//fw decl
class SignalTransition;
class MouseEventTransition;
class KeyEventTransition;
class AutomaticTransition;

class State: public QObject{
Q_OBJECT
public:
  State();
  void addTransition(SignalTransition * transition);
  void addTransition(MouseEventTransition * transition);
  void addTransition(KeyEventTransition* transition);
  void addTransition(AutomaticTransition* transition);
  bool isActive();
  void setActive();
  bool checkEvent(QKeyEvent* event);
  bool checkEvent(QMouseEvent* event);
  bool checkAutomaticTransitions();
  virtual void onEntry();
  virtual void onExit();
signals:
  void entered();
  void exited();
protected:
  QList<MouseEventTransition*> mouse_event_transitions_;
  QList<KeyEventTransition*> key_event_transitions_;
  QList<AutomaticTransition*> automatic_transitions_;
};


}}//ns
#endif // PYTHON_SHELL_STATE_HH
