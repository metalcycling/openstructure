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

#include "transition.hh"
#include "state.hh"
#include "state_machine.hh"
#include <cassert>

namespace ost { namespace gui {

State::State():
  QObject(),
  mouse_event_transitions_(),
  key_event_transitions_(),
  automatic_transitions_()
{
}
void State::addTransition(SignalTransition* transition)
{
  transition->setParent(this);
}
void State::addTransition(MouseEventTransition* transition)
{
  transition->setParent(this);
  mouse_event_transitions_.append(transition);
}
void State::addTransition(KeyEventTransition* transition)
{
  transition->setParent(this);
  key_event_transitions_.append(transition);
}
void State::addTransition(AutomaticTransition* transition)
{
  transition->setParent(this);
  automatic_transitions_.append(transition);
}
bool State::isActive()
{
  return parent() && dynamic_cast<StateMachine*>(parent())->isActive(this);
}
void State::setActive()
{
  assert(parent());
  dynamic_cast<StateMachine*>(parent())->setActive(this);
}

bool State::checkEvent(QKeyEvent* event)
{
  for(QList<KeyEventTransition*>::iterator it=key_event_transitions_.begin();it!=key_event_transitions_.end();++it){
    std::pair<bool,bool> pair=(*it)->checkEvent(event);
    if(pair.first){
      return pair.second;
    }
  }
  return false;
}

bool State::checkAutomaticTransitions()
{
  for(QList<AutomaticTransition*>::iterator it=automatic_transitions_.begin();it!=automatic_transitions_.end();++it){
    if((*it)->checkTransition()){
      return true;
    }
  }
  return false;
}

bool State::checkEvent(QMouseEvent* event)
{
  for(QList<MouseEventTransition*>::iterator it=mouse_event_transitions_.begin();it!=mouse_event_transitions_.end();++it){
    std::pair<bool,bool> pair=(*it)->checkEvent(event);
    if(pair.first){
      return pair.second;
    }
  }
  return false;
}
void State::onEntry()
{
  emit entered();
}
void State::onExit()
{
  emit exited();
}

}}//ns
