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


#include "transition_guard.hh"
#include "python_shell_widget.hh"
#include "shell_history.hh"
#include <iostream>

namespace ost { namespace gui {

TransitionGuard::TransitionGuard():
  QObject()
{
}

bool TransitionGuard::check()
{
  return true;
}

HistoryGuard::HistoryGuard(ShellHistory* history, BlockEditMode mode):
  TransitionGuard(),
  history_(history),
  mode_(mode)
{
}

bool HistoryGuard::check()
{
  return history_->GetCommandMode()==mode_;
}

EditPositionGuard::EditPositionGuard(PythonShellWidget* shell, int flags1, int flags2):
  TransitionGuard(),
  shell_(shell),
  flags1_(flags1),
  flags2_(flags2)
{
}

bool EditPositionGuard::check()
{
    return check_flag(flags1_) &&  check_flag(flags2_);
}

bool EditPositionGuard::check_flag(int flags)
{
  if(flags==0){
    return true;
  }
  bool returnflag=false;
  QTextCursor tc=shell_->textCursor();
  if(flags & SMALLER){
    returnflag|=tc.position()<shell_->GetEditStartBlock().position();
  }
  if(flags & EQUAL){
    returnflag|=tc.position()==shell_->GetEditStartBlock().position();
  }
  if(flags & BIGGER){
    returnflag |= tc.position()> shell_->GetEditStartBlock().position();
  }
  if(flags & ANCHORSMALLER){
    returnflag|=tc.anchor()<shell_->GetEditStartBlock().position();
  }
  if(flags & ANCHOREQUAL){
    returnflag|=tc.anchor()==shell_->GetEditStartBlock().position();
  }
  if(flags & ANCHORBIGGER){
    returnflag |= tc.anchor()> shell_->GetEditStartBlock().position();
  }
  return returnflag;
}

BlockStatusGuard::BlockStatusGuard(PythonShellWidget* shell, int status):
  TransitionGuard(),
  shell_(shell),
  status_(status)
{
}

bool BlockStatusGuard::check()
{
  QString command=shell_->GetCommand();
  return PythonInterpreter::Instance().GetCodeBlockStatus(command) & status_;
}



}} //ns
