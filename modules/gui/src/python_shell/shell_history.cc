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
/*
  Author: Andreas Schenk
 */
#include "shell_history.hh"

namespace ost { namespace gui {


ShellHistory::ShellHistory() :
  history_list_(),
  current_(""),
  current_mode_(EDITMODE_SINGLELINE),
  index_(0)
{
}
void ShellHistory::SetCurrentCommand(const QString& command, BlockEditMode block_edit_mode)
{
  current_=command;
  if(block_edit_mode==EDITMODE_MULTILINE_ACTIVE){
    current_mode_=EDITMODE_MULTILINE_INACTIVE;
  }else{
    current_mode_=block_edit_mode;
  }
}

void ShellHistory::AppendCommand(const QString& command, BlockEditMode block_edit_mode)
{
  if(block_edit_mode==EDITMODE_MULTILINE_ACTIVE){
  history_list_.push_back(std::pair<QString,BlockEditMode>(command,EDITMODE_MULTILINE_INACTIVE));
  }else{
  history_list_.push_back(std::pair<QString,BlockEditMode>(command,block_edit_mode));
  }
  index_=history_list_.size();
  current_="";
  current_mode_=EDITMODE_SINGLELINE;
}

QString ShellHistory::GetCommand()
{
  if(index_>=history_list_.size()){
    return current_;
  }else{
    return history_list_[index_].first;
  }
}

BlockEditMode ShellHistory::GetCommandMode()
{
  if(index_>=history_list_.size()){
    return current_mode_;
  }else{
    return history_list_[index_].second;
  }
}

bool ShellHistory::AtEnd()
{
  return index_>=history_list_.size();
}

void ShellHistory::MoveToPreviousMatch()
{
  int newindex=index_-1;
  while(newindex>=0){
    if(history_list_[newindex].first.indexOf(current_)==0 || current_.size()==0){
      index_=newindex;
      break;
    }
    --newindex;
  }
}

void ShellHistory::MoveToNextMatch()
{
  int newindex=index_+1;
  while(static_cast<unsigned int>(newindex)<history_list_.size()){
    if(history_list_[newindex].first.indexOf(current_)==0 || current_.size()==0){
      index_=newindex;
      break;
    }
    ++newindex;
  }
  if(static_cast<unsigned int>(newindex)!=index_){
    index_=history_list_.size();
  }
}


}}//ns

