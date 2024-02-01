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
/*  Author: Marco Biasini
*/

#include "python_context_parser.hh"

namespace ost { namespace gui {
  

void PythonContext::Clear()
{
  for (size_t i = 0; i < sub_contexts_.size(); ++i)
    delete sub_contexts_[i];
  sub_contexts_.clear();
  token_list_.clear();
  range_ = Range(String::npos, 0);
}

bool PythonContext::FindLastToken(PythonToken::Type type,
                                  const QString& value, 
                                  PythonToken& token)
{
  std::vector<PythonToken>::reverse_iterator i = token_list_.rbegin();
  for( ; i != token_list_.rend(); ++i ) {
    const PythonToken& t = *i;
    if (t.GetType() == type && (value == "" || t.GetValue() == value)) {
      token = t;
      return true;
    }
  }
  return false;
}

bool PythonContext::FindFirstToken(PythonToken::Type type, 
                                   const QString& value, 
                                   PythonToken& token)
{
  std::vector<PythonToken>::iterator i = token_list_.begin();
  for( ; i != token_list_.end(); ++i ) {
    const PythonToken& t = *i;
    if (t.GetType() == type && (value == "" || t.GetValue() == value)) {
      token = t;
      return true;
    }
  }
  return false;
}

PythonContextParser::PythonContextParser() 
  : command_(""),tokenizer_(),root_(new PythonContext)
{
  root_->SetRange(Range(0,0));
}

void PythonContextParser::Parse(const QString& command)
{
  command_ = command;
  tokenizer_.SetCommand(command_);
  // Remove old tree structure before building a new one.
  root_->Clear();
  root_->SetRange(Range(0,0));
  this->Parse(root_);
  this->root_->SetRange(Range(0,command_.size()));    
}

void PythonContextParser::Parse(PythonContext*& parent_context)
{
  PythonToken t;
  PythonContext* current_context = new PythonContext;
  current_context->SetRange(Range(parent_context->GetRange().End(),0));
  while ((t = tokenizer_.NextToken()).GetType() != PythonToken::END) {
    if (t.GetType() == PythonToken::GROUPING) {
      if (t.GetValue() == "(" || t.GetValue() == "[" || t.GetValue() == "{") {
        current_context->AddToken(t);        
        Parse(current_context);
        if (tokenizer_.CurrentToken().GetType() != PythonToken::END)
          current_context->AddToken(tokenizer_.CurrentToken());   
        else if (! current_context->GetSubContexts().empty()) {
          Range last_sub = current_context->GetSubContexts().back()->GetRange();
          Range curr_range(current_context->GetRange().location,
                           last_sub.End()-current_context->GetRange().location);
          current_context->SetRange(curr_range);          
        }

      }
      if (t.GetValue() == "," || t.GetValue() == ";") {
        parent_context->AddContext(current_context);
        current_context = new PythonContext;
        current_context->SetRange(Range(t.GetRange().End(),0));
      }
      // FIXME: Distinguish between the different closing symbols. ")" does 
      //        terminate a context started with"[". Boerk
      if (t.GetValue() == ")" || t.GetValue() == "]" || t.GetValue() == "}") {
        break;
      }
    } else if (t.GetType() == PythonToken::OPERATOR && t.GetValue() != ".") {
       parent_context->AddContext(current_context);
       current_context = new PythonContext;
       current_context->SetRange(t.GetRange());
       parent_context->AddContext(current_context);       
       current_context = new PythonContext;       
       current_context->SetRange(Range(t.GetRange().End(),0));
    } else {
      current_context->AddToken(t);
    } 
  }
  parent_context->AddContext(current_context);
}

bool PythonContext::GetContextAtLocation(size_t location, 
                                         PythonContext*& context)
{
  if ( location < range_.location || location > range_.End())
    return false;
  for (size_t i = 0; i < sub_contexts_.size(); ++i) {
    if (sub_contexts_[i]->GetContextAtLocation(location,context) == true)
      return true;
  }
  context = this;
  return true;
}

}}
