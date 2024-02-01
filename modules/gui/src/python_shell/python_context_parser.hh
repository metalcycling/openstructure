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
#ifndef PYTHON_CONTEXT_PARSER_HH
#define PYTHON_CONTEXT_PARSER_HH

/*
  Author: Marco Biasini
*/

#include "python_tokenizer.hh"
#include <vector>
#include <iostream>


namespace ost { namespace gui {
    
class DLLEXPORT_OST_GUI PythonContext
{
public:
  PythonContext() 
    : parent_(NULL)
  { }
  
  ~PythonContext()
  { this->Clear(); }
  
  bool IsRoot() const
  {
    return parent_ == NULL;
  }
  
  const Range& GetRange() const
  { return range_; }
  
  const std::vector<PythonContext*>& GetSubContexts() const
  { return sub_contexts_; }
  
  void SetRange( const Range& range)
  { range_ = range; }
  
  const PythonContext* GetParent() const
  { return parent_; }
  
  void AddToken(const PythonToken& token)
  { 
    range_.length = token.GetRange().End()-range_.location; 
    token_list_.push_back(token);
  }
  
  void Print()
  {
    for (size_t i = 0; i < token_list_.size(); ++i) {
      std::cout << token_list_[i].GetValue().toStdString() << std::endl;
    }
  }
  
  bool FindLastToken(PythonToken::Type type, 
                     const QString& value, 
                     PythonToken& token);
                     
  bool FindFirstToken(PythonToken::Type type, 
                      const QString& value, 
                      PythonToken& token);
                      
  void AddContext(PythonContext* context)
  { 
    sub_contexts_.push_back(context); 
    sub_contexts_.back()->parent_ = this;
  }  
  
  bool GetContextAtLocation(size_t location, PythonContext*& context);  
  
  void Clear();  
private:

  Range                           range_;
  std::vector<PythonContext*>     sub_contexts_;
  PythonContext*                  parent_;
  std::vector<PythonToken>        token_list_;
};

class DLLEXPORT_OST_GUI PythonContextParser
{
public:
  PythonContextParser();
  
  ~PythonContextParser()
  {
    delete root_;
  }
  
  void Print()
  {
    this->PrintRecursively(root_,0);
  }
  
  void Parse(const QString& command);
  
  const QString& GetCommand() const
  { return command_; }  
  
  bool GetContextAtLocation(size_t location, 
                            PythonContext*& context) 
  {
    return root_->GetContextAtLocation(location,context);
  }
  
protected:
  void Parse(PythonContext*& context);
  

  
private:
  void PrintRecursively(const PythonContext* context, int level)
  {
    typedef std::vector<PythonContext*>::const_iterator const_it;
    if (context->GetRange().length > 0) {
      for (int j = 0; j < level; ++j)
        std::cout<< "-" ;      
      try {
        std::cout << command_.mid(context->GetRange().location, 
                                  context->GetRange().length).toStdString() 
                  << std::endl;;      
      } catch(...) { }
    }

    const_it i = context->GetSubContexts().begin();
    for (;i!=context->GetSubContexts().end(); ++i) {
      this->PrintRecursively(*i,level+1);
    }
  }
  QString                       command_;
  PythonTokenizer               tokenizer_;
  PythonContext*                root_;
};

}}
#endif
