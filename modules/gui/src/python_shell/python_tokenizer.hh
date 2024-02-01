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
#ifndef PYTHON_TOKENIZER_HH
#define PYTHON_TOKENIZER_HH

/*
  Author: Marco Biasini
*/


#include <ost/gui/module_config.hh>
#include <QString>

namespace ost { namespace gui {

struct DLLEXPORT_OST_GUI Range
{
  size_t location, length;
  
  Range(size_t loc, size_t len):
    location(loc),
    length(len)
  { 
  }
  
  Range():
    location(-1),
    length(0)
  {      
  }

  size_t End() const
  {
    return location+length;
  }
};

class DLLEXPORT_OST_GUI PythonToken
{
public:
  enum Type {
    KEYWORD,
    STRING_LITERAL,
    STRING_DELIM,
    IDENTIFIER,
    NUMBER,
    OPERATOR,
    GROUPING,
    COMMENT,
    END
  };
  PythonToken();
  PythonToken(Type type,const Range& range, const QString& value);
  Type GetType() const;
  QString GetTypeAsString() const;
  const QString& GetValue() const;
  const Range& GetRange() const;
  
private:
  Type       type_;
  Range      range_;
  QString    value_;
};

class DLLEXPORT_OST_GUI PythonTokenizer
{
public:
  PythonTokenizer();
  PythonTokenizer(const QString& command, int string_state=0);
  void SetCommand(const QString& command);
  PythonToken NextToken();
  const PythonToken& CurrentToken() const;
  void SetInString(bool in_string);
  bool InString() const { return in_string_; }
  const QString& GetDelim() const { return string_delim_; }
private:
  void EatWhities();
  QString GetSubString(const Range& range);
  PythonToken GetIdentifierSeq();
  PythonToken GetNumberToken();
  PythonToken GetKeywordToken();
  PythonToken GetOperatorToken();
  PythonToken GetGroupingToken();
  PythonToken GetStringLiteral();
  PythonToken GetStringDelim();
  QString        command_;
  int            current_pos_;
  bool           in_string_;
  QString        string_delim_;
  PythonToken    current_token_;  
};

  
}}

#endif
