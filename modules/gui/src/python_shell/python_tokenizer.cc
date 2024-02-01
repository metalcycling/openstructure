//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
  Authors: Marco Biasini
*/

#include "python_tokenizer.hh"

namespace ost { namespace gui {

PythonToken::PythonToken():
  type_(END),
  range_(String::npos,0),
  value_("")
{
}

PythonToken::PythonToken(Type type,const Range& range, const QString& value):
  type_(type),
  range_(range),
  value_(value)
{
}

PythonToken::Type PythonToken::GetType() const
{
  return type_;
}

QString PythonToken::GetTypeAsString() const
{
  switch(type_){
  case STRING_LITERAL:
    return "String literal";
  case STRING_DELIM:
    return "String delimiter";
  case NUMBER:
    return "number";
  case OPERATOR:
    return "operator";
  case IDENTIFIER:
    return "primary expression";
  case KEYWORD:
    return "keyword";
  case END:
    return "EOF";
  case GROUPING:
    return "grouping character";
  case COMMENT:
    return "comment";
  }
  return "unknown";
}
const QString& PythonToken::GetValue() const
{
  return value_;
}

const Range& PythonToken::GetRange() const
{
  return range_;
}


namespace {
PythonToken NullToken()
{
  static PythonToken NULL_TOKEN = PythonToken(PythonToken::END,Range(),"");
  return NULL_TOKEN;
}
}//anon ns

PythonTokenizer::PythonTokenizer():
  command_(""),
  current_pos_(0)
{
}  

PythonTokenizer::PythonTokenizer(const QString& command, int string_state):
  command_(command),
  current_pos_(0),
  in_string_(string_state >0)
{
  if (string_state==2) {
    string_delim_="'''";
  } else if (string_state==1) {
    string_delim_="\"\"\"";
  }
}


const PythonToken& PythonTokenizer::CurrentToken() const
{
  return current_token_;
}  

void PythonTokenizer::SetInString(bool in_string)
{
  in_string_=in_string;
}

QString PythonTokenizer::GetSubString(const Range& range)
{
  return command_.mid(range.location,range.length);
}

void PythonTokenizer::SetCommand(const QString& command)
{
  command_ = command;
  current_pos_ = 0;
  current_token_ = NullToken();
}

void PythonTokenizer::EatWhities()
{
  while(current_pos_ < command_.size() && command_[current_pos_].isSpace()){
    current_pos_++;
  }
}
PythonToken PythonTokenizer::GetOperatorToken()
{
  static QString operators[] = {
    QString("+="), QString("-="), QString("*="), 
    QString("**"), QString("&="), QString("|="),
    QString("/="), QString("+"), QString("-"), 
    QString("/"), QString("*"), QString("."),
    QString("%"), QString("="), QString("!="), 
    QString("!"), QString("<="), QString(">="), 
    QString("<"), QString(":"), QString("?"), 
    QString("|"), QString("&"), QString("^"),
    QString("")
  };
  size_t index = -1;
  while (operators[++index] != QString("")) {
    if (command_.indexOf( operators[index], current_pos_) == current_pos_){
      QString op = operators[index];
      Range range( current_pos_, op.size() );
      current_pos_+=op.size();
     return PythonToken(PythonToken::OPERATOR,range,op);
    }
  }
  return NullToken();  
}

PythonToken PythonTokenizer::GetKeywordToken()
{
  static QString keywords[] = { 
    QString("for"), QString("while"), QString("def"), QString("pass"), 
    QString("class"), QString("if"), QString("else"), QString("print"), 
    QString("self"), QString("float"), QString("import"), QString("from"),    
    QString("and"), QString("or"), QString("not"), QString("str"), 
    QString("False"), QString("True"), QString("None"), 
    QString("in"),  QString("elif"), QString("") };
    
  size_t index = -1;
  while (keywords[++index]!=QString("")) {
    if (command_.indexOf( keywords[index], current_pos_) == current_pos_) {
      bool keyword_accepted = false;
      QString kw = keywords[index];
      if (current_pos_+kw.size() < command_.size()){
        QChar char_after_kw = command_[current_pos_+kw.size()];
        if(!char_after_kw.isLetterOrNumber() && char_after_kw!='_'){
          keyword_accepted=true;
        }
      }
      else{
        keyword_accepted = true;
      }
      if (keyword_accepted){
        Range range( current_pos_, kw.size() );
        current_pos_+=kw.size();
        return PythonToken(PythonToken::KEYWORD,range,kw);
      }
    }
  }
  return NullToken();
}

PythonToken PythonTokenizer::GetIdentifierSeq()
{
  Range range(current_pos_,0);
  while (current_pos_ < command_.size() && (command_[current_pos_].isLetterOrNumber() || 
    command_[current_pos_] == '_')) {
    current_pos_++;
  }
  range.length = current_pos_-range.location;
  return PythonToken(PythonToken::IDENTIFIER,range,GetSubString(range));
}

PythonToken PythonTokenizer::GetStringDelim()
{
  Range range(current_pos_,0);
  // Find out which type of delimiters are used for the String. i.e. <">, 
  // <""">, <'''> or <'>
  QString delimiter = "\"";
  if (command_[current_pos_] == '\'')
    delimiter = "'";
  if (command_[current_pos_] == '"') {
    if (current_pos_+2 < command_.size()) {
      QString sub = this->GetSubString(Range(current_pos_, 3));
      if (sub=="\"\"\"") {
        delimiter="\"\"\"";
      }
    }
  }
  range.length=delimiter.length();
  current_pos_=range.End();
  string_delim_=delimiter;
  return PythonToken(PythonToken::STRING_DELIM,
                     range,string_delim_);
}

PythonToken PythonTokenizer::GetStringLiteral() {
  Range range(current_pos_,0);  
  int str_end =command_.indexOf(string_delim_, current_pos_);
  if (str_end!=-1) {
    range.length=str_end-range.location;
  } else {
    range.length=command_.size()-range.location;
  }
  current_pos_=range.End();  
  return PythonToken(PythonToken::STRING_LITERAL,
                     range,GetSubString(range));
}

PythonToken PythonTokenizer::GetGroupingToken()
{
  QString groupCharacters[] = { QString("("), 
                               QString(")"), 
                               QString("{"), 
                               QString("}"), 
                               QString("["), 
                               QString("]"), 
                               QString(","), 
                               QString(";"), QString("") };
  int index = -1;
  Range range(current_pos_,0);
  while( groupCharacters[++index] != QString("") ) {
    if (command_.indexOf(groupCharacters[index],current_pos_) == current_pos_) {
      QString groupCharacter = groupCharacters[index];
      current_pos_+=groupCharacter.size();
      range.length = groupCharacter.size();
      return PythonToken(PythonToken::GROUPING,range,groupCharacter);
    }
  }
  return NullToken();
}

PythonToken PythonTokenizer::NextToken()
{
  if (command_.length()>current_pos_ && in_string_) {
    QChar current_char = command_[current_pos_];
    if (current_char == '"' || current_char == '\'') {
      current_token_ = GetStringDelim();
      in_string_=false;
      return current_token_;
    }    
    current_token_=GetStringLiteral();
    return current_token_;
  } else {
    this->EatWhities();
    if (current_pos_ < command_.size()) {
      QChar current_char = command_[current_pos_];
      if (current_char=='#') {
        current_token_=PythonToken(PythonToken::COMMENT, 
                                   Range(current_pos_, 
                                         command_.size()-current_pos_),
                                   command_.mid(current_pos_, 
                                                command_.size()-current_pos_));
        current_pos_=command_.size();
        return current_token_;
      }
      if (current_char.isDigit()) {
        current_token_ = GetNumberToken();
        return current_token_;
      }    
      if (current_char.isLetterOrNumber() || current_char == '_') {
        PythonToken keywordToken = this->GetKeywordToken();
        if (keywordToken.GetType() == PythonToken::END){
          current_token_ = this->GetIdentifierSeq();        
          return current_token_;
        }
        current_token_ = keywordToken;
        return current_token_;
      }
      PythonToken opToken = this->GetOperatorToken();
      if (opToken.GetType() != PythonToken::END) {
        current_token_ = opToken;
        return current_token_;
      }
      if (current_char == '"' || current_char == '\'') {
        current_token_ = this->GetStringDelim();
        in_string_=true;
        return current_token_;
      }
      current_token_ = this->GetGroupingToken();
      return current_token_;
    } else{
      current_token_ = NullToken();
      return current_token_;
    }    
  }
  current_token_=NullToken();
  return current_token_;
}

PythonToken PythonTokenizer::GetNumberToken()
{
  Range range(current_pos_,0);
  bool dot_present = false;
  while(current_pos_ < command_.size()){
    if (command_[current_pos_].isDigit()){
      current_pos_++;
    } else if (command_[current_pos_] == '.'){
      if (dot_present) {
        return NullToken();        
      }
      dot_present = true;
      current_pos_++;
    } else {
      break;
    }
  }
  range.length=current_pos_-range.location;
  return PythonToken(PythonToken::NUMBER,range,this->GetSubString(range));
}

}}

