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
#include "query_error.hh"

namespace ost { namespace mol {
  
String QueryError::GetFormattedMessage() const throw() {
  String formatted =query_+"\n";
  if (desc_.range.Loc != String::npos) {
    for (size_t i = 0; i < desc_.range.Loc;++i) {
      formatted+=" ";
    }
    if ( desc_.range.Length == 0) {
      formatted+="^";
    }
    for ( size_t i = 0; i < desc_.range.Length;++i) {
      formatted+="~";
    }
    formatted+="\n";
  }
  return "ERROR: "+desc_.msg+"\n"+formatted;
} 

QueryError::QueryError(const String& query, 
                       const QueryErrorDesc& desc) throw()
  : ost::Error(""),query_(query),desc_(desc) {
}
  
QueryError::~QueryError() throw() { 
}

const String& QueryError::GetMessage() const throw() { 
  return desc_.msg; 
}

const Range& QueryError::GetRange() const { 
  return desc_.range; 
}  


const String& QueryError::GetQueryString() const throw() {
  return query_;
} 

QueryErrorDesc::QueryErrorDesc() {
  
}

QueryErrorDesc& QueryErrorDesc::operator=(const QueryErrorDesc& _rhs) {
  this->msg = _rhs.msg;
  this->range = _rhs.range;
  return *this;
}

QueryErrorDesc::QueryErrorDesc(const QueryErrorDesc& _rhs) 
  : msg(_rhs.msg),range(_rhs.range) { 
}

}} //ns
