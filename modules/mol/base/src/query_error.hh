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
#ifndef OST_QUERY_ERROR_HH
#define OST_QUERY_ERROR_HH

#include <ost/message.hh>
#include <ost/range.hh>
#include <ost/mol/module_config.hh>

namespace ost { namespace mol {
  
struct DLLEXPORT_OST_MOL QueryErrorDesc {
  QueryErrorDesc(const QueryErrorDesc& _rhs);

  QueryErrorDesc();

  QueryErrorDesc& operator=(const QueryErrorDesc& _rhs);

  String      msg;
  Range       range;
};

class DLLEXPORT_OST_MOL QueryError : public ost::Error {
public:
  QueryError(const String& query, const QueryErrorDesc& desc) throw();    

  ~QueryError() throw();

  //! Get error description.
  const String& GetMessage()const throw();
  
  // Get range of where the error occurred
  const Range& GetRange() const;

  //! Get nicely formatted error message
  String GetFormattedMessage() const throw();
  //! Get original query String
  const String& GetQueryString() const throw();  
private:
  String             query_;
  QueryErrorDesc     desc_;
};

}} //ns
#endif
