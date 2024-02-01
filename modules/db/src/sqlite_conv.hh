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
#ifndef OST_IO_SQLITE_CONV_HH
#define OST_IO_SQLITE_CONV_HH

#include <sqlite3.h>

#include <ost/db/module_config.hh>

namespace ost { namespace db {
  
inline void SqlConvert(sqlite3_stmt* stmt, int& value, int col=0)
{
  value=sqlite3_column_int(stmt, col);
}

inline void SqlConvert(sqlite3_stmt* stmt, String& value, int col=0)
{
  value=String(reinterpret_cast<const char*>(sqlite3_column_text(stmt, col)));
}

inline void SqlConvert(sqlite3_stmt* stmt, sqlite3_int64& value, int col=0)
{
  value=sqlite3_column_int64(stmt, col);
}

template <typename T>
void SqlConvert(sqlite3_stmt* stmt, T& value)
{
  value.FromSql(stmt);
}

}}

#endif
