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
#ifndef OST_IO_SQLITE_WRAP_HH
#define OST_IO_SQLITE_WRAP_HH
/*
 Author: Marco Biasini
 
 lightweight sqlite wrapper
*/
#include <ost/message.hh>
#include <boost/shared_ptr.hpp>

#include <ost/base.hh>
#include <ost/message.hh>

#include <sqlite3.h>
#include "sqlite_conv.hh"

namespace ost { namespace db {

/// \brief Database connection
class Database;
typedef boost::shared_ptr<Database> DatabasePtr;

class PreparedStatement;
typedef boost::shared_ptr<PreparedStatement> PreparedStatementPtr;

class DLLEXPORT_OST_DB PreparedStatement {
  friend class Database;
public:
  PreparedStatement();
  ~PreparedStatement();
    
  sqlite3_stmt* Handle() 
  {
    return statement_;
  }
  /// \brief submit query after binding all a parameters.
  bool Submit()
  {
    return sqlite3_step(statement_)==SQLITE_DONE;
  }
  /// \brief bind integer value to placeholder
  void BindInt(int col, int value) 
  {
    sqlite3_bind_int(statement_, col, value);
  }
  /// \brief bind 64 bit integer value to placeholder
  void BindInt64(int col, sqlite3_int64 value) 
  {
    sqlite3_bind_int64(statement_, col, value);
  }  
  /// \brief bind String value to placeholder  
  void BindText(int col, const String& text) 
  {
    sqlite3_bind_text(statement_, col, text.c_str(),
                      static_cast<int>(text.size()), NULL);
  }  
  /// \brief bind binary blob  to placeholder  
  void BindBlob(int col, const char* start, size_t size)
  {
    sqlite3_bind_blob(statement_, col, start, size, NULL);
  }
private:
  sqlite3_stmt* statement_;
};

class DLLEXPORT DatabaseError : public Error {
public:
  DatabaseError(const String& msg):
    Error(msg)
  {
    
  }
};

/// \brief  Row set iterator
/// 
/// Row set iterators let you iterate over a set of rows as returned by a
/// select statements. The parametrisable type R is the type that should be used
/// for each row. For the standard types such as int and strings, conversion
/// routines are available. For your custom types, you have to write an
/// overloaded version of SqlConvert()
/// 
/// \code
/// PreparedStatementPtr ps=db->Prepare("select * from animals where num_legs=?");
/// ps->BindInt(num_legs)
/// 
/// for (RowSet<Animal> ani_it(ps); !ani_it.AtEnd(); ++ani_it) {
///   Animal& a=*ani_it;
///   std::cout << a.GetName() << std::endl;
/// }
/// \endcode
template <typename R>
class DLLEXPORT_OST_DB RowSet {
public:
  RowSet(PreparedStatementPtr stmt): 
    stmt_(stmt), value_()
  {
    this->Step();
  }
  RowSet():
    stmt_(), value_(), cur_val_(SQLITE_DONE) 
  { }
  R& operator*() {
    return value_;
  }
  RowSet<R>& operator++() 
  {
    this->Step();
    return *this;
  }
  
  bool AtEnd() 
  {
    return cur_val_!=SQLITE_ROW;
  }
private:
  void Step() 
  {
    cur_val_=sqlite3_step(stmt_->Handle());
    if (cur_val_==SQLITE_ROW) {
      SqlConvert(stmt_->Handle(), value_);
    } 
  }
  PreparedStatementPtr stmt_;
  R                    value_;
  int                  cur_val_;
};

/// \brief SQLite database handle
class DLLEXPORT_OST_DB Database {
public:
  /// \brief open database
  /// 
  /// \return a valid database handle
  /// \throw DatabaseError if the database could not be opened.
  static DatabasePtr Open(const String& database);
  
  ~Database();
  /// \brief prepare statement for execution
  /// \throw DatabaseError if the statement could not be prepared.
  /// \note execution of multiple statemens is not supported.
  PreparedStatementPtr Prepare(const String& query);
  
  /// \brief Get SQLite handle
  sqlite3* Handle() 
  {
    return conn_;
  }
  
  sqlite3_int64 LastRowID();
  /// \brief clone database
  /// 
  /// useful for saving an in-memory database to disk
  DatabasePtr Copy(const String& filename);
private:
  sqlite3*     conn_;
};
 
}}

#endif
