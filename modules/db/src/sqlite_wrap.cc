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
#include "sqlite_wrap.hh"

namespace ost { namespace db {
  
  
DatabasePtr Database::Open(const String& database) 
{
  DatabasePtr db(new Database);
  int retval=sqlite3_open(database.c_str(), &db->conn_);
  if (SQLITE_OK==retval) {
    return db;
  }
  throw DatabaseError(String(sqlite3_errmsg(db->conn_)));
}
  
PreparedStatement::PreparedStatement():
 statement_(NULL)
{    
}

PreparedStatementPtr Database::Prepare(const String& statement) 
{
  PreparedStatementPtr stmt(new PreparedStatement());
  int retval=sqlite3_prepare_v2(conn_, statement.c_str(), 
                                static_cast<int>(statement.length()),
                                &stmt->statement_, NULL);
  if (SQLITE_OK!=retval) {
    throw DatabaseError(String(sqlite3_errmsg(conn_)));
  }
  return stmt;
}

Database::~Database() 
{
  if (conn_) {
    sqlite3_close(conn_);
  }
}

PreparedStatement::~PreparedStatement()
{
  if (statement_)
    sqlite3_finalize(statement_);
}

DatabasePtr Database::Copy(const String& filename)
{
  DatabasePtr clone=Database::Open(filename);
  sqlite3_backup* backup;
  
  backup=sqlite3_backup_init(clone->conn_, "main", conn_, "main");
  if(backup){
    sqlite3_backup_step(backup, -1);
    sqlite3_backup_finish(backup);
  }
  int rc=sqlite3_errcode(clone->conn_);
  if (rc!=SQLITE_OK) {
    throw DatabaseError(sqlite3_errmsg(clone->conn_));    
  }
  return clone;
}

sqlite3_int64 Database::LastRowID()
{
  return sqlite3_last_insert_rowid(conn_);
}

}}
