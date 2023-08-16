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
  Author: Marco Biasini
 */

#include <iostream>
#include <sstream>
#include <boost/format.hpp>
#include <ost/log.hh>
#include "compound_lib.hh"
#include <ost/version.hh>
#include <ost/string_ref.hh>
#include <sqlite3.h>

using boost::format;

namespace ost { namespace conop {

namespace {

/*
This is the oldest version (GetCreationDate) of compound libraries we support.
*/
const String COMPAT_VERSION = "1.5.0";

/*
COMMENT ON CREATE_CMD

CREATE_CMD specifies so called affinities. e.g. VARCHAR(64) where common sense
interprets the number 64 as max length of entries in that column.
However, sqlite 3 totally ignores this an interprets it as TEXT without any
limits. Long story short, don't worry about formulas longer than 64 characters
or longer names etc.
*/

const char* CREATE_CMD[]={
"CREATE TABLE IF NOT EXISTS chemlib_info (                                      "
"  creation_date     TIMESTAMP,                                                 "
"  ost_version_used  VARCHAR(64) NOT NULL);",
"CREATE TABLE IF NOT EXISTS chem_compounds (                                    "
"  id                INTEGER PRIMARY KEY AUTOINCREMENT,                         "
"  tlc               VARCHAR(5) NOT NULL,                                       "
"  olc               VARCHAR(1) NOT NULL,                                       "
"  dialect           VARCHAR(1) NOT NULL,                                       "
"  chem_class        VARCHAR(1),                                                "
"  chem_type         VARCHAR(1),                                                "
"  formula           VARCHAR(64) NOT NULL,                                      "
"  pdb_initial       TIMESTAMP,                                                 "
"  pdb_modified      TIMESTAMP,                                                 "
"  name              VARCHAR(256),                                              "
"  inchi_code        TEXT,                                                      "
"  inchi_key         TEXT,                                                      "
"  smiles            TEXT,                                                      "
"  obsolete          BOOL,                                                      "
"  replaced_by       VARCHAR(5)                                                 "
");",
" CREATE UNIQUE INDEX IF NOT EXISTS compound_tlc_index ON chem_compounds        "
"                                  (tlc, dialect)",
" CREATE INDEX IF NOT EXISTS compound_smiles_index ON chem_compounds            "
"                                  (smiles, dialect)",
" CREATE INDEX IF NOT EXISTS compound_inchi_code_index ON chem_compounds        "
"                                  (inchi_code, dialect)",
" CREATE INDEX IF NOT EXISTS compound_inchi_key_index ON chem_compounds         "
"                                  (inchi_key, dialect)",
"CREATE TABLE IF NOT EXISTS atoms (                                             "
" id                 INTEGER PRIMARY KEY AUTOINCREMENT,                         "
" compound_id        INTEGER REFERENCES chem_compounds (id) ON DELETE CASCADE,  "
" name               VARCHAR(4) NOT NULL,                                       "
" alt_name           VARCHAR(4) NOT NULL,                                       "
" element            VARCHAR(2) NOT NULL,                                       "
" is_aromatic        VARCHAR(1) NOT NULL,                                       "
" is_leaving         VARCHAR(1) NOT NULL,                                       "
" ordinal            INT,                                                       "
" charge             INT                                                        "
");",
" CREATE INDEX IF NOT EXISTS atom_name_index ON atoms                           "
"                                  (compound_id, name, alt_name)",
" CREATE TABLE IF NOT EXISTS bonds  (                                           "
"   id               INTEGER PRIMARY KEY AUTOINCREMENT,                         "
"   compound_id      INTEGER REFERENCES chem_compounds (id) ON DELETE CASCADE,  "
"   atom_one         INTEGER REFERENCES atoms (id) ON DELETE CASCADE,           "
"   atom_two         INTEGER REFERENCES atoms (id) ON DELETE CASCADE,           "
"   bond_order       INT                                                        "
" );",
" CREATE INDEX IF NOT EXISTS bond_index ON bonds (compound_id)",
" CREATE TRIGGER delete_related_objects                                         "
" BEFORE DELETE ON chem_compounds                                               "
" FOR EACH ROW BEGIN                                                            "
"     DELETE FROM bonds WHERE compound_id=OLD.id;                               "
"     DELETE FROM atoms WHERE compound_id=OLD.id;                               "
" END",
" DROP TRIGGER IF EXISTS delete_related_objects", NULL};    


const char* INSERT_COMPOUND_STATEMENT="INSERT INTO chem_compounds               "
"        (tlc, olc, dialect, chem_class, chem_type, formula, pdb_initial,       "
"         pdb_modified, name, inchi_code, inchi_key, smiles, obsolete,          "
"         replaced_by) "
" VALUES (?, ?, ?, ?, ?, ?, DATE(?), DATE(?), ?, ?, ?, ?, ?, ?)";

const char* INSERT_ATOM_STATEMENT="INSERT INTO atoms                            "
"        (compound_id, name, alt_name, element, is_aromatic,                    "
"         is_leaving, ordinal, charge)                                          "
" VALUES (?, ?, ?, ?, ?, ?, ?, ?)";

const char* INSERT_BOND_STATEMENT="insert into bonds                            "
"        (compound_id, atom_one, atom_two, bond_order)                          "
" VALUES (?, ?, ?, ?)";

const char* INSERT_CHEMLIB_INFO_STATEMENT="insert into chemlib_info             "
"        (creation_date, ost_version_used)                                      "
" VALUES (DATE(?), ?)";

}

struct CompoundLib::Database {

  Database(): ptr(NULL) { }

  Database(sqlite3* p): ptr(p) { }

  ~Database() {
    if (ptr) {
      int retval = sqlite3_close(ptr);
      if (retval != SQLITE_OK) {
        LOG_ERROR("Problem while closing SQLite db for CompoundLib: "
                  << sqlite3_errmsg(ptr));
      }
    }
  }

  sqlite3* ptr;
};


void CompoundLib::SetChemLibInfo(void){
  sqlite3_stmt* stmt=NULL;
  //~ if (!db_->ptr) {
    //~ LOG_ERROR(sqlite3_errmsg("Connection to DB not made"));
  //~ }
  int retval=sqlite3_prepare_v2(db_->ptr, INSERT_CHEMLIB_INFO_STATEMENT, 
                                strlen(INSERT_CHEMLIB_INFO_STATEMENT), &stmt, NULL);
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  Date date = Date(1900+timeinfo->tm_year, 1+timeinfo->tm_mon, timeinfo->tm_mday);
  String date_str=date.ToString();
  if (SQLITE_OK==retval) {

    sqlite3_bind_text(stmt, 1, date_str.c_str(), 
                      strlen(date_str.c_str()), NULL);
    const char* ost_version = OST_VERSION_STRING;
    sqlite3_bind_text(stmt, 2, ost_version, 
                      strlen(ost_version), NULL);
  } else {
    std::cout << "failed" <<std::endl;
  }
  retval=sqlite3_step(stmt);
  if (SQLITE_DONE!=retval) {
    if (sqlite3_errcode(db_->ptr)==SQLITE_CONSTRAINT) {
      LOG_ERROR("chemlib info already exists");
    } else {
      LOG_ERROR(sqlite3_errmsg(db_->ptr));
    }
  }
  sqlite3_finalize(stmt);  
}

Date CompoundLib::GetCreationDate(void){
  
  String query="SELECT creation_date FROM chemlib_info";
  sqlite3_stmt* stmt;
  int retval=sqlite3_prepare_v2(db_->ptr, query.c_str(), 
                                static_cast<int>(query.length()),
                                &stmt, NULL);
  if (SQLITE_OK==retval) {
    int ret=sqlite3_step(stmt);
    if (SQLITE_DONE==ret) {
      sqlite3_finalize(stmt);        
      return Date();
    }
    if (SQLITE_ROW==ret) {
      const char* strr = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
      Date date = Date::FromString(StringRef(strr, strlen(strr))); 
      sqlite3_finalize(stmt);
      return date;   
    }
    assert(SQLITE_DONE==sqlite3_step(stmt));
  } else {
    sqlite3_finalize(stmt);      
    return Date();
  }
  sqlite3_finalize(stmt);    
  return Date();
}


String CompoundLib::GetOSTVersionUsed() {
  String query="SELECT ost_version_used FROM chemlib_info";
  sqlite3_stmt* stmt;
  String version;

  int retval=sqlite3_prepare_v2(db_->ptr, query.c_str(), 
                                static_cast<int>(query.length()),
                                &stmt, NULL);
  if (SQLITE_OK==retval) {
    int ret=sqlite3_step(stmt);
    if (SQLITE_DONE==ret) {
      sqlite3_finalize(stmt);        
      return String();
    }
    if (SQLITE_ROW==ret) {
      version = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0))); 
      sqlite3_finalize(stmt);
      return version;   
    }
    assert(SQLITE_DONE==sqlite3_step(stmt));
  } else {
    LOG_WARNING("your compound library might be outdated.");
    sqlite3_finalize(stmt);      
    return String();
  }
  sqlite3_finalize(stmt);    
  return String();
}


void CompoundLib::AddCompound(const CompoundPtr& compound)
{
  sqlite3_stmt* stmt=NULL;  
  int retval=sqlite3_prepare_v2(db_->ptr, INSERT_COMPOUND_STATEMENT, 
                                strlen(INSERT_COMPOUND_STATEMENT), &stmt, NULL);
  String crea_date_str, modi_date_str;
  if (SQLITE_OK==retval) {
    sqlite3_bind_text(stmt, 1, compound->GetID().c_str(), 
                      compound->GetID().length(), NULL);
    char olc=compound->GetOneLetterCode();
    sqlite3_bind_text(stmt, 2, &olc, 1, NULL);
    char chem_class=compound->GetChemClass();
    char chem_type=compound->GetChemType();
    char dialect=compound->GetDialect();
    sqlite3_bind_text(stmt, 3, &dialect, 1, NULL);
    sqlite3_bind_text(stmt, 4, &chem_class, 1, NULL);
    sqlite3_bind_text(stmt, 5, &chem_type, 1, NULL);
    sqlite3_bind_text(stmt, 6, compound->GetFormula().c_str(),
                      compound->GetFormula().length(), NULL);
    Date crea_date=compound->GetCreationDate();
    Date modi_date=compound->GetModificationDate();
    crea_date_str=crea_date.ToString();
    modi_date_str=modi_date.ToString();
    sqlite3_bind_text(stmt, 7, crea_date_str.c_str(), crea_date_str.length(),
                      NULL);
    sqlite3_bind_text(stmt, 8, modi_date_str.c_str(), modi_date_str.length(),
                      NULL);
    sqlite3_bind_text(stmt, 9, compound->GetName().c_str(),
                      compound->GetName().length(), NULL);
    sqlite3_bind_text(stmt, 10, compound->GetInchi().c_str(),
                      compound->GetInchi().length(), NULL);
    sqlite3_bind_text(stmt, 11, compound->GetInchiKey().c_str(),
                      compound->GetInchiKey().length(), NULL);
    sqlite3_bind_text(stmt, 12, compound->GetSMILES().c_str(),
                      compound->GetSMILES().length(), NULL);
    sqlite3_bind_int(stmt, 13, compound->GetObsolete());
    if (compound->GetReplacedBy() != "") {
      sqlite3_bind_text(stmt, 14, compound->GetReplacedBy().c_str(),
                        compound->GetReplacedBy().length(), NULL);
    }
  } else {
    LOG_ERROR(sqlite3_errmsg(db_->ptr));
    sqlite3_finalize(stmt);
    return;
  }
  retval=sqlite3_step(stmt);
  if (SQLITE_DONE!=retval) {
    if (sqlite3_errcode(db_->ptr)==SQLITE_CONSTRAINT) {
      LOG_ERROR("Compound '" << compound->GetID() << "' already exists for the "
                << compound->GetDialectAsString() << " dialect.");
    } else {
      LOG_ERROR(sqlite3_errmsg(db_->ptr));
    }
  }
  sqlite3_finalize(stmt);  
  sqlite3_int64 compound_id=sqlite3_last_insert_rowid(db_->ptr);
  // insert atoms
  const AtomSpecList& al=compound->GetAtomSpecs();
  std::vector<sqlite3_int64> atom_ids(al.size(), 0);
  for (AtomSpecList::const_iterator i=al.begin(), e=al.end(); i!=e; ++i) {
    const AtomSpec& a=*i;
    retval=sqlite3_prepare_v2(db_->ptr, INSERT_ATOM_STATEMENT, 
                              strlen(INSERT_ATOM_STATEMENT), &stmt, NULL);
    if (SQLITE_OK==retval) {
      sqlite3_bind_int64(stmt, 1, compound_id);
      sqlite3_bind_text(stmt, 2, a.name.c_str(), a.name.length(), NULL);
      sqlite3_bind_text(stmt, 3, a.alt_name.c_str(), 
                        a.alt_name.length(), NULL);
      sqlite3_bind_text(stmt, 4, a.element.c_str(), a.element.length(), NULL);                        
      sqlite3_bind_int(stmt, 5, a.is_aromatic);
      sqlite3_bind_int(stmt, 6, a.is_leaving);
      sqlite3_bind_int(stmt, 7, a.ordinal);
      sqlite3_bind_int(stmt, 8, a.charge);
      retval=sqlite3_step(stmt);
      if (retval != SQLITE_DONE) {
        LOG_ERROR(sqlite3_errmsg(db_->ptr));
      }
      atom_ids[a.ordinal]=sqlite3_last_insert_rowid(db_->ptr);
    } else {
      LOG_ERROR(sqlite3_errmsg(db_->ptr));
    }
    sqlite3_finalize(stmt);
  }
  const BondSpecList& bl=compound->GetBondSpecs();
  for (BondSpecList::const_iterator i=bl.begin(), e=bl.end(); i!=e; ++i) {
    const BondSpec& b=*i;
    retval=sqlite3_prepare_v2(db_->ptr, INSERT_BOND_STATEMENT, 
                              strlen(INSERT_BOND_STATEMENT), &stmt, NULL);    
    if (SQLITE_OK==retval) {
      sqlite3_bind_int64(stmt, 1, compound_id);
      sqlite3_bind_int64(stmt, 2, atom_ids[b.atom_one]);
      sqlite3_bind_int64(stmt, 3, atom_ids[b.atom_two]);      
      sqlite3_bind_int(stmt, 4, b.order);
      retval=sqlite3_step(stmt);
      if (retval != SQLITE_DONE) {
        LOG_ERROR(sqlite3_errmsg(db_->ptr));
      };
    } else {
      LOG_ERROR(sqlite3_errmsg(db_->ptr));
    }
    sqlite3_finalize(stmt);    
  }
}

CompoundLibPtr CompoundLib::Copy(const String& filename) const
{
  CompoundLibPtr clone=CompoundLibPtr(new CompoundLib);
  int retval=sqlite3_open(filename.c_str(), &clone->db_->ptr);
  if (SQLITE_OK==retval) {
    sqlite3_backup* backup;

    backup=sqlite3_backup_init(clone->db_->ptr, "main", db_->ptr, "main");
    if (backup){
      sqlite3_backup_step(backup, -1);
      sqlite3_backup_finish(backup);
    }
    int rc=sqlite3_errcode(clone->db_->ptr);
    if (rc!=SQLITE_OK) {
      LOG_ERROR(sqlite3_errmsg(clone->db_->ptr));
      return CompoundLibPtr();
    }
    return clone;      
  }
  LOG_ERROR(sqlite3_errmsg(clone->db_->ptr));
  return CompoundLibPtr();
}



CompoundLibPtr CompoundLib::Create(const String& database)
{
  CompoundLibPtr lib(new CompoundLib);
  int retval=sqlite3_open(database.c_str(), &lib->db_->ptr);
  if (SQLITE_OK==retval) {
    const char** cmd=CREATE_CMD;
    while (*cmd) {
      sqlite3_stmt* stmt;
      retval=sqlite3_prepare_v2(lib->db_->ptr, *cmd, strlen(*cmd), &stmt, NULL);
      if (SQLITE_OK==retval) {
        retval=sqlite3_step(stmt);
        sqlite3_finalize(stmt);        
        assert(SQLITE_DONE==retval);
      } else {
        LOG_ERROR(sqlite3_errmsg(lib->db_->ptr));
        sqlite3_finalize(stmt);
        return CompoundLibPtr();
      }      
      ++cmd;
    }
    return lib;
  }
  LOG_ERROR(sqlite3_errmsg(lib->db_->ptr));
  return CompoundLibPtr();  
}


CompoundLibPtr CompoundLib::Load(const String& database, bool readonly) 
{
  int flags=readonly ? SQLITE_OPEN_READONLY : SQLITE_OPEN_READWRITE;
  CompoundLibPtr lib(new CompoundLib);
  int retval=sqlite3_open_v2(database.c_str(), &lib->db_->ptr, flags, NULL);
  if (SQLITE_OK!=retval) {
    LOG_ERROR(sqlite3_errmsg(lib->db_->ptr));
    return CompoundLibPtr();
  }
  String aq;
  sqlite3_stmt* stmt;
  // check if SMILES are available
  aq="SELECT smiles FROM chem_compounds LIMIT 1";
  retval=sqlite3_prepare_v2(lib->db_->ptr, aq.c_str(),
                            static_cast<int>(aq.length()),
                            &stmt, NULL);
  lib->smiles_available_ = retval==SQLITE_OK;
  sqlite3_finalize(stmt);

  // check if obsolete info are available
  aq="SELECT obsolete, replaced_by FROM chem_compounds LIMIT 1";
  retval=sqlite3_prepare_v2(lib->db_->ptr, aq.c_str(),
                            static_cast<int>(aq.length()),
                            &stmt, NULL);
  lib->obsolete_available_ = retval==SQLITE_OK;
  sqlite3_finalize(stmt);

  // check if charges are available
  aq="SELECT charge FROM atoms LIMIT 1";
  retval=sqlite3_prepare_v2(lib->db_->ptr, aq.c_str(),
                            static_cast<int>(aq.length()),
                            &stmt, NULL);
  lib->charges_available_ = retval==SQLITE_OK;
  sqlite3_finalize(stmt);

  lib->creation_date_ = lib->GetCreationDate();
  lib->ost_version_used_ = lib->GetOSTVersionUsed();

  // Report compatibility issues
  if (lib->ost_version_used_.compare(COMPAT_VERSION) < 0) {
    std::stringstream ss;
    ss << "Compound lib was created with an unsupported version of OST: "
       << lib->ost_version_used_
       << ". Please update your compound library.";
    throw ost::Error(ss.str());
  }
  if (!lib->smiles_available_) {
    LOG_WARNING("SMILES not available in compound library v."
                << lib->ost_version_used_
                << ". Only empty strings will be returned.");
  }
  if (!lib->obsolete_available_) {
    LOG_WARNING("Obsolete information not available in compound library v."
                << lib->ost_version_used_
                << ". No compound will be marked as obsolete.");
  }
  if (!lib->charges_available_) {
    LOG_WARNING("Charges not available in compound library v."
                << lib->ost_version_used_
                << ". All charges will be 0.");
  }
  return lib;
}

void CompoundLib::LoadAtomsFromDB(CompoundPtr comp, int pk) const {
  String aq="SELECT name, alt_name, element, ordinal, is_leaving";
  if (charges_available_) {
    aq+=", charge";
  }
  aq = str(format(aq +
                  " FROM atoms WHERE compound_id=%d"
                  " ORDER BY ordinal ASC") % pk);
  sqlite3_stmt* stmt;
  int retval=sqlite3_prepare_v2(db_->ptr, aq.c_str(), 
                                static_cast<int>(aq.length()),
                                &stmt, NULL);
  if (SQLITE_OK==retval) {
      int ret=0;
      while (SQLITE_ROW==(ret=sqlite3_step(stmt))) {
        AtomSpec atom_sp;
        atom_sp.name=String(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0)));
        atom_sp.alt_name=String(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1)));        
        atom_sp.element=String(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2))); 
        atom_sp.ordinal=sqlite3_column_int(stmt, 3);  
        atom_sp.is_leaving=bool(sqlite3_column_int(stmt, 4)!=0);
        if (charges_available_) {
          atom_sp.charge=sqlite3_column_int(stmt, 5);
        }
        comp->AddAtom(atom_sp);
      }
  } else {
    LOG_ERROR(sqlite3_errmsg(db_->ptr));
  }
  sqlite3_finalize(stmt);
}
void CompoundLib::ClearCache()
{
  compound_cache_.clear();
}

void CompoundLib::LoadBondsFromDB(CompoundPtr comp, int pk) const {
  sqlite3_stmt* stmt;
  String aq=str(format("SELECT a1.ordinal, a2.ordinal, b.bond_order FROM bonds AS b "
                       "LEFT JOIN atoms AS a1 ON b.atom_one=a1.id "
                       "LEFT JOIN atoms as a2 ON b.atom_two=a2.id "
                       "WHERE b.compound_id=%d") % pk);
  int retval=sqlite3_prepare_v2(db_->ptr, aq.c_str(), 
                                static_cast<int>(aq.length()),
                                &stmt, NULL);
  if (SQLITE_OK==retval) {
      int ret=0;
      while (SQLITE_ROW==(ret=sqlite3_step(stmt))) {
        BondSpec bond_sp;
        bond_sp.atom_one=sqlite3_column_int(stmt, 0);
        bond_sp.atom_two=sqlite3_column_int(stmt, 1);
        bond_sp.order=sqlite3_column_int(stmt, 2);
        comp->AddBond(bond_sp);
      }
  } else {
    LOG_ERROR(sqlite3_errmsg(db_->ptr));
  } 
  sqlite3_finalize(stmt);
}

String CompoundLib::BuildFindCompoundQuery(const String& id,
                                           Compound::Dialect dialect,
                                           const String& by) const {

  // Build the query
  String query="SELECT id, tlc, olc, chem_class, dialect, formula, chem_type, name, inchi_code, inchi_key";
  if(smiles_available_) {
    query+=", smiles";
  }
  if(obsolete_available_) {
    query+=", obsolete, replaced_by";
  }
  query+=" FROM chem_compounds"
         " WHERE " + by + "=? AND dialect='"+String(1, char(dialect))+"'";

  return query;
}

CompoundPtr CompoundLib::LoadCompoundFromDB(sqlite3_stmt* stmt) const {
  int pk=sqlite3_column_int(stmt, 0);
  const char* id=reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
  CompoundPtr compound(new Compound(id));
  compound->SetOneLetterCode((sqlite3_column_text(stmt, 2))[0]);
  compound->SetChemClass(mol::ChemClass(sqlite3_column_text(stmt, 3)[0]));
  compound->SetDialect(Compound::Dialect(sqlite3_column_text(stmt, 4)[0]));
  const char* f=reinterpret_cast<const char*>(sqlite3_column_text(stmt, 5));
  compound->SetFormula(f);
  compound->SetChemType(mol::ChemType(sqlite3_column_text(stmt, 6)[0]));
  const char* name=reinterpret_cast<const char*>(sqlite3_column_text(stmt, 7));
  compound->SetName(name);
  const char* inchi_code=reinterpret_cast<const char*>(sqlite3_column_text(stmt, 8));
  if (inchi_code) {
    compound->SetInchi(inchi_code);
  }
  const char* inchi_key=reinterpret_cast<const char*>(sqlite3_column_text(stmt, 9));
  if (inchi_key) {
    compound->SetInchiKey(inchi_key);
  }
  int next_column = 10;
  if (smiles_available_) {
    const char* smiles=reinterpret_cast<const char*>(sqlite3_column_text(stmt, next_column));
    next_column++;
    if (smiles) {
      compound->SetSMILES(smiles);
    }
  }
  if (obsolete_available_) {
    bool obsolete=sqlite3_column_int(stmt, next_column);
    compound->SetObsolete(obsolete);
    next_column++;
    const char* replaced_by=reinterpret_cast<const char*>(sqlite3_column_text(stmt, next_column));
    next_column++;
    if (replaced_by) {
      compound->SetReplacedBy(replaced_by);
    }
  }

  // Load atoms and bonds
  this->LoadAtomsFromDB(compound, pk);
  this->LoadBondsFromDB(compound, pk);

  return compound;
}


CompoundPtrList CompoundLib::FindCompounds(const String& query,
                                                    const String& by,
                                                    Compound::Dialect dialect) const {
  CompoundPtrList compounds_vec;

  // Validate "by" argument
  std::set<std::string> allowed_keys{"inchi_code", "inchi_key", "formula"};
  if(smiles_available_) {
     allowed_keys.insert("smiles");
  }
  if (allowed_keys.find(by) == allowed_keys.end()) {
    std::stringstream msg;
    msg << "Invalid 'by' key: " << by;
    throw ost::Error(msg.str());
  }

  String sql_query = BuildFindCompoundQuery(query, dialect, by);

  if(obsolete_available_) {
    // Prefer active compounds, then the ones with a replacement
    sql_query += " ORDER BY obsolete, replaced_by IS NULL";
  }

  // Run the query
  sqlite3_stmt* stmt;
  int retval=sqlite3_prepare_v2(db_->ptr, sql_query.c_str(),
                                static_cast<int>(sql_query.length()),
                                &stmt, NULL);
  sqlite3_bind_text(stmt, 1, query.c_str(),
                      strlen(query.c_str()), NULL);

  if (SQLITE_OK==retval) {
    int ret=sqlite3_step(stmt);
    if (SQLITE_DONE==ret) {
      sqlite3_finalize(stmt);
      return compounds_vec;  // Empty
    }
    while (SQLITE_ROW==ret) {
      CompoundPtr compound = LoadCompoundFromDB(stmt);
      compounds_vec.push_back(compound);
      // next row
      ret=sqlite3_step(stmt);
    }
    assert(SQLITE_DONE==sqlite3_step(stmt));
  } else {
    LOG_ERROR("ERROR: " << sqlite3_errmsg(db_->ptr));
    sqlite3_finalize(stmt);
    return compounds_vec;  // empty
  }
  sqlite3_finalize(stmt);
  return compounds_vec;
}

CompoundPtr CompoundLib::FindCompound(const String& id, 
                                      Compound::Dialect dialect) const {
  // Check cache
  String cache_key = id;
  CompoundMap::const_iterator i=compound_cache_.find(cache_key);
  if (i!=compound_cache_.end()) {
    LOG_DEBUG("Retrieved compound " << cache_key << " from cache");
    return i->second;
  }

  String query = BuildFindCompoundQuery(id, dialect, "tlc");

  // Run the query
  sqlite3_stmt* stmt;
  int retval=sqlite3_prepare_v2(db_->ptr, query.c_str(), 
                                static_cast<int>(query.length()),
                                &stmt, NULL);
  sqlite3_bind_text(stmt, 1, id.c_str(),
                      strlen(id.c_str()), NULL);

  if (SQLITE_OK==retval) {
    int ret=sqlite3_step(stmt);
    if (SQLITE_DONE==ret) {
      sqlite3_finalize(stmt);        
      return CompoundPtr();
    }
    if (SQLITE_ROW==ret) {
      CompoundPtr compound = LoadCompoundFromDB(stmt);
      compound_cache_.insert(std::make_pair(cache_key, compound));
      sqlite3_finalize(stmt);
      return compound;   
    }
    assert(SQLITE_DONE==sqlite3_step(stmt));
  } else {
    LOG_ERROR("ERROR: " << sqlite3_errmsg(db_->ptr));
    sqlite3_finalize(stmt);      
    return CompoundPtr();
  }
  sqlite3_finalize(stmt);    
  return CompoundPtr();
}

CompoundLib::CompoundLib():
  CompoundLibBase(),
  db_(new Database),
  compound_cache_(),
  smiles_available_(),
  obsolete_available_(),
  charges_available_(),
  creation_date_(),
  ost_version_used_() { }

CompoundLib::~CompoundLib() {
  delete db_;
}

}}
