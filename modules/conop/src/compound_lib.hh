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
#ifndef OST_CONOP_COMPOUND_LIB_HH
#define OST_CONOP_COMPOUND_LIB_HH

#include <map>
#include <boost/shared_ptr.hpp>

#include "module_config.hh"
#include "compound.hh"
#include "compound_lib_base.hh"

namespace ost { namespace conop {

class CompoundLib;

typedef boost::shared_ptr<CompoundLib> CompoundLibPtr;

class DLLEXPORT_OST_CONOP CompoundLib : public CompoundLibBase {
public:
  static CompoundLibPtr Load(const String& database, bool readonly=true);
  static CompoundLibPtr Create(const String& database);
  ~CompoundLib();
  
  virtual CompoundPtr FindCompound(const String& id, 
                                   Compound::Dialect dialect,
                                   const String& by="tlc") const;
  void AddCompound(const CompoundPtr& compound);
  CompoundLibPtr Copy(const String& filename) const;
  void ClearCache();
  Date GetCreationDate(void);
  String GetOSTVersionUsed(void);
  void SetChemLibInfo(void);
private:
    CompoundLib();

    void LoadAtomsFromDB(CompoundPtr comp, int pk) const;
    void LoadBondsFromDB(CompoundPtr comp, int pk) const;    
private:
  struct Database;
  Database* db_;
  mutable CompoundMap       compound_cache_;
  bool                      smiles_available_; //whether smiles are available in db - introduced in 2.6.0
  bool                      obsolete_available_; //whether obsolete info is available in db - introduced in 2.6.0
  bool                      charges_available_; //whether atom charges are available in db - introduced in 2.6.0
  Date                      creation_date_;
  String                    ost_version_used_;
};

}}

#endif
