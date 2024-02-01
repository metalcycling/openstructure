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
#ifndef OST_CONOP_HH
#define OST_CONOP_HH

#include <map>
#include <ost/conop/module_config.hh>
#include <ost/mol/entity_handle.hh>
#include "processor.hh"
#include "compound_lib.hh"
namespace ost { namespace conop {


typedef enum {
  NO_PEPTIDE_BONDS=1,
} ConopFlag;

  
class DLLEXPORT_OST_CONOP Conopology {
public:
  // singleton
  static Conopology& Instance();

  // returns the default compound library (if any)
  CompoundLibPtr GetDefaultLib() const { return lib_; }
  void SetDefaultLib(const CompoundLibPtr& lib)
  { 
    lib_ = lib; 
  }

  bool IsValidElement(const String& element) const;

private:
  Conopology();
  Conopology(const Conopology&):
    known_elements_(),
    lib_()
  {}
  Conopology& operator=(const Conopology&) {return *this;}

  std::set<String> known_elements_;
  CompoundLibPtr lib_;
};

}} //

#endif
