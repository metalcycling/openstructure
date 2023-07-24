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
#ifndef OST_CONOP_CHEMDICT_PARSER_HH
#define OST_CONOP_CHEMDICT_PARSER_HH

/*
  Author: Marco Biasini
 */
 
 
#include <ost/mol/chem_class.hh>
#include <ost/mol/chem_type.hh>
#include <ost/io/mol/star_parser.hh>
#include <ost/conop/compound_lib.hh>
 
namespace ost { namespace io {

typedef enum {
  ATOM_SPEC,
  BOND_SPEC,
  DESC_SPEC,
  DONT_KNOW
} LoopType;

class DLLEXPORT_OST_IO ChemdictParser : public StarParser {
public:
  ChemdictParser(std::istream& stream, conop::Compound::Dialect dialect): 
    StarParser(stream), compound_(new conop::Compound("UNK")), 
    last_(0), loop_type_(DONT_KNOW), dialect_(dialect)
  {
    this->InitTypeMap();
    this->InitPDBXTypeMap();
  }

  virtual bool OnBeginData(const StringRef& data_name);

  virtual bool OnBeginLoop(const StarLoopDesc& header);

  virtual void OnDataRow(const StarLoopDesc& header, 
                         const std::vector<StringRef>& columns);

  virtual void OnDataItem(const StarDataItem& item);

  virtual void OnEndData();

  void SetCompoundLib(const conop::CompoundLibPtr& lib) 
  {
    lib_=lib;
  }
private:
  void InitTypeMap();  
  void InitPDBXTypeMap();
  conop::CompoundLibPtr                   lib_;
  conop::CompoundPtr                      compound_;
  typedef enum {
    ATOM_NAME=0,
    ALT_ATOM_NAME=1,
    IS_AROMATIC=2,
    ORDINAL=3,
    IS_LEAVING=4,
    ELE=5,
    STEREO_CONF=6,
    ATOM_ID1=0,
    ATOM_ID2=1,
    BOND_ORDER=2,
    DESC_TYPE=0,
    DESC=1,
    PROGRAM=2,
  } PropIndex;
  char                                    last_;
  int                                     indices_[10];
  bool                                    insert_;
  static std::map<String, mol::ChemClass> tm_;  
  static std::map<String, mol::ChemType>  xtm_;
  std::map<String, int>                   atom_map_;
  LoopType                                loop_type_;  
  conop::AtomSpec                         atom_;
  conop::Compound::Dialect                dialect_;
};


}}


#endif
