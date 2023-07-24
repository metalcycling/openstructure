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
#ifndef OST_CONOP_COMPOUND_HH
#define OST_CONOP_COMPOUND_HH

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>
#include <ost/string_ref.hh>
#include <ost/conop/module_config.hh>

#include <ost/mol/chem_class.hh>
#include <ost/mol/chem_type.hh>

namespace ost { namespace conop {

struct DLLEXPORT_OST_CONOP Date {
  Date(int y, int m, int d):
    year(y), month(m), day(d)
  { }
  Date():
    year(1900), month(1), day(1)
  { }
  bool operator<(const Date& date) const
  {
    return year<date.year && month<date.month && day<date.day;
  }
  bool operator==(const Date& date) const
  {
    return year==date.year && month==date.month && day==date.day;
  }  
  
  bool operator!=(const Date& date) const
  {
    return !this->operator==(date);
  }
  
  static Date FromString(const StringRef& str)
  {
    std::vector<StringRef> parts=str.split('-');
    assert(parts.size()==3);
    std::pair<bool, int> year=parts[0].to_int();
    std::pair<bool, int> month=parts[1].to_int();
    std::pair<bool, int> day=parts[2].to_int();
    assert(year.first); assert(month.first); assert(day.first);
    return Date(year.second, month.second, day.second);
  }
  
  String ToString() const;

  int year;
  int month;
  int day;  
};

struct DLLEXPORT_OST_CONOP AtomSpec {
  AtomSpec():
    ordinal(0),
    name(),
    alt_name(),
    element(),
    is_leaving(false),
    is_aromatic(),
    charge(0)
  {
  }
  AtomSpec(int o, const String& n, const String& a, const String& e,
           bool l, bool r, int c=0):
    ordinal(o),
    name(n),
    alt_name(a),
    element(e),
    is_leaving(l),
    is_aromatic(r),
    charge(c)
  {}
  int    ordinal;
  String name;
  String alt_name;
  String element;
  bool   is_leaving;
  bool   is_aromatic;
  int    charge;
  bool operator==(const AtomSpec& rhs) const {
    return ordinal==rhs.ordinal && name==rhs.name && alt_name==rhs.alt_name &&
           element==rhs.element && is_leaving==rhs.is_leaving && 
           is_aromatic==rhs.is_aromatic;

  }
  bool operator!=(const AtomSpec& rhs) const {
    return !this->operator==(rhs);
  }
};

struct DLLEXPORT_OST_CONOP BondSpec {
  BondSpec():
    atom_one(0),
    atom_two(0),
    order(1)
  {
  }
  
  BondSpec(int a, int b, int o): atom_one(a), atom_two(b), order(o) {}
  bool operator==(const BondSpec& rhs) const {
    return atom_one==rhs.atom_one && atom_two==rhs.atom_two;
  }

  bool operator!=(const BondSpec& rhs) const {
    return !this->operator==(rhs);
  }
  int atom_one;
  int atom_two;
  int order;
};

typedef std::vector<AtomSpec> AtomSpecList;
typedef std::vector<BondSpec> BondSpecList;
class Compound;
typedef boost::shared_ptr<Compound> CompoundPtr;


/// \brief Knows about the atoms and bonds of a chemical compounds
class DLLEXPORT_OST_CONOP Compound {
public:
  typedef enum {
    PDB    ='P',
    CHARMM ='C',
    OPLS   ='O',
    AMBER  ='A',
  } Dialect;
  
  Compound(const String& id):
    olc_('?'),
    tlc_(id),
    formula_(),
    name_(),
    inchi_(),
    inchi_key_(),
    smiles_(),
    atom_specs_(),
    bond_specs_(),
    chem_class_(),
    chem_type_(),
    dialect_(Compound::PDB),
    creation_date_(),
    mod_date_()
  {
  }

  /// \brief three-letter code that is unique for every compound 
  const String& GetID() const {
    return tlc_;
  }
  Dialect GetDialect() const { return dialect_; }
  
  String GetDialectAsString() const { 
    switch (dialect_) {
      case CHARMM:
        return "CHARMM";
      case PDB:
        return "PDB";
      case OPLS:
        return "OPLS";
      case AMBER:
        return "AMBER";
      default:
        return "";
    }
  }
  void SetDialect(Dialect dialect) { dialect_=dialect; }
  
  void SetOneLetterCode(char olc) {
    olc_=olc;
  }

  /// \brief one letter code, if available.
  ///
  /// For compounds without a defined one-letter code, the code is set to "?".
  /// One letter codes are ambigous. For instance, both glycine and guanosine
  /// have G as their one-letter code.
  char GetOneLetterCode() const {
    return olc_;
  }

  void SetChemClass(mol::ChemClass chem_class) {
    chem_class_=chem_class;
  }

  mol::ChemClass GetChemClass() const {
    return chem_class_;
  }

  void SetChemType(mol::ChemType chem_type) {
    chem_type_=chem_type;
  }

  /// \brief PDB ligand classification from component dictionary
  ///
  /// The PDB classifies all compounds into 7 categories. This classification
  /// is extracted from the PDB component dictionary (field: pdbx_type)
  mol::ChemType GetChemType() const {
    return chem_type_;
  }

  bool IsPeptideLinking() const {
    return chem_class_.IsPeptideLinking();
  }

  bool IsNucleotideLinking() const {
    return chem_class_.IsNucleotideLinking();
  }

  void AddAtom(const AtomSpec& atom) {
    atom_specs_.push_back(atom);
  }

  void AddBond(const BondSpec& bond) {
    bond_specs_.push_back(bond);
  }

  const AtomSpecList& GetAtomSpecs() const {
    return atom_specs_;
  }
  
  int GetAtomSpecIndex(const String& name) const;
  
  const String& GetName() { return name_; }
  
  void SetName(const String& name) { name_=name; }

  void SetFormula(const String& formula) { formula_=formula; }

  const String& GetFormula() { return formula_; }

  void SetInchi(const String& inchi) { inchi_=inchi; }

  const String& GetInchi() { return inchi_; }

  void SetInchiKey(const String& inchikey) { inchi_key_=inchikey; }

  const String& GetInchiKey() { return inchi_key_; }

  void SetSMILES(const String& smiles) { smiles_=smiles; }

  const String& GetSMILES() { return smiles_; }

  const BondSpecList& GetBondSpecs() const {
    return bond_specs_;
  }
  const Date& GetModificationDate() const 
  { 
    return mod_date_; 
  }
  const Date& GetCreationDate() const 
  { 
    return creation_date_; 
  }  
  
  void SetModificationDate(const Date& mod_date) 
  { 
    mod_date_=mod_date; 
  }
  
  void SetCreationDate(const Date& creation_date) 
  { 
    creation_date_=creation_date; 
  }  
private:
  Compound();
  char                         olc_;
  String                       tlc_;
  String                       formula_;
  String                       name_;
  String                       inchi_;
  String                       inchi_key_;
  String                       smiles_;
  AtomSpecList                 atom_specs_;
  BondSpecList                 bond_specs_;
  mol::ChemClass               chem_class_;
  mol::ChemType                chem_type_;
  Dialect                      dialect_;  
  Date                         creation_date_;
  Date                         mod_date_;
};

typedef std::map<String, CompoundPtr> CompoundMap;

}}

#endif
