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
#ifndef OST_BASE_CHEM_CLASS_HH
#define OST_BASE_CHEM_CLASS_HH

#include <ost/mol/module_config.hh>


namespace ost { namespace mol {

struct DLLEXPORT ChemClass {

  static constexpr char PEPTIDE_LINKING   ='P';
  static constexpr char D_PEPTIDE_LINKING ='D';
  static constexpr char L_PEPTIDE_LINKING ='L';
  static constexpr char RNA_LINKING       ='R';
  static constexpr char DNA_LINKING       ='S';
  static constexpr char NON_POLYMER       ='N';
  static constexpr char L_SACCHARIDE      ='X';
  static constexpr char D_SACCHARIDE      ='Y';
  static constexpr char SACCHARIDE        ='Z';
  static constexpr char WATER             ='W';
  static constexpr char UNKNOWN           ='U';
  
  // for backward compatibility to 1.1 and earlier
  static constexpr char PeptideLinking   =PEPTIDE_LINKING;
  static constexpr char DPeptideLinking  =D_PEPTIDE_LINKING;
  static constexpr char LPeptideLinking  =L_PEPTIDE_LINKING;
  static constexpr char RNALinking       =RNA_LINKING;  
  static constexpr char DNALinking       =DNA_LINKING;    
  static constexpr char NonPolymer       =NON_POLYMER;
  static constexpr char LSaccharide      =L_SACCHARIDE;
  static constexpr char DSaccharide      =D_SACCHARIDE;    
  static constexpr char Saccharide       =SACCHARIDE;
  static constexpr char Water            =WATER;
  static constexpr char Unknown          =UNKNOWN;
  explicit ChemClass(char chem_class)
    : chem_class_(chem_class) {
  }

  ChemClass()
    : chem_class_(UNKNOWN) {
  }
  bool operator==(const ChemClass& cc) const {
    return cc.chem_class_==chem_class_;
  }

  bool operator!=(const ChemClass& cc) const {
    return !this->operator==(cc);
  }

  bool IsPeptideLinking() const {
    return (chem_class_==ChemClass::PEPTIDE_LINKING ||
            chem_class_==ChemClass::D_PEPTIDE_LINKING ||
            chem_class_==ChemClass::L_PEPTIDE_LINKING);
  }
  bool IsNucleotideLinking() const {
    return (chem_class_==ChemClass::DNA_LINKING || 
            chem_class_==ChemClass::RNA_LINKING);
  }
  
  bool IsWater() const { return chem_class_==ChemClass::WATER; }
  operator char() const {
    return chem_class_;
  }

  bool IsSaccharide() const {
    return (chem_class_==ChemClass::SACCHARIDE ||
            chem_class_==ChemClass::L_SACCHARIDE ||
            chem_class_==ChemClass::D_SACCHARIDE);
  }

private:
  char chem_class_;
};

}} // ns
#endif
