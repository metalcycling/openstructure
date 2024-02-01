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

  const static char PEPTIDE_LINKING   ='P';
  const static char D_PEPTIDE_LINKING ='D';
  const static char L_PEPTIDE_LINKING ='L';
  const static char RNA_LINKING       ='R';
  const static char DNA_LINKING       ='S';
  const static char NON_POLYMER       ='N';
  const static char L_SACCHARIDE      ='X';
  const static char D_SACCHARIDE      ='Y';
  const static char SACCHARIDE        ='Z';
  const static char WATER             ='W';
  const static char UNKNOWN           ='U';
  
  // for backward compatibility to 1.1 and earlier
  const static char PeptideLinking   =PEPTIDE_LINKING;
  const static char DPeptideLinking  =D_PEPTIDE_LINKING;
  const static char LPeptideLinking  =L_PEPTIDE_LINKING;
  const static char RNALinking       =RNA_LINKING;  
  const static char DNALinking       =DNA_LINKING;    
  const static char NonPolymer       =NON_POLYMER;
  const static char LSaccharide      =L_SACCHARIDE;
  const static char DSaccharide      =D_SACCHARIDE;    
  const static char Saccharide       =SACCHARIDE;
  const static char Water            =WATER;
  const static char Unknown          =UNKNOWN;
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
private:
  char chem_class_;
};

}} // ns
#endif
