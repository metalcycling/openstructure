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

  typedef enum {
    PEPTIDE_LINKING   ='P',
    D_PEPTIDE_LINKING ='D',
    L_PEPTIDE_LINKING ='L',
    RNA_LINKING       ='R',
    DNA_LINKING       ='S',
    NON_POLYMER       ='N',
    L_SACCHARIDE      ='X',
    D_SACCHARIDE      ='Y',
    SACCHARIDE        ='Z',
    WATER             ='W',
    UNKNOWN           ='U'
  } Type;
  
  explicit ChemClass(Type chem_class): chem_class_(chem_class) { }

  explicit ChemClass(char type): chem_class_(Type(type)) { }

  ChemClass(): chem_class_(UNKNOWN) { }

  bool operator==(const ChemClass& cc) const {
    return cc.chem_class_ == chem_class_;
  }

  bool operator!=(const ChemClass& cc) const {
    return !this->operator == (cc);
  }

  bool IsPeptideLinking() const {
    return (chem_class_ == PEPTIDE_LINKING ||
            chem_class_ == D_PEPTIDE_LINKING ||
            chem_class_ == L_PEPTIDE_LINKING);
  }
  bool IsNucleotideLinking() const {
    return (chem_class_ == DNA_LINKING || 
            chem_class_ == RNA_LINKING);
  }
  
  bool IsWater() const {
    return chem_class_ == WATER;
  }

  operator char() const {
    return chem_class_;
  }

  bool IsSaccharide() const {
    return (chem_class_ == SACCHARIDE ||
            chem_class_ == L_SACCHARIDE ||
            chem_class_ == D_SACCHARIDE);
  }

private:
  Type chem_class_;
};

}} // ns
#endif
