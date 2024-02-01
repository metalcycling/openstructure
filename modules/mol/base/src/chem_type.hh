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
#ifndef OST_BASE_CHEM_TYPE_HI
#define OST_BASE_CHEM_TYPE_HI
#include <vector>

#include <boost/shared_ptr.hpp>

#include <ost/mol/module_config.hh>


namespace ost { namespace mol {

struct ChemType {
  const static char IONS                  ='I';
  const static char NONCANONICALMOLS      ='M';
  const static char SACCHARIDES           ='S';
  const static char NUCLEOTIDES           ='N';
  const static char AMINOACIDS            ='A';
  const static char COENZYMES             ='E';
  const static char WATERCOORDIONS        ='C';
  const static char DRUGS                 ='D';
  const static char WATERS                ='W';
  const static char UNKNOWN               ='U';

  explicit ChemType(char chem_type)
    : chem_type_(chem_type) {
  }

  ChemType()
    : chem_type_(UNKNOWN) {
  }

  bool operator==(const ChemType& cc) const {
    return cc.chem_type_==chem_type_;
  }

  bool operator!=(const ChemType& cc) const {
    return cc.chem_type_!=chem_type_;
  }

  bool IsIon() const {
    return (chem_type_==ChemType::IONS ||
            chem_type_==ChemType::WATERCOORDIONS);
  }

  bool IsNucleotide() const {
    return (chem_type_==ChemType::NUCLEOTIDES);
  }

  bool IsSaccharide() const {
    return (chem_type_==ChemType::SACCHARIDES);
  }

  bool IsAminoAcid() const {
    return (chem_type_==ChemType::AMINOACIDS);
  }

  bool IsCoenzyme() const {
    return (chem_type_==ChemType::COENZYMES);
  }

  bool IsDrug() const {
    return (chem_type_==ChemType::DRUGS);
  }

  bool IsWater() const {
    return (chem_type_==ChemType::WATERS);
  }

  bool IsNonCanonical() const {
    return (chem_type_==ChemType::NONCANONICALMOLS);
  }

  bool IsKnown() const {
    return (chem_type_!=ChemType::UNKNOWN);
  }

  operator char() const {
    return chem_type_;
  }

private:
  char chem_type_;
};

}} // ns
#endif
