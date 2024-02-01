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
#ifndef OST_CONOP_AMINO_ACIDS_HH
#define OST_CONOP_AMINO_ACIDS_HH
/*
  Author: Marco Biasini
 */
#include <set>
#include <ost/dllexport.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/torsion_handle.hh>

#include <ost/conop/module_config.hh>

namespace ost { namespace conop {

// convenience enum for standard set of 20 amino acids
// Xxx signifies an unknown amino acid.
typedef enum  {
  ALA, ARG, ASN,
  ASP, GLN, GLU,
  LYS, SER, CYS,
  MET, TRP, TYR,
  THR, VAL, ILE,
  LEU, GLY, PRO,
  HIS, PHE, XXX
} AminoAcid;

/// \brief from residue name to amino acid.
/// 
/// for non standard amino acids Xxx is returned.
DLLEXPORT_OST_CONOP AminoAcid ResidueToAminoAcid(const mol::ResidueHandle& r);

// \brief from amino acid to residue name
DLLEXPORT_OST_CONOP String AminoAcidToResidueName(AminoAcid aa);

DLLEXPORT_OST_CONOP String OneLetterCodeToResidueName(char olc);

DLLEXPORT_OST_CONOP AminoAcid OneLetterCodeToAminoAcid(char olc);

char DLLEXPORT_OST_CONOP ResidueNameToOneLetterCode(String rn);

AminoAcid DLLEXPORT_OST_CONOP ResidueNameToAminoAcid(String rn);

class AminoAcidSetIterator : public std::iterator<std::forward_iterator_tag,
                                                  AminoAcid> {
public:
  AminoAcidSetIterator(unsigned int bits, int start):
    bits_(bits), curr_(start)
  {
  }
  AminoAcidSetIterator& operator++()
  {
    this->Advance();
    return *this;
  }
    
  AminoAcid operator*() const 
  {
    return AminoAcid(curr_);
  }
  
  bool operator==(const AminoAcidSetIterator& rhs) const
  {
    return curr_==rhs.curr_;
  }
  bool operator!=(const AminoAcidSetIterator& rhs) const
  {
    return !this->operator==(rhs);
  }
private:
  void Advance()
  {
    ++curr_;
    while (curr_<=XXX && !(bits_ & (1 << curr_))) { ++curr_; }
  }
  unsigned int bits_;
  int          curr_;
};

/// \brief Amino acid bit set
class DLLEXPORT_OST_CONOP AminoAcidSet {
public:
  typedef AminoAcidSetIterator Iterator;
  
  static AminoAcidSet CreatePolarSet();
  
  static AminoAcidSet CreateAromaticSet();
  
  static AminoAcidSet CreateApolarSet();
  
  static AminoAcidSet CreateSet(AminoAcid aa);

  static std::vector<AminoAcidSet> CreateCompleteSet();

  static std::vector<AminoAcidSet> CreateThreeStateSet();

  static std::vector<AminoAcidSet> CreatePseudoSet();

  AminoAcidSet();
  
  /// \brief add another amino acid
  void Add(AminoAcid amino_acid);
  
  //// \brief remove an amino acid
  void Remove(AminoAcid amino_acid);
  
  /// \brief whether the set contains a given amino acid
  bool Contains(AminoAcid amino_acid) const;
  
  /// \brief whether the set is empty
  bool Empty() const;
  
  /// \brief iterator pointing to the beginning of the amino acid set
  /// \code 
  /// for (AminoAcidSet::iterator i=aa_set.Begin(), e=aa_set.End(); i!=e; ++i) {
  ///   . . .
  /// }
  /// \endcode
  /// will iterate over all the amino acids present in the set.
  Iterator Begin() const;
  
  /// \brief iterator pointing to the end of the amino acids.
  Iterator End() const;
  
  bool operator==(const AminoAcidSet& rhs) const;
  
  bool operator!=(const AminoAcidSet& rhs) const;
  
  template <typename DS>
  void Serialize(DS& ds)
  {
    ds & bits_;
  }
private:
  unsigned int bits_;
};

DLLEXPORT_OST_CONOP std::ostream& operator<<(std::ostream& os, 
                         const AminoAcidSet& aa_set);

typedef std::vector<AminoAcidSet> AminoAcidAlphabet;

}}

#endif
