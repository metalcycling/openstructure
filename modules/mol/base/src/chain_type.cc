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

#include <ost/message.hh>
#include "chain_type.hh"

namespace ost { namespace mol {

ChainType ChainTypeFromString(StringRef identifier)
{

  // chain types as found in the entity category of a mmcif file
  if (StringRef("polymer", 7) == identifier) {
    return CHAINTYPE_POLY;
  } else if (StringRef("non-polymer", 11) == identifier) {
    return CHAINTYPE_NON_POLY;
  } else if (StringRef("water", 5) == identifier) {
    return CHAINTYPE_WATER;
  } else if (StringRef("macrolide", 9) == identifier) {
    return CHAINTYPE_MACROLIDE;
  } else if (StringRef("branched", 8) == identifier) {
    return CHAINTYPE_BRANCHED;
  // chain types as found in the entity_poly category of a mmcif file
  } else if (StringRef("polypeptide(D)", 14) == identifier) {
    return CHAINTYPE_POLY_PEPTIDE_D;
  } else if (StringRef("polypeptide(L)", 14) == identifier) {
    return CHAINTYPE_POLY_PEPTIDE_L;
  } else if (StringRef("polydeoxyribonucleotide", 23) == identifier) {
    return CHAINTYPE_POLY_DN;
  } else if (StringRef("polyribonucleotide", 18) == identifier) {
    return CHAINTYPE_POLY_RN;
  } else if (StringRef("polysaccharide(D)", 17) == identifier) {
    return CHAINTYPE_POLY_SAC_D;
  } else if (StringRef("polysaccharide(L)", 17) == identifier) {
    return CHAINTYPE_POLY_SAC_L;
  } else if (StringRef("polydeoxyribonucleotide/polyribonucleotide hybrid",
                        49) == identifier) {
    return CHAINTYPE_POLY_DN_RN;
  } else if (StringRef("cyclic-pseudo-peptide", 21) == identifier) {
    return CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE;
  } else if (StringRef("peptide nucleic acid", 20) == identifier) {
    return CHAINTYPE_POLY_PEPTIDE_DN_RN;
  } else if (StringRef("oligosaccharide", 15) == identifier) {
    return CHAINTYPE_OLIGOSACCHARIDE;
  } else if (StringRef("other", 5) == identifier) {
    // According to the mmCIF dictionary, "other" only exists in
    // _entity_poly.type. Therefore, "other" can only be a generic polymer.
    return CHAINTYPE_POLY;
  }

  throw Error("Unrecognised chain type descriptor found: '" +
                           identifier.str() +"'!");
}

ChainType ChainTypeFromString(const String& identifier)
{
  return ChainTypeFromString(StringRef(identifier.c_str(),
                                       identifier.length()));
}

String StringFromChainType(ChainType type)
{
  // chain types as found in the entity category of a mmcif file
  if (CHAINTYPE_POLY == type) {
    // "polymer" in _entity.type
    // "other" in _entity_poly.type
    return "polymer";
  } else if (CHAINTYPE_NON_POLY == type) {
    return "non-polymer";
  } else if (CHAINTYPE_WATER == type) {
    return "water";
  } else if (CHAINTYPE_MACROLIDE == type) {
    return "macrolide";
  } else if (CHAINTYPE_BRANCHED == type) {
    return "branched";
  // chain types as found in the entity_poly category of a mmcif file
  } else if (CHAINTYPE_POLY_PEPTIDE_D == type) {
    return "polypeptide(D)";
  } else if (CHAINTYPE_POLY_PEPTIDE_L == type) {
    return "polypeptide(L)";
  } else if (CHAINTYPE_POLY_DN == type) {
    return "polydeoxyribonucleotide";
  } else if (CHAINTYPE_POLY_RN == type) {
    return "polyribonucleotide";
  } else if (CHAINTYPE_POLY_SAC_D == type) {
    return "polysaccharide(D)";
  } else if (CHAINTYPE_POLY_SAC_L == type) {
    return "polysaccharide(L)";
  } else if (CHAINTYPE_POLY_DN_RN == type) {
    return "polydeoxyribonucleotide/polyribonucleotide hybrid";
  } else if (CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE == type) {
    return "cyclic-pseudo-peptide";
  } else if (CHAINTYPE_POLY_PEPTIDE_DN_RN == type) {
    return "peptide nucleic acid";
    // chain types as found in the pdbx_entity_branch category of a mmcif file
  } else if (CHAINTYPE_OLIGOSACCHARIDE == type) {
    return "oligosaccharide";
    // other...
  } else if (CHAINTYPE_UNKNOWN == type) {
    return "other";
  }

  std::stringstream ss("Unknown ChainType item found: '");
  ss << type << "'!";
  throw Error(ss.str());
}

String EntityTypeFromChainType(ChainType type) {
  switch(type) {
     case ost::mol::CHAINTYPE_POLY: return "polymer";
     case ost::mol::CHAINTYPE_NON_POLY: return "non-polymer";
     case ost::mol::CHAINTYPE_WATER: return "water";
     case ost::mol::CHAINTYPE_POLY_PEPTIDE_D: return "polymer";
     case ost::mol::CHAINTYPE_POLY_PEPTIDE_L: return "polymer";
     case ost::mol::CHAINTYPE_POLY_DN: return "polymer";
     case ost::mol::CHAINTYPE_POLY_RN: return "polymer";
     case ost::mol::CHAINTYPE_POLY_SAC_D: return "polymer";
     case ost::mol::CHAINTYPE_POLY_SAC_L: return "polymer";
     case ost::mol::CHAINTYPE_POLY_DN_RN: return "polymer";
     case ost::mol::CHAINTYPE_MACROLIDE: return "macrolide";         
     case ost::mol::CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE: return "polymer";
     case ost::mol::CHAINTYPE_POLY_PEPTIDE_DN_RN: return "polymer";
     case ost::mol::CHAINTYPE_BRANCHED: return "branched";
     case ost::mol::CHAINTYPE_OLIGOSACCHARIDE: return "branched";
     default: {
       std::stringstream ss;
       ss <<"Unknown ChainType item found: '" << type << "'!";
       throw Error(ss.str());
     }
  }

}

String EntityPolyTypeFromChainType(ChainType type) {
  switch(type) {
     case ost::mol::CHAINTYPE_POLY: return "other";
     case ost::mol::CHAINTYPE_POLY_PEPTIDE_D: return "polypeptide(D)";
     case ost::mol::CHAINTYPE_POLY_PEPTIDE_L: return "polypeptide(L)";
     case ost::mol::CHAINTYPE_POLY_DN: return "polydeoxyribonucleotide";
     case ost::mol::CHAINTYPE_POLY_RN: return "polyribonucleotide";
     case ost::mol::CHAINTYPE_POLY_SAC_D: return "other"; // older dictionaries have "polysaccharide(D)"
     case ost::mol::CHAINTYPE_POLY_SAC_L: return "other"; // older dictionaries have "polysaccharide(L)"
     case ost::mol::CHAINTYPE_POLY_DN_RN: return "polydeoxyribonucleotide/polyribonucleotide hybrid";
     case ost::mol::CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE: return "cyclic-pseudo-peptide";
     case ost::mol::CHAINTYPE_POLY_PEPTIDE_DN_RN: return "peptide nucleic acid";
     default: {
       std::stringstream ss;
       ss << "Cannot return entity poly type from chain of type: '" << type << "'!";
       throw Error(ss.str());
     }
  }
}

String BranchedTypeFromChainType(ChainType type) {
  switch(type) {
     case ost::mol::CHAINTYPE_BRANCHED: return "oligosaccharide"; // the only one
     case ost::mol::CHAINTYPE_OLIGOSACCHARIDE: return "oligosaccharide";
     default: {
       std::stringstream ss;
       ss << "Cannot return branched type from chain of type: '" << type << "'!";
       throw Error(ss.str());
     }
  }

}

}} //ns
