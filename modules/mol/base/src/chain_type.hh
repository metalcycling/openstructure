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
#ifndef OST_CHAIN_TYPE_HH
#define OST_CHAIN_TYPE_HH

#include <ost/mol/module_config.hh>
#include <ost/base.hh>
#include <ost/string_ref.hh>

#include "module_config.hh"

namespace ost { namespace mol {

/// \enum different kinds of chains
///
/// Warning: this class mixes vocabulary from _entity.type and
// _entity_poly.type, which is more detailed. As a result it is not a 1:1
// mapping and cannot be used to to read/write mmCIF entity types accurately.

typedef enum {
  CHAINTYPE_POLY,           ///< polymer
  CHAINTYPE_NON_POLY,       ///< non-polymer
  CHAINTYPE_WATER,          ///< water
  CHAINTYPE_POLY_PEPTIDE_D, ///< (D) amino acid sequence
  CHAINTYPE_POLY_PEPTIDE_L, ///< (L) amino acid sequence
  CHAINTYPE_POLY_DN,        ///< polydeoxyribonucleotide
  CHAINTYPE_POLY_RN,        ///< polyribonucleotide
  CHAINTYPE_POLY_SAC_D,     ///< polysaccharide(D)
  CHAINTYPE_POLY_SAC_L,     ///< polysaccharide(L)
  CHAINTYPE_POLY_DN_RN,     ///< polydeoxyribonucleotide/ -ribonucleotide hybrid
  CHAINTYPE_UNKNOWN,        ///< guess what
  // new chain types
  CHAINTYPE_MACROLIDE,              ///< macrolide
  CHAINTYPE_CYCLIC_PSEUDO_PEPTIDE,  ///< cyclic-pseudo-peptide
  CHAINTYPE_POLY_PEPTIDE_DN_RN,     ///< peptide nucleic acid
  CHAINTYPE_BRANCHED,               ///< carbohydrate
  CHAINTYPE_OLIGOSACCHARIDE,        ///< oligosaccharide (branched carbohydrate,
                                    ///< i.e. _entity.type is strictly 'branched')
  CHAINTYPE_N_CHAINTYPES    ///< no. of chain types
} ChainType;

/// \brief Create a ChainType item for a given string
///
/// \param identifier StringRef to be translated
///
/// \return The ChainType corresponding to the input, throws a
///         ost:Error on unknown type
ChainType DLLEXPORT_OST_MOL ChainTypeFromString(const StringRef identifier);

/// \brief Create a ChainType item for a given string
///
/// \param identifier String to be translated
///
/// \return The ChainType corresponding to the input, throws a
///         ost::Error on unknown type
ChainType DLLEXPORT_OST_MOL ChainTypeFromString(const String& identifier);

/// \brief Return the String identifier for a given type
///
/// \param type ChainType to be translated
///
/// \return String corresponding to the input, throws a ost::Error on
///         unknown type
String DLLEXPORT_OST_MOL StringFromChainType(ChainType type);

/// \brief Return _entity.type consistent with respective mmCIF vocabulary
///        (mmcif_pdbx_v50):
///        - branched
///        - macrolide
///        - non-polymer
///        - polymer
///        - water
///
///        For consistency with older vocabularies, CHAINTYPE_POLY_SAC_D
///        and CHAINTYPE_POLY_SAC_L return "polymer"
///
/// \param type ChainType to be translated
///
/// \return String corresponding to the input, throws a ost::Error on
///         unknown type
String DLLEXPORT_OST_MOL EntityTypeFromChainType(ChainType type);

/// \brief Return _entity_poly.type consistent with mmCIF dictionary
///        (mmcif_pdbx_v50):
///        - cyclic-pseudo-peptide 	
///        - other 	
///        - peptide nucleic acid 	
///        - polydeoxyribonucleotide 	
///        - polydeoxyribonucleotide/polyribonucleotide hybrid 	
///        - polypeptide(D) 	
///        - polypeptide(L) 	
///        - polyribonucleotide
///
///        For consistency with older dictionaries, CHAINTYPE_POLY_SAC_D
///        and CHAINTYPE_POLY_SAC_L are still accepted but return "other".
///        Older dictionaries still had "polysaccharide(D)" and
///        "polysaccharide(L)""
///
/// \param type ChainType to be translated
///
/// \return String corresponding to the input, throws a ost::Error on
///         unknown type or if it's not of _entity.type polymer
String DLLEXPORT_OST_MOL EntityPolyTypeFromChainType(ChainType type);

/// \brief Return pdbx_entity_branch.type consistent with mmCIF dictionary
///        (mmcif_pdbx_v50):
///        - oligosaccharide	
///
/// \param type ChainType to be translated
///
/// \return String corresponding to the input, throws a ost::Error on
///         unknown type or if it's not of _entity.type branched
String DLLEXPORT_OST_MOL BranchedTypeFromChainType(ChainType type);

}} //ns

#endif
