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
#ifndef OST_MOL_ALG_DISTANCE_TEST_COMMON_HH
#define OST_MOL_ALG_DISTANCE_TEST_COMMON_HH

#include <ost/mol/entity_view.hh>
#include <ost/mol/alg/module_config.hh>

namespace ost { namespace mol { namespace alg {

/// \brief Returns the name of the symmetry equivalent atom in residues with symmetric side-chains
///
/// If the atom does belongs to a residue with a symmetric side-chain and if the atom has a symmetry
/// equivalent, the function returns the name of the symmetry equivalent atom, otherwise it returns
///  the name of the original atom
String DLLEXPORT_OST_MOL_ALG SwappedName(const String& name);

/// \brief Checks if an atom in a residue has a symmetry equivalent
///
/// Returns true if the atom belongs to a residue with a symmetric side-chain and a symmetry equivalent
/// atom exists. Returns false otherwise 
bool DLLEXPORT_OST_MOL_ALG Swappable(const String& rname, const String& aname);

/// \brief Contains the infomation needed to uniquely identify an atom in a structure
///
/// Used by the the Local Distance Difference Test classes and functions
class DLLEXPORT_OST_MOL_ALG UniqueAtomIdentifier
{

public:
  /// \brief Constructor with all the relevant information
  UniqueAtomIdentifier(const String& chain,const ResNum& residue,const String& residue_name, const String& atom): chain_(chain),residue_(residue),residue_name_(residue_name),atom_(atom) {}

  // to make the compiler happy (boost python map suite)
  UniqueAtomIdentifier(): chain_(""),residue_(ResNum(1)),residue_name_(""),atom_("") {}

  /// \brief Returns the name of the chain to which the atom belongs, as a String
  String GetChainName() const { return chain_; }

  /// \brief Returns the ResNum of the residue to which the atom belongs
  ResNum GetResNum() const { return residue_; }

  /// \brief Returns the name of the residue to which the atom belongs, as a String
  String GetResidueName() const { return residue_name_; }

  /// \brief Returns the name of the atom, as a String
  String GetAtomName() const { return atom_; }

  /// \brief Returns the qualified name of the atom, as a String
  String GetQualifiedAtomName() const ;

  // required because UniqueAtomIdentifier is used as a key for a std::map
  bool operator==(const UniqueAtomIdentifier& rhs) const ;

  // required because UniqueAtomIdentifier is used as a key for a std::map
  bool operator<(const UniqueAtomIdentifier& rhs) const;

private:

  String chain_;
  ResNum residue_;
  String residue_name_;
  String atom_;
};

// typedef used to make the code cleaner
typedef std::pair<UniqueAtomIdentifier,UniqueAtomIdentifier> UAtomIdentifiers;

/// \brief Residue distance list.
///
/// Container for all the interatomic distances that are checked in a Local Distance Difference Test
/// and are originating from a single specific residue
typedef std::map<std::pair<UniqueAtomIdentifier,UniqueAtomIdentifier>,std::pair<Real,Real> > ResidueRDMap;

/// \brief Global distance list.
///
/// Container for all the residue-based interatomic distance lists that are checked in a Local Distance Difference Test
/// and  belong to the same structure
typedef std::map<ost::mol::ResNum,ResidueRDMap> GlobalRDMap;

// used by the multi-reference distance-list generator function
typedef std::map<UniqueAtomIdentifier,int> ExistenceMap;

}}}

#endif
