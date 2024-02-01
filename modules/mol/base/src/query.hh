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
#ifndef OST_QUERY_HF
#define OST_QUERY_HF

/*
 Author: Marco Biasini
 */

#include <string>
#include <vector>
#include <set>

#include <boost/variant.hpp>
#include <boost/logic/tribool.hpp>

#include <ost/range.hh>

#include <ost/mol/module_config.hh>
#include <ost/mol/impl/query_impl_fw.hh>
#include <ost/mol/query_error.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/query_state.hh>

namespace ost { namespace mol {

/// \brief query flags
/// 
/// Query flags influence the evaluation of Select statements for EntityView
/// and EntityHandle. By default, only inclusive bonds, i.e. bonds whose 
/// involved atoms are both part of the resulting view are included in the view.
/// This can be changed by passing EXCLUSIVE_BONDS or NO_BONDS to the Select
/// call. EXCLUSIVE_BONDS adds bonds as soon as one of the two bonding parters
/// evaluates to true. NO_BONDS suppresses the bond inclusion step. This is 
/// mostly used as a speed optimisation.
/// 
/// If MATCH_RESIDUES is present, the matching behaviour is changed: Instead of
/// including only atoms for which the query evaluates to true, all atoms of
/// residues for which at least one atom matches the query are added to the 
/// resulting view.
///
/// For matching GenericProperties, use gapropname, grpropname, gcpropname for
/// the GenericProp 'propname' at the atom, residue and chain level
/// respectively. The g marks the property as generic.
struct DLLEXPORT_OST_MOL QueryFlag {
  typedef enum {
    EXCLUSIVE_BONDS=1,
    NO_BONDS=2,
    MATCH_RESIDUES=4
  } Flag;
};

typedef uint QueryFlags;

/// \brief Selection Query
///
/// The query system is described \ref query "here" in more detail.
class DLLEXPORT_OST_MOL Query {
public:
  /// Construct new query object.
  Query(const String& query_string="");
  
  /// \brief internally used constructor
  Query(const impl::QueryImplP& impl);
  
  /// \brief whether atom  fulfills the query
  bool IsAtomSelected(const AtomHandle& ah) const;

  /// \brief whether chain fulfills the query
  /// 
  /// True is returned if at least one of the atoms of the chain is 
  /// selected by the query
  bool IsChainSelected(const ChainHandle& chain) const;
  
  /// \brief whether residue fulfills the query
  /// 
  /// True is returned if at least one of the atoms of the residue is 
  /// selected by the query
  bool IsResidueSelected(const ResidueHandle& residue) const;  
  
  //! Get String
  const String& GetQueryString() const;
  
  const QueryErrorDesc& GetErrorDescription() const;

  /// \brief create new query state to be used for a given entity
  QueryState CreateQueryState(const EntityView& view) const;
  
  /// \brief create new query state to be used for a given entity
  QueryState CreateQueryState(const EntityHandle& view) const;
  
  // only valid if parsed without error
  bool IsValid() const;
  // an empty query String defaults to matching everything
  bool MatchAll() const;
  
  // for internal use only
  impl::QueryImplP Impl();
  const impl::QueryImplP& Impl() const;
private:  
  impl::QueryImplP   impl_;
};

// inlined helper function to quote strings for use in queries (e.g. cname=..).
// throws Error if string cannot be quoted
inline String DLLEXPORT_OST_MOL QueryQuoteName(const String& name) {
  // check what quotation marks to use
  char quote = '\'';
  if (name.find('\'') != String::npos) {
    if (name.find('"') != String::npos) {
      throw Error("Cannot quote chain name " + name + " because it contains "
                  "both ' and \" in its name.");
    }
    quote = '"';
  }
  // check problematic \ at end (escapes quotation mark and breaks logic)
  if (name[name.length() - 1] == '\\') {
    throw Error("Cannot quote chain name " + name + " because it ends in \\.");
  }
  return quote + name + quote;
}

}} // ns

#endif
