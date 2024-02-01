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
#ifndef OST_QUERY_STATE_HH
#define OST_QUERY_STATE_HH

/*
 Author: Marco Biasini
 */

#include <ost/geom/geom.hh>
 
#include <ost/mol/impl/query_impl_fw.hh>
#include <ost/mol/impl/query_ast_fw.hh>

#include <ost/mol/handle_type_fw.hh>
#include <boost/logic/tribool.hpp>
#include <ost/mol/impl/entity_impl_fw.hh>
#include <ost/mol/impl/chain_impl_fw.hh>
#include <ost/mol/impl/residue_impl_fw.hh>
#include <ost/mol/impl/atom_impl_fw.hh>

#include "property_id.hh"

namespace ost { namespace mol {

class EntityHandle;
class EntityView;
struct LazilyBoundData;
struct LazilyBoundRef;

/// \brief Holds state used during query evaluation
/// 
/// The usage of this class is documented  \ref query "here"
class DLLEXPORT_OST_MOL QueryState {
  friend class impl::QueryImpl;
public:
  bool IsAtomSelected(const AtomHandle& atom);

  bool IsResidueSelected(const ResidueHandle& residue);

  bool IsChainSelected(const ChainHandle& chain);  
  
  /// Evaluate chain
  boost::logic::tribool EvalChain(const impl::ChainImplPtr& c);    

  /// Evaluate residue. Assumes that EvalChain has been called
  boost::logic::tribool EvalResidue(const impl::ResidueImplPtr& r);
  
  /// Evaluate atom. Assumes that EvalAtom has been called
  boost::logic::tribool EvalAtom(const impl::AtomImplPtr& a);
  
  void Reset(Prop::Level level);
    
  ~QueryState();
protected:
  /// \internal
  QueryState(const impl::QueryImpl& query, const EntityHandle& ref);
  QueryState(const impl::QueryImpl& query, const EntityView& ref);  
  QueryState();  

  const LazilyBoundRef& GetBoundObject(int i) const;  
  /// Evaluate stack.
  boost::logic::tribool EvalStack(Prop::Level level);

private:
  bool do_within(const geom::Vec3& pos, const impl::WithinParam& p, 
                 impl::CompOP op);
  std::vector<boost::logic::tribool>   s_;
  std::vector<boost::logic::tribool>   value_stack_;
  boost::shared_ptr<LazilyBoundData>   r_;
  const impl::QueryImpl&               q_;
};

}} //ns

#endif
