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
/*
 Author: Marco Biasini
 */

#include <vector>
#include <algorithm>

#include <boost/logic/tribool.hpp>

#include "query.hh"
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/chain_impl.hh>
#include <ost/mol/impl/query_impl.hh>

#include "atom_handle.hh"
#include "residue_handle.hh"
#include "chain_handle.hh"

namespace ost { namespace mol {

Query::Query(const String& query_string):
  impl_(impl::QueryImplP(new impl::QueryImpl(query_string))) {
}

Query::Query(const impl::QueryImplP& impl):
  impl_(impl) {
}

bool Query::IsAtomSelected(const AtomHandle& atom) const
{
  QueryState qs = impl_->CreateQueryState(atom.GetEntity());
  return qs.IsAtomSelected(atom);
}

bool Query::IsResidueSelected(const ResidueHandle& residue) const
{
  QueryState qs = impl_->CreateQueryState(residue.GetEntity());
  return qs.IsResidueSelected(residue);
}

bool Query::IsChainSelected(const ChainHandle& chain) const
{
  QueryState qs = impl_->CreateQueryState(chain.GetEntity());
  return qs.IsChainSelected(chain);
}

const String& Query::GetQueryString() const 
{
  return impl_->GetQueryString();
}

const QueryErrorDesc& Query::GetErrorDescription() const {
  return impl_->GetErrorDescription();
}

bool Query::IsValid() const {
  return impl_->IsValid();
}

bool Query::MatchAll() const {
  return impl_->IsEmpty();
}

impl::QueryImplP Query::Impl()
{
  return impl_;
}

const impl::QueryImplP& Query::Impl() const 
{
  return impl_;
}


QueryState Query::CreateQueryState(const EntityView& view) const
{
  return Impl()->CreateQueryState(view);
}

QueryState Query::CreateQueryState(const EntityHandle& handle) const
{
  return Impl()->CreateQueryState(handle);
}

}} // ns
