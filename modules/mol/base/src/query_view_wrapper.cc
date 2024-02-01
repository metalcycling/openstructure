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
#include "query_view_wrapper.hh"

#include <ost/mol/view_type_fw.hh>

namespace ost { namespace mol {

QueryViewWrapper::QueryViewWrapper():
        view_set_(false),
        entity_handle_(),
        entity_view_(),
        query_(),
        flags_(0){}

QueryViewWrapper::QueryViewWrapper(const EntityHandle& entity_handle):
    view_set_(false),
    entity_handle_(entity_handle),
    entity_view_(),
    query_(),
    flags_(0){}

QueryViewWrapper::QueryViewWrapper(const EntityView& entity_view):
    view_set_(true),
    entity_handle_(),
    entity_view_(entity_view),
    query_(),
    flags_(0){}

QueryViewWrapper::QueryViewWrapper(const Query& query, const EntityHandle& entity_handle):
    view_set_(false),
    entity_handle_(entity_handle),
    entity_view_(),
    query_(query),
    flags_(0){}

QueryViewWrapper::QueryViewWrapper(const Query& query, QueryFlags flags,
                                   const EntityHandle& entity_handle):
    view_set_(false),
    entity_handle_(entity_handle),
    entity_view_(),
    query_(query),
    flags_(flags){}

QueryViewWrapper::QueryViewWrapper(const Query& query, const EntityView& view):
    view_set_(true),
    entity_handle_(),
    entity_view_(view),
    query_(query),
    flags_(0){}

QueryViewWrapper::QueryViewWrapper(const Query& query, QueryFlags flags,
                                   const EntityView& view):
    view_set_(true),
    entity_handle_(),
    entity_view_(view),
    query_(query),
    flags_(flags){}

EntityView QueryViewWrapper::GetEntityView() const
{
  if(view_set_) {
    return query_.MatchAll() ? entity_view_ : entity_view_.Select(query_,flags_);
  } else {
    if(entity_handle_.IsValid()){
      return entity_handle_.Select(query_,flags_);
    }
    return entity_view_;
  }
}

const Query& QueryViewWrapper::GetQuery() const
{
  return query_;
}

void QueryViewWrapper::SetQuery(const Query& query)
{
  query_ = query;
}

void QueryViewWrapper::SetFlags(QueryFlags flags){
  flags_ = flags;
}
QueryFlags QueryViewWrapper::GetFlags() const{
  return flags_;
}

bool QueryViewWrapper::IsDataValid() const{
  if(view_set_){
    return entity_view_.IsValid();
  }
  return entity_handle_.IsValid();
}

bool QueryViewWrapper::DependsOnQuery() const
{
  return query_.MatchAll();
}

}} // ns
