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
#ifndef OST_QUERY_VIEW_WRAPPER_HH
#define OST_QUERY_VIEW_WRAPPER_HH

#include <ost/mol/view_type_fw.hh>
#include <ost/mol/module_config.hh>
#include <ost/mol/entity_view.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/query.hh>

namespace ost { namespace mol {

/*!
  helper class to manage relationship between entity, query and view
*/
class DLLEXPORT_OST_MOL QueryViewWrapper {
public:
  explicit QueryViewWrapper(const EntityHandle& entity_handle);
  explicit QueryViewWrapper(const EntityView& entity_view);
  QueryViewWrapper();
  QueryViewWrapper(const Query& query, const EntityHandle& handle);
  QueryViewWrapper(const Query& query, const EntityView& view = mol::EntityView());
  QueryViewWrapper(const Query& query, QueryFlags flags, const EntityHandle& handle);
  QueryViewWrapper(const Query& query, QueryFlags flags, const EntityView& view = mol::EntityView());

  /// returns true if view will be recalculated from query
  bool DependsOnQuery() const;
  bool IsDataValid() const;

  /// regenerates current view based on internal state
  EntityView GetEntityView() const;

  /// explicitely sets current query
  void SetQuery(const Query& query);
  /// retrieve current query
  const Query& GetQuery() const;

  /// set query flags to be used internally for the Select statement
  void SetFlags(QueryFlags flags);
  /// returns internally used query flags
  QueryFlags GetFlags() const;

  /// return underlying handle
  EntityHandle GetEntity() const {return view_set_ ? entity_view_.GetHandle() : entity_handle_;}

private:
  bool view_set_;
  EntityHandle entity_handle_;
  EntityView entity_view_;
  Query query_;
  QueryFlags flags_;
};

} } // ns

#endif
