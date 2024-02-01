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
#ifndef OST_GFX_IMPL_BACKBONE_TRACE_HH
#define OST_GFX_IMPL_BACKBONE_TRACE_HH

/*
  Authors: Ansgar Philippsen, Marco Biasini
 */

#include <ost/mol/entity_view.hh>

#include <ost/gfx/module_config.hh>
#include <ost/gfx/impl/entity_detail.hh>

namespace ost { namespace gfx { namespace impl {

class BackboneTraceBuilder;

/// \internal
/// \brief a smoothed backbone trace
/// 
/// This class can be used to calculate a smoothed backbone trace for an 
/// entity view.
class DLLEXPORT_OST_GFX BackboneTrace {
public:

  // empty trace
  BackboneTrace();

  // initialize with a view, and build
  BackboneTrace(const mol::EntityView& ent);

  // number of node-lists
  int GetListCount() const;

  // grab a list
  const NodeEntryList& GetList(int index) const;
  // grab a list
  NodeEntryList& GetList(int index);

  // reset the view and rebuild
  void ResetView(const mol::EntityView& ent);

  // used internally - adds a finished nodelist
  void AddNodeEntryList(const NodeEntryList& entries);  
  
  // used internally - calculates some derived values for a nodelist
  void PrepList(NodeEntryList& nelist) const;

  // re-creates internal nodelist-list based on view
  /*
    seq_hack will apply an additional hackish N/N+1 rnum check
    to determine if two consecutive residues are connected
  */
  void Rebuild();

  // entity has new positions
  void OnUpdatedPositions();
  
  // extract portions of this backbone trace for a subview
  // this is faster then re-generating a trace
  BackboneTrace CreateSubset(const mol::EntityView& subview);

  // experimental
  void SetSeqHack(bool f);
  bool GetSeqHack() const {return seq_hack_;}
  
  // experimental
  void SetTwistHack(bool f);
  bool GetTwistHack() const {return twist_hack_;}

private:  
  mol::EntityView      view_;
  NodeEntryListList    node_list_list_;
  bool seq_hack_;
  bool twist_hack_;

};

}}}

#endif
