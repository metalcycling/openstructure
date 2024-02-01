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
#ifndef OST_GFX_SCENE_OBSERVER_HH
#define OST_GFX_SCENE_OBSERVER_HH

/*
  Author: Ansgar Philippsen, Marco Biasini
*/

#include <boost/shared_ptr.hpp>

#include <ost/mol/view_type_fw.hh>

#include "gfx_node_fw.hh"
#include "gfx_object_fw.hh"

namespace ost { namespace gfx {

/// \brief Interface for observing graphical objects
/// 
/// Similar to the EntityObserver, SceneObserver defines an interface for 
/// objects interested in changes of the Scene. The interface provides callbacks
/// for different state changes. All these callbacks are invoked after the state
/// has changed.
class SceneObserver {
public:
  virtual ~SceneObserver() {}

  /// \brief Called when a new node has been added to the scene.
  virtual void NodeAdded(const GfxNodeP& node) {}

  /// \brief Called when an object is removed from the scene.
  virtual void NodeRemoved(const GfxNodeP& node) {}
  
  /// \brief Called when an object is updated
  virtual void NodeChanged(const GfxNodeP& node) {}
  
  /// \brief Called when the selection mode changed
  virtual void SelectionModeChanged(int sel_mode) {}
  
  /// \brief Called when selection of an object changed
  virtual void SelectionChanged(const GfxObjP& obj, 
                                const mol::EntityView& sel) {}

  /// \brief Called when a nodes transform changes
  virtual void NodeTransformed(const GfxNodeP& node) {}

  /// \brief Called when the render mode of a node changed
  virtual void RenderModeChanged(const GfxNodeP& node) {}
};

}}

#endif
