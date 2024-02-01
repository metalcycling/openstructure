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
#ifndef OST_GFX_NODE_HH
#define OST_GFX_NODE_HH

/*
  Author: Ansgar Philippsen
*/

#include <vector>

#include <boost/enable_shared_from_this.hpp>

#include <ost/gfx/module_config.hh>

#include "render_pass.hh"
#include "gfx_node_fw.hh"
#include "gfx_object_fw.hh"
#include "gfx_node_visitor.hh"
#include "povray_fw.hh"
#include "exporter_fw.hh"

namespace ost { namespace gfx {

typedef std::vector<GfxNodeP> GfxNodeVector;

class DLLEXPORT_OST_GFX GfxNode: public boost::enable_shared_from_this<GfxNode>
{
 public:

  // initialize with a name
  GfxNode(const String& name);

  virtual ~GfxNode();

  // deep copy interface
  virtual GfxNodeP Copy() const;

  virtual void DeepSwap(GfxNode& n);

  // render all child leaves and nodes
  virtual void RenderGL(RenderPass pass);

  // render all child leaves and nodes into POVray state
  virtual void RenderPov(PovState& pov);

  virtual void Export(Exporter* ex);

  // visitor interface
  virtual void Apply(GfxNodeVisitor& v, GfxNodeVisitor::Stack st);

  virtual int GetType() const;

  // return name
  String GetName() const;

  // remove all child nodes
  void RemoveAll();
  // change name
  void Rename(const String& name);

  // add a graphical object - leaf
  void Add(GfxObjP obj);
  
  
  /// \brief returns true if no scene node of the given name is a child
  ///   of this node.
  bool IsNameAvailable(const String& name) const;
  
  
  // remove given graphical object
  void Remove(GfxObjP obj);

  // add another node - branch
  void Add(GfxNodeP node);

  // remove given node
  void Remove(GfxNodeP node);

  // remove obj or node of given name (or several if multiple matches)
  void Remove(const String& name);

  size_t GetChildCount() const;
  
  // hide all child leafs and nodes
  void Hide();
  // display all child leafs and nodes
  void Show();
  // return visibility state
  bool IsVisible() const;

  virtual void ContextSwitch();
  
  /// \brief whether the node (or one of it's parents) has been added to the
  ///    scene
  bool IsAttachedToScene() const;
  
  
  gfx::GfxNodeP GetParent() const;
  
  const GfxNodeVector& GetChildren() const { return node_vector_; }
  GfxNodeVector& GetChildren() { return node_vector_; }
 private:
  GfxNode(const GfxNode& o);
  GfxNode& operator=(const GfxNode&);


  String                   name_;
  bool                     show_;
  GfxNodeVector            node_vector_;
  boost::weak_ptr<GfxNode> parent_;
};

}}

#endif
