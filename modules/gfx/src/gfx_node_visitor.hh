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
#ifndef OST_GFX_NODE_VISITOR_HH
#define OST_GFX_NODE_VISITOR_HH

/*
  Author: Ansgar Philippsen
*/

#include <stack>

#include <ost/gfx/module_config.hh>
#include "gfx_node_fw.hh"
#include "gfx_object_fw.hh"

namespace ost { namespace gfx {

class DLLEXPORT_OST_GFX GfxNodeVisitor
{
public:
  // in the below methods, contains the 'path' to the current 
  // node or obj, where top() returns the direct parent
  typedef std::stack<GfxNode*> Stack;

  virtual void OnStart() {}
  virtual void OnEnd() {}

  // is called for each base node in the hierarchy
  // expected to return true if descent into child nodes is wanted
  virtual bool VisitNode(GfxNode* node, const Stack& st) {return true;}

  // is called for each leaf
  virtual void VisitObject(GfxObj* obj, const Stack& st) {}

  // is called just before returning to parent hierarchy
  virtual void LeaveNode(GfxNode* node, const Stack& st) {}
  
  virtual ~GfxNodeVisitor() {}
};

}} // ns

#endif
