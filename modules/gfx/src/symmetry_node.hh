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
#ifndef OST_GFX_SYMMETRY_NODE_HH
#define OST_GFX_SYMMETRY_NODE_HH

#include <ost/geom/geom.hh>
#include <ost/gfx/module_config.hh>
#include <ost/gfx/gfx_node.hh>

namespace ost { namespace gfx {

// temporary symmetry op class. This will be moved to base (or geom).
class  DLLEXPORT_OST_GFX SymmetryOp {
public:
  SymmetryOp(const geom::Mat3& rot,
             const geom::Vec3& trans=geom::Vec3(0, 0, 0));
  SymmetryOp(const geom::Mat4& transform);
  // push onto GL_MODELVIEW stack
  void Push();
  // pop from GL_MODELVIEW stack
  void Pop();

  bool operator==(const SymmetryOp& rhs) const { return tf_==rhs.tf_; }
  bool operator!=(const SymmetryOp& rhs) const { return tf_!=rhs.tf_; }
private:
  geom::Mat4 tf_;
};

typedef std::vector<SymmetryOp> SymmetryOpList;

/// \brief renders all child nodes and the symmetry related copies
class DLLEXPORT_OST_GFX SymmetryNode : public GfxNode {
public:
  SymmetryNode(const String& name, const SymmetryOpList& sym_ops);

  virtual void RenderGL(RenderPass pass);
private:
  SymmetryOpList  sym_ops_;
};

///  \example gfx_symmetry.py
///
/// Uses the gfx::SymmetryNode class to draw the symmetry equivalents of a 
/// small protein fragment. 
/// 
/// By using rigid body manipulator, the relative orientation of the fragments
/// can be adjusted.
}}

#endif
