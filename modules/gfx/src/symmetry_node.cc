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

#include "gl_helper.hh"
#include "symmetry_node.hh"

namespace ost { namespace gfx {

SymmetryOp::SymmetryOp(const geom::Mat3& rot, const geom::Vec3& trans):
  tf_(rot)
{
  geom::Mat4 t;
  t.PasteTranslation(trans);
  tf_=Transpose(tf_*t);
}

SymmetryOp::SymmetryOp(const geom::Mat4& transform):
  tf_(transform)
{
}

void SymmetryOp::Push()
{
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glMultMatrix(tf_.Data());
}

void SymmetryOp::Pop()
{
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

SymmetryNode::SymmetryNode(const String& name, const SymmetryOpList& sym_ops):
  GfxNode(name), sym_ops_(sym_ops)
{ }

void SymmetryNode::RenderGL(RenderPass pass)
{
  for (SymmetryOpList::iterator i=sym_ops_.begin(),
       e=sym_ops_.end(); i!=e; ++i) {
     SymmetryOp& op=*i;
     op.Push();
     GfxNode::RenderGL(pass);
     op.Pop();
  }
}


}}
