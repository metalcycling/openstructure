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
#include "xcs_editor.hh"
#include <ost/mol/atom_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/bond_handle.hh>

#include "impl/atom_impl.hh"
#include "impl/entity_impl.hh"

#include <ost/invalid_handle.hh>

/*
  Author: Marco Biasini
 */
namespace ost { namespace mol {


XCSEditor::XCSEditor():
  EditorBase(EntityHandle(), UNBUFFERED_EDIT)
{
}

XCSEditor::XCSEditor(const EntityHandle& ent, EditMode mode):
  EditorBase(ent, mode)
{
  if (ent.IsValid()) {
    ent.Impl()->IncXCSEditorCount();
  }
}
  
XCSEditor::~XCSEditor()
{
  if (ent_.IsValid()) {
    ent_.Impl()->DecXCSEditorCount();
  }
}

XCSEditor::XCSEditor(const XCSEditor& rhs):
  EditorBase(rhs.ent_, rhs.mode_)
{
  if (ent_.IsValid()) {
    ent_.Impl()->IncXCSEditorCount();
  }  
}

XCSEditor& XCSEditor::operator=(const XCSEditor& rhs)
{
  ent_=rhs.ent_;
  if (ent_.IsValid()) {
    ent_.Impl()->IncXCSEditorCount();
  }
  return *this;
}

void XCSEditor::SetAtomTransformedPos(const AtomHandle& atom,
                                      const geom::Vec3& position)
{
  CheckHandleValidity(atom);
  impl::EntityImplPtr eip=ent_.Impl();
  atom.Impl()->TransformedPos()=position;
  atom.Impl()->OriginalPos() = eip->HasTransform() ? eip->GetTransform().ApplyInverse(position) : position;
  eip->MarkICSDirty();
  eip->MarkOrganizerDirty();
  this->Update();
}

namespace {
  template<typename T>
  void set_transformed_pos(impl::EntityImpl* ent, const AtomHandleList& alist, T *positions)
  {
    bool has_tf=ent->HasTransform();
    for(AtomHandleList::const_iterator ait=alist.begin();ait!=alist.end();++ait) {
      if(ait->IsValid()) {
        geom::Vec3& tpos=ait->Impl()->TransformedPos();
        tpos[0]=static_cast<Real>(positions[0]);
        tpos[1]=static_cast<Real>(positions[1]);
        tpos[2]=static_cast<Real>(positions[2]);
        ait->Impl()->OriginalPos()=has_tf ? ent->GetTransform().ApplyInverse(tpos) : tpos;
      }
      positions+=3;
    }
    ent->MarkICSDirty();
    ent->MarkOrganizerDirty();
  }
} // anon ns

void XCSEditor::SetAtomTransformedPos(const AtomHandleList& alist, float *positions)
{
  set_transformed_pos<float>(ent_.Impl().get(),alist,positions);
  this->Update();
}

void XCSEditor::SetAtomTransformedPos(const AtomHandleList& alist, double *positions)
{
  set_transformed_pos<double>(ent_.Impl().get(),alist,positions);
  this->Update();
}

void XCSEditor::SetAtomOriginalPos(const AtomHandle& atom,
                                   const geom::Vec3& position)
{
  CheckHandleValidity(atom);
  impl::EntityImplPtr eip=ent_.Impl();
  atom.Impl()->OriginalPos()=position;
  atom.Impl()->TransformedPos() = eip->HasTransform() ? eip->GetTransform().Apply(position) : position;
  eip->MarkICSDirty();
  eip->MarkOrganizerDirty();
  this->Update();
}

namespace {
  template<typename T>
  void set_original_pos(impl::EntityImpl* ent, const AtomHandleList& alist, T *positions)
  {
    bool has_tf=ent->HasTransform();
    for(AtomHandleList::const_iterator ait=alist.begin();ait!=alist.end();++ait) {
      if(ait->IsValid()) {
        geom::Vec3& opos=ait->Impl()->OriginalPos();
        opos[0]=static_cast<Real>(positions[0]);
        opos[1]=static_cast<Real>(positions[1]);
        opos[2]=static_cast<Real>(positions[2]);
        ait->Impl()->TransformedPos()= has_tf ? ent->GetTransform().Apply(opos) : opos;
      }
      positions+=3;
    }
    ent->MarkICSDirty();
    ent->MarkOrganizerDirty();
  }
} // anon ns

void XCSEditor::SetAtomOriginalPos(const AtomHandleList& alist, float *positions)
{
  set_original_pos<float>(ent_.Impl().get(),alist,positions);
  this->Update();
}

void XCSEditor::SetAtomOriginalPos(const AtomHandleList& alist, double *positions)
{
  set_original_pos<double>(ent_.Impl().get(),alist,positions);
  this->Update();
}

void XCSEditor::SetAtomPos(const AtomHandle& atom, const geom::Vec3& position)
{
  this->SetAtomTransformedPos(atom,position);
}

void XCSEditor::SetAtomPos(const AtomHandleList& alist, float *positions)
{
  this->SetAtomTransformedPos(alist,positions);
}

void XCSEditor::SetAtomPos(const AtomHandleList& alist, double *positions)
{
  this->SetAtomTransformedPos(alist,positions);
}

void XCSEditor::ApplyTransform(const geom::Mat4& transform)
{
  geom::Transform tf;
  tf.SetMatrix(transform);
  this->ApplyTransform(tf);
}

void XCSEditor::ApplyTransform(const geom::Transform& transform)
{
  ent_.Impl()->ApplyTransform(transform);
  ent_.Impl()->UpdateTransformedPos();
  ent_.Impl()->MarkICSDirty();
  ent_.Impl()->MarkOrganizerDirty();
  this->Update();
}

void XCSEditor::SetTransform(const geom::Mat4& transform)
{
  geom::Transform tf;
  tf.SetMatrix(transform);
  this->SetTransform(tf);
}

void XCSEditor::SetTransform(const geom::Transform& transform)
{
  ent_.Impl()->SetTransform(transform);
  ent_.Impl()->UpdateTransformedPos();
  ent_.Impl()->MarkICSDirty();
  ent_.Impl()->MarkOrganizerDirty();
  this->Update();
}

void XCSEditor::FixTransform()
{
  ent_.Impl()->FixTransform();
}

void XCSEditor::Update()
{
  if (GetMode()==UNBUFFERED_EDIT) {
    ent_.Impl()->UpdateICSIfNeeded();
    ent_.Impl()->UpdateOrganizerIfNeeded();
  }
}

void XCSEditor::ForceUpdate()
{
  ent_.Impl()->UpdateICSIfNeeded();
  ent_.Impl()->UpdateOrganizerIfNeeded();
}

void XCSEditor::UpdateICS()
{
  ent_.Impl()->UpdateICSIfNeeded();
}

}} // ns
