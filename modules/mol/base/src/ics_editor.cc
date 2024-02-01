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
#include "ics_editor.hh"
#include "mol.hh"
#include "impl/torsion_impl.hh"
#include <ost/invalid_handle.hh>
#include "impl/torsion_impl.hh"
#include "impl/entity_impl.hh"
#include "impl/connector_impl.hh"
#include "impl/atom_impl.hh"
/*
  Author: Marco Biasini
 */

namespace ost { namespace mol {

ICSEditor::ICSEditor(const EntityHandle& ent, EditMode mode):
  EditorBase(ent, mode)
{
  if (ent.IsValid()) {
    ent.Impl()->IncICSEditorCount();
  }  
}

ICSEditor::~ICSEditor()
{
  if (ent_.IsValid()) {
    ent_.Impl()->DecICSEditorCount();
  }
}

void ICSEditor::UpdateXCS()
{
  ent_.Impl()->UpdateXCSIfNeeded();
  ent_.Impl()->UpdateOrganizerIfNeeded();  
}

void ICSEditor::Update()
{
  if (mode_==UNBUFFERED_EDIT) {
    this->UpdateXCS();
  }
}

ICSEditor::ICSEditor(const ICSEditor& rhs):
  EditorBase(rhs.ent_, rhs.mode_)
{
  if (ent_.IsValid()) {
    ent_.Impl()->IncICSEditorCount();
  }  
}

ICSEditor& ICSEditor::operator=(const ICSEditor& rhs)
{
  ent_=rhs.ent_;
  if (ent_.IsValid()) {
    ent_.Impl()->IncICSEditorCount();
  }
  return *this;
}

void ICSEditor::SetBondLength(const BondHandle& bond, Real length)
{
    CheckHandleValidity(bond);
    ent_.Impl()->MarkXCSDirty();
    bond.Impl()->SetLength(length);
    this->Update();
}

bool ICSEditor::SetAngle(const AtomHandle& atom_a, const AtomHandle& atom_b, 
                         const AtomHandle& atom_c, Real angle)
{

    CheckHandleValidity(atom_a);
    CheckHandleValidity(atom_b);
    CheckHandleValidity(atom_c);        
    ent_.Impl()->MarkXCSDirty();    
    bool return_value = ent_.Impl()->SetAngle(atom_a.Impl(), atom_b.Impl(), atom_c.Impl(),
                                 angle);
    this->Update();
    return return_value;
}

void ICSEditor::SetTorsionAngle(TorsionHandle torsion, Real angle,
                                bool update_others)
{
  CheckHandleValidity(torsion);  
  torsion.Impl()->SetAngle(angle, update_others);
  ent_.Impl()->MarkXCSDirty();
  this->Update(); 
}

void ICSEditor::RotateTorsionAngle(TorsionHandle torsion, Real angle,
                                   bool update_others)
{
  CheckHandleValidity(torsion);
  torsion.Impl()->RotateAngle(angle, update_others);
  ent_.Impl()->MarkXCSDirty();  
  this->Update(); 
}


void ICSEditor::RotateTorsionAngle(const AtomHandle& atom_a,
                                   const AtomHandle& atom_b,
                                   const AtomHandle& atom_c,
                                   const AtomHandle& atom_d,
                                   Real angle,
                                   bool update_others)
{
  CheckHandleValidity(atom_a);
  CheckHandleValidity(atom_b);
  CheckHandleValidity(atom_c);
  CheckHandleValidity(atom_d);
  impl::Dihedral d(atom_a.Impl(), atom_b.Impl(), atom_c.Impl(), atom_d.Impl());
  Real cur_angle=d.GetAngleICS();
  d.SetAngleICS(cur_angle+angle, update_others);
  ent_.Impl()->MarkXCSDirty();
  this->Update();
}


void ICSEditor::SetTorsionAngle(const AtomHandle& atom_a, 
                                const AtomHandle& atom_b,
                                const AtomHandle& atom_c, 
                                const AtomHandle& atom_d,
                                Real angle,
                                bool update_others)
{
  CheckHandleValidity(atom_a);
  CheckHandleValidity(atom_b);
  CheckHandleValidity(atom_c);  
  CheckHandleValidity(atom_d);  
  impl::Dihedral d(atom_a.Impl(), atom_b.Impl(), atom_c.Impl(), atom_d.Impl());
  d.SetAngleICS(angle, update_others);
  ent_.Impl()->MarkXCSDirty();
  this->Update();   
}


}} // ns
