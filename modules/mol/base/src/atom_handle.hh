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
#ifndef OST_ATOM_HANDLE_HR
#define OST_ATOM_HANDLE_HR

#include <ost/mol/module_config.hh>
#include <ost/mol/atom_base.hh>
#include <ost/mol/entity_visitor_fw.hh>
#include <ost/mol/handle_type_fw.hh>

namespace ost { namespace  mol {

/// \brief Handle to atom datatype
///
/// The bonds of an atom are accessible via #GetBondList(). For a list of atoms
/// that this atom is bonded to, use #GetBondPartners().
///
/// New bonds are formed with one of the EditorBase::Connect() methods.
/// 
/// A valid atom handle is always part of a residue.
class DLLEXPORT_OST_MOL AtomHandle : public AtomBase {
public:
  AtomHandle();
  AtomHandle(const impl::AtomImplPtr& impl);

public:  

  ResidueHandle GetResidue() const;
  
  EntityHandle GetEntity() const;  
  
  /// \name Bonds
  //@{
  /// \brief Get number of bonds this atom forms.
  ///
  /// Count of all bonds where this atom is either the first or second atom.
  /// \sa BondHandle::GetFirst, BondHandle::GetSecond
  int GetBondCount() const;
  
  /// \brief get bond list
  ///
  /// Returns all bonds where this atom is either the first or second atom.
  /// \sa BondHandle::GetFirst, BondHandle::GetSecond
  BondHandleList GetBondList() const;
  
  /// \brief Find bond to specific atom
  /// 
  /// \return The bond or an invalid bond handle if no such bond exists.
  BondHandle FindBondToAtom(const AtomHandle& handle) const;  
  
  /// \brief Get all atoms this atom is bonded to
  AtomHandleList GetBondPartners() const;  
  //@}
  
  /// \name Traversing the hierarchy
  //@{
  /// \brief Apply visitor.
  /// 
  /// This method is used by the EntityVisitor interface to traverse the
  /// entity - chain - residue -atom hierarchy.
  /// \sa EntityVisitor
  void Apply(EntityVisitor&);
  //@}

  /// \brief get this handle
  /// 
  /// Useful for duck-typing in Python and templated code.
  AtomHandle GetHandle() const;
  
  /// \brief Get unique identifier for atom
  /// 
  /// Get hash code that uniquely identifies every atom. The hash code is
  /// identical for all atom views pointing to a given atom.
  unsigned long GetHashCode() const;
  
  bool operator==(const AtomHandle& ref) const;
  bool operator!=(const AtomHandle& ref) const;

};

}} // ns

#endif // OST_ATOM_HANDLE_HH
