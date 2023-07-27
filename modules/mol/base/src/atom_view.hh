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
#ifndef OST_ATOM_VIEW_HH
#define OST_ATOM_VIEW_HH

/*
   Author: Marco Biasini
 */
#include <ost/mol/module_config.hh>
#include <ost/mol/atom_base.hh>
#include <ost/mol/entity_visitor_fw.hh>
#include <ost/mol/view_type_fw.hh>
#include <ost/mol/handle_type_fw.hh>

namespace ost { namespace mol {

/// \brief definition of AtomView
class DLLEXPORT_OST_MOL AtomView : public AtomBase {

  friend class EntityView;
public:
  
  /// \name View validity
  //@{
  /// \brief check validity of view
  /// 
  /// check, whether the atom view points to a valid atom.
  /// \note It is an error to use any method other than #IsValid, Impl and 
  ///       #operator bool() when the handle is invalid. An InvalidHandle
  ///       exception will be thrown.
  operator bool() const { return this->IsValid(); }
  /// \brief check validity of handle
  /// \sa #operator bool()
  bool IsValid() const;
  //@}
  // constructors  
  AtomView();
  AtomView(const ResidueView& residue_view,
           const AtomHandle& atom);

  /// \brief Get parent residue view
  ResidueView GetResidue() const;
  
  /// \brief Shortcut to entity
  EntityView GetEntity() const;
  
  /// \brief Get handle this view points to
  AtomHandle GetHandle() const;

  /// \brief Get bond list
  /// Only bonds included in the view will be returned. To get all bonds of the
  /// the atom, use `atom.GetHandle().GetBondList()`.
  BondHandleList GetBondList() const;
  
  /// \brief Get bond count
  /// Only bonds included in the view will be counted. To count all bonds of the
  /// the atom, use `atom.GetHandle().GetBondCount()`.
  int GetBondCount() const;
  
  /// \brief Get all atoms this atom is bonded to
  AtomViewList GetBondPartners() const;  

  /// \brief Remove all bonds this atom is involved in
  void RemoveBonds();
  
  /// \brief Apply entity visitor to this atom
  void Apply(EntityVisitor& visitor);
  void Apply(EntityViewVisitor& visitor);

  /// \brief get unique id
  /// 
  /// The unique id is the same for all AtomViews pointing to the same atom 
  /// view data. 
  unsigned long GetHashCode() const;
  bool operator==(const AtomView& rhs) const;
  bool operator!=(const AtomView& rhs) const;
protected:
  /// \brief Add bond to atom.
  /// If the bond is already present, it will not be added again.
  /// \param bond The bond to add
  /// \return true if the bond has been added, false if the bond was already 
  /// present in the list
  bool AddBond(const BondHandle& bond);
  
  void RemoveBondInternal(const BondHandle& bond);
private:
  AtomViewDataPtr  data_;
};

}} // ns

#endif // OST_ATOM_VIEW_HH

