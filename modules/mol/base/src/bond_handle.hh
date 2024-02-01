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
#ifndef OST_BOND_HANDLE_HE
#define OST_BOND_HANDLE_HE

#include <ost/mol/impl/connector_impl_fw.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/generic_property.hh>

#include <ost/mol/entity_visitor_fw.hh>

namespace ost { namespace mol {

/// \brief chemical bond
///
/// Represents a chemical bond between two atoms (#GetFirst(), #GetSecond()).
/// New bonds are created with EditorBase::Connect().
class DLLEXPORT_OST_MOL BondHandle: 
    public GenericPropContainer<BondHandle> {
public:
  friend class ConstGenericPropContainer<BondHandle>;  
  /// necessary dummy ctor, creates invalid handle
  BondHandle();
  /// ctor for internal use, in public interface for convenience purposes
  BondHandle(const impl::ConnectorImplP& im);

  /// \name Handle validity
  //@{
  /// \brief check validity of handle
  /// 
  /// check, whether the bond handle points to a valid bond.
  /// \note It is an error to use any method other than #IsValid, #Impl and 
  ///       #operator bool() when the handle is invalid. An InvalidHandle
  ///       exception will be thrown.
  operator bool() const;

  /// \brief check validity of handle
  /// \sa #operator bool()
  bool IsValid() const;
  //@}
  
  //! entry point for entity visitor
  void Apply(EntityVisitor& v);
  void Apply(EntityViewVisitor& v);

  /// \name Atoms involved in bond
  //@{
  /// These methods return the atoms involved in the bond. No assumptions about
  /// the order should be made. With the internal coordinate system enabled,
  /// first and second may even be swapped when rebuilding the internal 
  /// connectivity tree.
  ///
  /// \brief Get first atom
  /// In python also available as the property first
  AtomHandle GetFirst() const;

  /// \brief get second atom
  /// In python also available as the property second
  AtomHandle GetSecond() const;

  /// \brief get other atom
  /// Returns one of the two atoms that does not match the given one.
  AtomHandle GetOther(const AtomHandle& a) const;
  //@}

  /// \brief return midpoint between the two atoms (transformed coordinates)
  geom::Vec3 GetPos() const;

  /// \brief return midpoint between the two atoms (original coordinates)
  geom::Vec3 GetOriginalPos() const;

  /// return the length of the bond
  Real GetLength() const;
  
  /// \brief get the bond order
  ///        ((1) single, (2) double, (3) triple, (4) aromatic)
  unsigned char GetBondOrder() const;
  
  /// \brief set the bond order
  ///        ((1) single, (2) double, (3) triple, (4) aromatic)
  void SetBondOrder(unsigned char bo);

  /// unique bond id
  unsigned long GetHashCode() const;
  
  bool operator==(const BondHandle& rhs) const;
  bool operator!=(const BondHandle& rhs) const;
  
  const impl::ConnectorImplP& Impl() const {
    return impl_;
  }
  
  impl::ConnectorImplP& Impl() {
    return impl_;
  }
  
protected:
  
  GenericPropContainerImpl* GpImpl();
  
  const GenericPropContainerImpl* GpImpl() const;  
  void CheckValidity() const;
private:
  impl::ConnectorImplP impl_;
};

/// \brief check whether a bond exists, that connects the two atom handles.
///
/// This function checks, if a bond between \c a and \c b exists and 
/// returns true if this is the case. In any other case, including when one
/// of \c a and \c b is an invalid AtomHandle, the method returns false.
/// The function does \b not check whether there may be a connection involving
/// more than one bond connecting \c a and \c b
///
/// \relates BondHandle
bool DLLEXPORT_OST_MOL BondExists(const AtomHandle& a, const AtomHandle& b);

typedef std::vector<BondHandle>  BondHandeList;

DLLEXPORT_OST_MOL std::ostream& operator<<(std::ostream& os,
                                            const BondHandle& bond);


}} // ns

#endif
