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
#ifndef OST_TORSION_HANDLE_HH
#define OST_TORSION_HANDLE_HH

#include <ost/mol/impl/torsion_impl_fw.hh>
#include <ost//mol/atom_handle.hh>

#include <ost/mol/entity_visitor_fw.hh>

namespace ost { namespace mol {

/// \brief Named torsion (dihedral) angle
/// 
/// A torsion (angle) is determined by three consecutive connections, which 
/// means by four consecutive atoms, and is defined by the two planes, one 
/// formed by A1-A2-A3, the other by A2-A3-A4
/// 
/// Torsion are named, such that they can be retrieved later by name. They are 
/// usually added to the entity when loaded from file. For ab-initio model 
/// building, the torsions may be added manually or by an  appropriate 
/// builder from the conop module such as conop::HeuristicBuilder.
/// 
/// The angle of the torsion may be adjusted with via a call to 
/// XCSEditor::SetTorsionAngle.
/// 
/// \sa \ref roll_helix.py "Roll Helix Example"
class DLLEXPORT_OST_MOL TorsionHandle {
public:
  //! necessary dummy ctor, creates invalid handle
  TorsionHandle();
  //! ctor for internal use, in public interface for convenience purposes
  TorsionHandle(const impl::TorsionImplP& im);

  /// \name Handle validity
  //@{
  /// \brief check validity of handle
  /// 
  /// check, whether the torsion handle points to a valid torsion.
  /// \note It is an error to use any method other than #IsValid, #Impl() and 
  ///       #operator bool() when the handle is invalid. An InvalidHandle
  ///       exception will be thrown.
  operator bool() const;
  /// \brief check validity of handle
  /// \sa #operator bool()
  bool IsValid() const;
  //@}
  
  //! entry point for entity visitor
  void Apply(EntityVisitor& v);

  const String& GetName() const;

  /// \name Accessing atoms
  //@{
  /// \brief  return first atom
  AtomHandle GetFirst() const;
  /// \brief  return second atom
  AtomHandle GetSecond() const;
  /// \brief  return third atom
  AtomHandle GetThird() const;
  /// \brief  return fourth atom
  AtomHandle GetFourth() const;
  //@}
    
  /// \brief get dihedral angle based on whatever coordinate system is 
  ///        available  
  Real GetAngle() const;
  /// \brief get position, i.e midpoint between the second and third atom (transformed coordinates)
  geom::Vec3 GetPos() const;
    /// \brief get position, i.e midpoint between the second and third atom (original coordinates)
  geom::Vec3 GetOriginalPos() const;

  impl::TorsionImplP& Impl() {
    return impl_;
  }
  
  const impl::TorsionImplP& Impl() const {
    return impl_;
  }
protected:
    void CheckValidity() const; 
private:
  impl::TorsionImplP impl_;
};
                      
Real DLLEXPORT_OST_MOL DihedralAngle(const AtomHandle& a1, 
                                        const AtomHandle& a2,
                                        const AtomHandle& a3, 
                                        const AtomHandle& a4);

DLLEXPORT_OST_MOL std::ostream& operator<<(std::ostream& os,
                                            const TorsionHandle& torsion);

}} // ns

#endif
