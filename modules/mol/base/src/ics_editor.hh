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
#ifndef OST_ICS_EDITOR_HE
#define OST_ICS_EDITOR_HE
/*
  Author: Marco Biasini
 */

#include <ost/mol/module_config.hh>
#include <ost/geom/geom.hh>

#include "handle_type_fw.hh"
#include "editor_base.hh"
namespace ost { namespace mol {

/// \brief editor for internal molecule coordinates
/// \sa \ref editors "Entity Editors", \ref roll_helix.py "Roll Helix Example"
class DLLEXPORT_OST_MOL ICSEditor : public EditorBase {
public:
  friend class EntityHandle;
    ~ICSEditor();
protected:
  ICSEditor(const EntityHandle& ent, EditMode mode);
  
  void Update();
public:  
  ICSEditor(const ICSEditor& rhs);
  
  ICSEditor& operator=(const ICSEditor& rhs);
      
  /// \brief  set the lengh of the bond, ie the distance of the second atom 
  ///          from the first
  ///
  /// Calling BondHandle::SetLength will bring the external coordinate system
  /// out of sync. It will not be updated until EntityHandle::UpdateICS() is
  /// called.
  ///
  /// \sa    EntityHandle::UpdateICS  
  void SetBondLength(const BondHandle& bond, Real length);
  
  /// \brief    Set angle of vectors formed by a12 and a13 to the given value.
  ///
  /// The angle is adjusted in such a way that the plane given by the 3 atom
  /// positions is kept invariant (the vectors are rotated around the plane
  /// normal).
  /// 
  /// Setting the angle requires the 3 atoms to be connected with bond. If they
  /// are not, the method will fail. The method is solely based on the internal
  /// coordinate system.

  bool SetAngle(const AtomHandle& atom_a, const AtomHandle& atom_b, 
                const AtomHandle& atom_c, Real angle);

  /// \brief  set absolute torsion angle
  /// 
  ///  All connectors at the third atom (A3) will be adjusted accordingly. If
  /// you only want to adjust the bond between A3 and A4, and leave the other
  /// bonds untouched, use the SetDihedralAngle() function
  /// 
  /// \sa     RotateTorsionAngle                 
  void SetTorsionAngle(TorsionHandle torsion, Real angle,
                       bool update_others=true);
  /// \brief  rotate torsion angle
  void RotateTorsionAngle(TorsionHandle torsion, Real delta,
                          bool update_others=true);

  void RotateTorsionAngle(const AtomHandle& atom_a,
                          const AtomHandle& atom_b,
                          const AtomHandle& atom_c,
                          const AtomHandle& atom_d,
                          Real angle,
                          bool update_others=true);

  void SetTorsionAngle(const AtomHandle& atom_a, const AtomHandle& atom_b,
                       const AtomHandle& atom_c, const AtomHandle& atom_d,
                       Real angle, bool update_others=true);
                       
  /// \brief immediately update external coordinate system
  void UpdateXCS();
};

}} // ns

/// \example roll_helix.py
/// 
/// Shows how to modify the internal coordinate system of a protein. The 
/// backtone torsion angles are all set to PHI=-60, PSI=-45; the protein is thus 
/// folded into a very long alpha helix.
#endif
