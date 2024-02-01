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
#ifndef OST_BASE_DIHEDRAL_HH
#define OST_BASE_DIHEDRAL_HH

#include <ost/mol/module_config.hh>
#include <ost/mol/impl/atom_impl_fw.hh>
#include <ost/geom/vec3.hh>

namespace ost { namespace mol { namespace impl {
  
/// \brief Low level dihedral implementation
/// 
/// \sa \ref roll_helix.py "Roll Helix Example"
class Dihedral {
public:
  /// \brief construct new dihedral as defined by the four consecutive atoms
  ///        atom1 to atom4
  Dihedral(const AtomImplPtr& atom1, 
           const AtomImplPtr& atom2, 
           const AtomImplPtr& atom3,
           const AtomImplPtr& atom4);
  /// \brief construct new dihedral as defined by the four consecutive atoms
  ///        stored in atoms
  Dihedral(const AtomImplList& atoms);
    
  /// \brief get dihedral angle base on atom position
  Real GetAngleXCS() const;
  
  /// \brief get dihedral angle based on internal coordinate system
  Real GetAngleICS() const;
  
  /// \brief set dihedral angle. Requires an intact internal coordinate system.
  /// \param angle 
  ///         new angle in radians
  /// \param update_other 
  ///         If true, connectors of atom3 are updated accordingly, otherwise 
  ///         they will be left untouched
  void SetAngleICS(Real angle, bool update_other);
  
  geom::Vec3 GetPos() const;
  geom::Vec3 GetOriginalPos() const;

  bool IsAtomInvolved(const AtomImplPtr& atom);
    
  AtomImplPtr GetFirst() const;
  
  AtomImplPtr GetSecond() const;
  
  AtomImplPtr GetThird() const;
  
  AtomImplPtr GetFourth() const;
  
  bool Matches(const AtomImplPtr& a1, const AtomImplPtr& a2, 
               const AtomImplPtr& a3, const AtomImplPtr& a4) const;  
    
private:
  AtomImplList atoms_;
};

}}}

#endif
