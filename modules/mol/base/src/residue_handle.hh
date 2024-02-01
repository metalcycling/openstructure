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
#ifndef OST_RESIDUE_HANDLE_HU
#define OST_RESIDUE_HANDLE_HU

#include <ost/geom/vec3.hh>
#include <ost/geom/aligned_cuboid.hh>

#include <ost/mol/query.hh>
#include <ost/mol/module_config.hh>
#include <ost/mol/residue_base.hh>
#include <ost/mol/entity_visitor_fw.hh>

namespace ost { namespace mol {

/// \brief  Residue handle
/// 
/// In a polymer a residue is a single unit such as a nucleotide or an amino
/// acid. In case of ligands, a residue ofter represents a completye molecule, 
/// such as a HEM group.
///
/// Every residue belongs to a \ref ChainHandle "chain" and has a defined 
/// position (#GetIndex) in that chain. In a peptide chain, the previous (towards 
/// N-terminus) and next (towards C-terminus) can be obtained with GetPrev() and
/// GetNext() respectively.
///
/// A residue consist of one or more \ref AtomHandle "atoms". The list of atoms
/// can be retrieved with #GetAtomList(), wheras #FindAtom() will locate a
/// particular atom in the residue.
/// 
/// To add atoms to the residue, you have to request an editor and call
/// EditorBase::InsertAtom() or EditorBase::InsertAltAtom()
/// 
/// \code
/// XCSEditor edi=residue.GetEntity().EditXCS();
/// edi.InsertAtom(residue, "CA", geom.Vec3(1.0, 0.0, 0.0));
/// \endcode
class DLLEXPORT_OST_MOL ResidueHandle : public ResidueBase {
public:
  ResidueHandle();

  //@}
  ResidueHandle(const impl::ResidueImplPtr& impl)
    : ResidueBase(impl) {}
  EntityHandle GetEntity() const;

  /// \brief Get residue's mass
  double GetMass() const;
  
  /// \brief Get residue's center of mass (mass weighted)
  geom::Vec3 GetCenterOfMass() const;
  
  /// \brief Get residue's center of atoms (not mass weighted)
  ///
  /// Returns the center of all the atoms in this residue. This is
  /// similar to GetCenterOfMass(), but the atoms are not mass weighted
  geom::Vec3 GetCenterOfAtoms() const;
  
  /// \brief get axis-aligned bounding box of residue
  geom::AlignedCuboid GetBounds() const;
  
  /// \brief The chain this residue is attached to
  ChainHandle GetChain() const;

  /// \brief  Get previous residue in chain
  /// \return Invalid handle, when there is no previous residue.
  ResidueHandle GetPrev() const;

  /// \brief  Get next residue in chain
  /// \return Invalid handle, when there is no next residue.
  ResidueHandle GetNext() const;

  /// \brief  Find atom by name.
  /// \param aname
  ///         is the atom name
  AtomHandle FindAtom(const String& aname) const;

  /// \brief  Apply entity visitor
  /// \param visitor
  ///          is the entity visitor to apply
  /// \sa     EntityVisitor
  void Apply(EntityVisitor& visitor);

  /// \brief  Get number of atoms of this residue
  ///
  ///  \see ResidueHandle::GetAtomCount
  AtomHandleList GetAtomList() const;

  /// \brief   Get number of atoms of this residue
  int GetAtomCount() const;

  /// \brief   Get number of bonds of this residue. This method might NOT work as expected (see warning below).
  /// \warning For a bond between two different residues, it is not clearly
  ///          defined to which residue it belongs. Therefore, this method
  ///          works only correctly if the residue has NO bonds to any other
  ///          residues!
  int GetBondCount() const;

  /// \brief  Find torsion by name
  /// \note   For the PHI, ISP and PSI torsions of the backbone use the designated
  ///         GetPhiTorsion and GetPsiTorsion methods, respectively.
  /// \param torsion_name
  ///         is the name of the torsion you are looking for. Torsions are named
  ///         with uppercase letters.
  /// \return The torsion, or an invalid torsion handle if the residues does not
  ///         have a torsion with the given name.
  /// \sa     ResidueBase::GetPhiTorsion, ResidueBase::GetPsiTorsion
  TorsionHandle FindTorsion(const String& torsion_name) const;

  /// \brief   Get list of all torsion angles of this residue
  TorsionHandleList GetTorsionList() const;

  /// \brief Get position of residue in chain.
  int GetIndex() const;

  /// \name     Alternative atom positions
  ///
  /// Alternative atom positions have traditionally been used to indicate that
  /// the position of an atom is not limited to one point, but that the atom may
  /// occur at more than one position. The following functions let you work with
  /// alternative atom locations on the level of residues.
  //@{
  /// \brief    Test if residue has group of alternative atoms with a certain
  ///           name.
  bool HasAltAtomGroup(const String& group) const;

  /// \brief    Name of active alternative atom group
  ///
  /// \return   The alternative atom group name or an empty String, if the atom
  ///           does not have any alternative atoms.
  const String& GetCurrentAltGroupName() const;
  
  /// \brief    Get names of alternative atom groups
  std::vector<String> GetAltAtomGroupNames() const;

  /// \brief    Test if the residue has any alternative atom positions
  bool HasAltAtoms() const;
  
  /// \brief   Get alternative atom position for given atom
  /// 
  /// If the atom is not in \c group, an IntegrityError is thrown.
  geom::Vec3 GetAltAtomPos(const AtomHandle& atom, const String& group) const;
  
  /// \brief  switch atoms that have alternative positions to the position  in
  ///         the group. Atoms that do not have an alternative position are left
  ///         untouched.
  /// \param group
  ///         is the group name to switch to.
  /// \return True, if such a group exists, false if not.
  bool SwitchAtomPos(const String& group);
  //@}

  /// \brief get this handle
  /// 
  /// This method exists as a convenience duck-typing in Python.
  ResidueHandle GetHandle() const;
  unsigned long GetHashCode() const;

  /// \brief return view based on a query object
  /// \sa Query
  EntityView Select(const Query& q, QueryFlags flags=0) const;

  /// \brief return view based on query String.
  /// \sa Query
  EntityView Select(const String& query_string, QueryFlags flags=0) const;
    
  bool operator==(const ResidueHandle& ref) const;
  bool operator!=(const ResidueHandle& ref) const;
  
  

};

/// \brief  Test if two residues are in next in sequence
///
/// Tests whether two residues are next in sequence and if they are connected by
/// a bond. At the moment, this is only implemented for peptides
bool DLLEXPORT_OST_MOL InSequence(const ResidueHandle& residue_one,
                                  const ResidueHandle& residue_two);
                                       

}} //ns

#endif // OST_RESIDUE_HANDLE_HH

