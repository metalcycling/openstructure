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
#ifndef OST_CHAIN_HANDLE_HS
#define OST_CHAIN_HANDLE_HS

#include <vector>

#include <ost/mol/module_config.hh>

#include <ost/geom/geom.hh>

#include "chain_base.hh"

#include "query.hh"
#include "view_type_fw.hh"
#include "handle_type_fw.hh"
#include "residue_prop.hh"
#include "entity_visitor_fw.hh"
#include "sec_structure.hh"

namespace ost { namespace mol {

/// \brief linear chain of residues
///
/// In peptide chain, the residues are ordered from N- to C-terminus. To obtain
/// the n-th residue in the chain, use #GetResidueByIndex. Alternatively, the
/// residues may be adressed by their residue number with #FindResidue. 
/// Residue numbers are unique within a chain. It is however not guaranteed that 
/// the residue numbers are in ascending order.
///
/// New residues are inserted via one of the editors obtained from 
/// EntityHandle::EditXCS() or EntityHandle::EditICS(). To
/// append a residue to the C-terminal end, use EditorBase::AppendResidue, to
/// insert after or before a specific residue, EditorBase::InsertResidueAfter
/// and ditorBase::InsertResidueABefore will do the job, respectively.
class DLLEXPORT_OST_MOL ChainHandle : public ChainBase {

public:
  ChainHandle();

  ChainHandle(const impl::ChainImplPtr& impl);
  

  /// \brief Get parent entity.
  EntityHandle GetEntity() const;

  /// \brief  Apply visitor
  /// \param visitor
  ///         is the visitor to apply
  void Apply(EntityVisitor& visitor);

  //@}
  
  /// \brief  Get residue before the given residue
  /// 
  ///  A previous residue is returned regardless whether the residues are
  /// connected by a bond.
  ResidueHandle GetPrev(const ResidueHandle& rh);
  ResidueHandle GetNext(const ResidueHandle& rh);

  /// \name Single item addressing
  //@{
  /// \brief  Find residue by number. If the chain does not contain a residue
  ///         with the given number, an invalid ResidueHandle will be returned.
  /// 
  /// \sa     ResidueBase::IsValid, #FindAtom, EntityHandle::FindResidue
  ResidueHandle FindResidue(const ResNum& num) const;

  /// \brief  Find atom of of specific residue
  /// 
  ///  Find atom with \var atom_name that belongs to the residue with residue
  /// number \var num.
  /// 
  /// \sa     ResidueHandle::IsValid, #FindResidue, EntityHandle::FindAtom
  /// \return The found atom, or an invalid atom handle if no such atom exists.
  AtomHandle FindAtom(const ResNum& num, 
                      const String& atom_name) const;
  /// \brief  Get residue by index
  /// 
  ///  Indices are start at zero
  ResidueHandle GetResidueByIndex(int index) const;
  //@}
  
  //! \brief   Get number of atoms of this chain.
  //
  // To count all atoms of the entity the chain belongs to, use 
  // \c EntityHandle::GetAtomCount(). For count of specific residues, use
  // \c ResidueHandle::GetResidueCount().
  int GetAtomCount() const;

  /// \brief   Get number of all inclusive bonds of this chain.
  ///
  /// \warning For a bond between two different chains, it is not clearly
  ///          defined to which chain it belongs. Therefore, this method will
  ///          only return the number of bonds where both bonded atoms are
  ///          within the same chain (all inclusive bonds).
  int GetBondCount() const;
  
  /// \brief   Get number of residues of this chain.
  ///
  /// \sa EntityHandle::GetResidueCount()
  int GetResidueCount() const;
  
  /// \brief Retrieve atoms of this chain as a list
  ///
  /// The atoms are ordered by atom number in ascending order. If you only
  /// want to determine the number of atoms, use the dedicated
  /// \c ChainHandle::GetAtomCount  
  
  
  AtomHandleList GetAtomList() const;
  
  /// \brief Get residues of this chain as a list
  ///
  /// The residues are ordered by residue number in ascending order. If you only
  /// want to determine the number of residues, use the dedicated
  /// \c ChainHandle::GetResidueCount
  ResidueHandleList GetResidueList() const;

  /// \brief Get entity's mass
  Real GetMass() const;
  
  /// \brief Get entity's center of mass (mass weighted)
  geom::Vec3 GetCenterOfMass() const;
  
  /// \brief Get entity's center of atoms (not mass weighted)
  ///
  /// Returns the center of all the atoms in an entity. This is
  /// similar to GetCenterOfMass(), but the atoms are not mass weighted
  geom::Vec3 GetCenterOfAtoms() const;
  
  geom::AlignedCuboid GetBounds() const;  
  /// \brief assign secondary structure to the inclusive residue range 
  ///        start, end
  void AssignSecondaryStructure(SecStructure ss, 
                                const ResNum& start, 
                                const ResNum& end);

  /// \brief return view based on a query object
  /// \sa Query
  EntityView Select(const Query& q, QueryFlags flags=0) const;

  /// \brief return view based on query String.
  /// \sa Query
  EntityView Select(const String& query_string, QueryFlags flags=0) const;

  /// \brief get this handle
  /// 
  /// Useful for duck-typing in Python and in templates
  ChainHandle GetHandle() const;

  unsigned long GetHashCode() const;
  
  /// \brief whether the residues form an ordered sequence with respect to their
  ///     reside numbers.
  bool InSequence() const;
  bool operator==(const ChainHandle& ref) const;
  bool operator!=(const ChainHandle& ref) const;

  /// \brief checks whether res breaks the in sequence property
  ///        and updates it accordingly
  void SetInSequence(const int index);
};

}} // ns

#endif // OST_CHAIN_HANDLE_HH

