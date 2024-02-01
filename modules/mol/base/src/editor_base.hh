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
#ifndef OST_EDITOR_BASE_HS
#define OST_EDITOR_BASE_HS
/*
  Author: Marco Biasini
 */
 
#include <ost/mol/module_config.hh>

#include <ost/geom/geom.hh>

#include "entity_handle.hh"
#include "residue_prop.hh"
#include "editor_type_fw.hh"

namespace ost { namespace mol {

/// \brief base class for entity editors
///
/// Entities, chains, residues and atoms are assembled and modified via 
/// editors.
/// Editors are usually used as one of the concrete subclasses of EditorBase:
/// XCSEditor and ICSEditor. The former operates on the external carthesian
/// coordinate system while the latter operates on the internal coordinate
/// system. To create a new editor, use EntityHandle::EditXCS() and
/// EntityHandle::EditICS(), respectively.
/// 
/// For an introduction to the concept of editors go \ref editors "here"
class DLLEXPORT_OST_MOL EditorBase {
public:
  /// \brief create new chain
  ///
  /// create new chain and attach it to entity.
  ///
  /// \param chain_name
  ///          The chain name. Can be an arbitrary String. However, if you
  ///          intend to use the PDB export funtionality, the \c chain_name
  ///          should be a single letter, preferably in the range A-Z.
  ///
  /// \return  The newly created chain handle
  ChainHandle InsertChain(const String& chain_name);

  /// \brief create new chain
  ///
  /// create new chain with properties of a provided chain handle and attach it to entity.
  ///
  /// \param chain_name
  ///          The chain name. Can be an arbitrary String. However, if you
  ///          intend to use the PDB export funtionality, the \c chain_name
  ///          should be a single letter, preferably in the range A-Z.
  /// \param chain
  ///          The newly created chain will take over all generic
  ///          attached to this handle.
  /// \param deep
  ///          If set to true, all residues and atoms of chain will be
  ///          completely copied into the created chain. No bonds and angles
  ///          are added.
  ///
  /// \return  The newly created chain handle
  ChainHandle InsertChain(const String& chain_name, ChainHandle chain, bool deep=false);

  /// \name Inserting, removing and modifying order of residues
  ///@{
  /// \brief Append residue to the end of the chain
  /// 
  /// Append residue with residue key \var k to the end of the chain. The 
  /// residue number will be set to the residue number of the last added residue
  /// plus one but with the alternative position identifier left unchanged.
  /// The following code will append a new residue with key GLY:
  /// \code
  /// ResidueHandle gly=chain.AppendResidue("GLY");
  /// \endcode
  /// \note     The created residue is empty and does not contain any atoms.
  /// \return   The newly created residue
  /// \sa       #InsertResidueBefore, EntityHandle::InsertChain, CreateEntity,
  ///           ResidueHandle::InsertAtom.
  ResidueHandle AppendResidue(ChainHandle chain, const ResidueKey& k);
  
  ResidueHandle AppendResidue(ChainHandle chain, const ResidueKey& k, 
                              const ResNum& num); 

  /// \brief Append residue to the ent of the chain
  ///
  /// Append residue with all properties of provided residue into chain.
  ///
  /// \param residue
  ///          All properties of this residue will be copied into the newly created
  ///          residue. The newly created residue will not contain any atoms, except
  ///          you set \var deep to true.
  ///
  /// \param deep
  ///          if set to true, all atoms from the source residue will be copied into
  ///          the newly created residue. No bonds and angles are added.
  ///
  /// \return  The newly created residue handle

  ResidueHandle AppendResidue(ChainHandle chain, ResidueHandle residue, bool deep=false);

  /// \brief  Insert residue into chain
  /// 
  /// Insert residue with residue number \var num and key \var k into the
  /// chain.
  /// \code
  /// EntityHandle e=CreateEntity();
  /// ChainHandle chain=e.InsertChain("A");
  /// ResidueHandle first=chain.AppendResidue("MET");
  /// ResNum number=first.GetNumber();  
  /// //the next two lines have the same effect  
  /// chain.AppendResidue("GLY");
  /// chain.InsertResidueBefore("GLY, "++number)
  /// \endcode
  /// \sa     #AppendResidue, EntityHandle::InsertChain, CreateEntity,
  ///         ResidueHandle::InsertAtom.
  ResidueHandle InsertResidueBefore(ChainHandle chain, int index, 
                                    const ResNum& num,
                                    const ResidueKey& k);

  /// \brief insert residue into chain
  /// 
  /// This method is identical to InsertResidueBefore() but inserts the
  /// residue after the specified position
  /// \sa #InsertResidueBefore
  ResidueHandle InsertResidueAfter(ChainHandle chain, int index,
                                   const ResNum& num,
                                   const ResidueKey& k);
  
  /// \name Adding/removing atoms
  //@{
  /// \brief  Insert new atom
  ///
  /// Insert new atom and add it to residue. For atoms with alternative atom
  /// locations use ResidueHandle::InsertAltAtom().
  /// \param residue is the parent residue
  /// \param name is the atom name. While free to choose a name, it is adviced 
  ///     to properly name the residues according to IUPAC rules as several
  ///     algorithms as well as most builders in the conop module rely on proper 
  ///     naming.
  /// \param pos is the position of the atom in global coordinates
  /// \param ele
  /// \param occupancy
  /// \param b_factor
  /// \param is_hetatm
  AtomHandle InsertAtom(ResidueHandle residue, const String& name,
                        const geom::Vec3& pos, const String& ele="",
                        Real occupancy=1.0, Real b_factor=0.0, 
                        bool is_hetatm=false);

  /// \brief Insert new atom
  ///
  /// Inserts new atom with all properties from the provided atom handle.
  ///
  /// \param residue is the parent residue
  /// \param atom from which all informations will be copied over to the
  ///        newly created atom
  ///
  /// \return the newly created AtomHandle

  AtomHandle InsertAtom(ResidueHandle residue, AtomHandle atom);

  /// \brief Insert new atom with alternative position indicator
  /// \sa EditorBase::AddAltAtomPos(), ResidueHandle
  AtomHandle InsertAltAtom(ResidueHandle residue, const String& name, 
                           const String& alt_group, const geom::Vec3& pos,
                           const String& ele="", Real occ=1.0,
                           Real b_factor=0.0);

  /// \brief Insert new atom with alternative position indicator
  ///
  /// All informations will be copied over from atom, except bonds

  AtomHandle InsertAltAtom(ResidueHandle residue, AtomHandle atom,
                           const String& alt_group);


  /// \brief  Add alternative atom position
  /// \param group is the name of the alternative atom position group. If no 
  ///     group of that name exists, it will be created.
  /// \param atom is the atom handle the additional location belongs to. It is
  ///     required that the atom has been inserted via InsertAltAtom, If the
  ///     atom is a conventional atom without alternative location, a
  ///     Error will be thrown.
  /// \param position
  ///         is the alternative position
  /// \param occ
  /// \param b_factor
  /// \sa EditorBase::InsertAltAtom(), ResidueHandle
  void AddAltAtomPos(const String& group, const AtomHandle& atom, 
                     const geom::Vec3& position, Real occ=1.0,
                     Real b_factor=0.0);
  //\}
  
  /// \brief connect two atoms with bond
  ///
  /// If there is already a bond between the two atoms, no new bond is created.
  /// Instead, the already existing bond is returned.
  BondHandle Connect(const AtomHandle& first, const AtomHandle& second);
  BondHandle Connect(const AtomHandle& first, const AtomHandle& second,
                     unsigned char bond_order);
  BondHandle Connect(const AtomHandle& first, const AtomHandle& second,
                     Real len, Real theta, Real phi);
  BondHandle Connect(const AtomHandle& first, const AtomHandle& second,
                     Real len, Real theta, Real phi,
                     unsigned char bond_order);

  void RenameResidue(ResidueHandle res, const String& new_name);
  
  void SetResidueNumber(ResidueHandle res, const ResNum& num);
   
  void RenameChain(ChainHandle chain, const String& new_name);

  /// \brief Assign type of chain according to ChainType.
  ///
  /// \param chain chain to assign to
  /// \param type type of the chain
  void SetChainType(ChainHandle chain, const ChainType type);

  /// \brief Assign a description to a chain.
  ///
  /// \param chain chain to assign to
  /// \param desc description
  void SetChainDescription(ChainHandle chain, const String desc);

  /// \brief   Delete all atoms of residue
  ///
  /// All associated torsions and bonds will also be removed
  void DeleteAllAtoms(ResidueHandle res);
  
  /// \brief  Delete specific atom from this residue
  ///
  ///  All associated torsions and bonds will also be removed as well.
  ///
  /// \param atom
  ///         is the atom to remove. If no such atom exists, this method will
  ///         have no effect
  void DeleteAtom(const AtomHandle& atom);
  /// \ brief Delete a set of atoms
  ///
  /// All associated torsions and bonds will be removed as well
  void DeleteAtoms(const AtomHandleList& atoms);

  /// \ brief Delete bond
  ///
  /// \param bond
  ///         Is the bond to remove. If no such bond exists, this method will
  ///         have no effect
  void DeleteBond(const BondHandle& bond);

  /// \ brief Delete a set of bond
  ///
  /// \param bonds
  ///         bonds to remove. If no such bonds exist, this method will
  ///         have no effect
  void DeleteBonds(const BondHandleList& bonds);
  
  /// \brief Add named torsion to entity
  TorsionHandle AddTorsion(const String& name, const AtomHandle& a1,
                           const AtomHandle& a2, const AtomHandle& a3,
                           const AtomHandle& a4);  
  /// \brief delete chain
  ///
  /// all associated atoms, residues, torsions and bonds will also be deleted.
  void DeleteChain(const ChainHandle& chain);

  /// \brief Delete residue
  ///
  /// Atoms, bonds and torsions involving this residue are removed as well
  ///
  /// \sa ChainHandle::DeleteResidue, DeleteAtom
  void DeleteResidue(const ResidueHandle& residue_handle);

  /// \brief reorder residues of given chain based on their residue number
  void ReorderResidues(const ChainHandle& chain);

  /// \brief reorder residues of all chains based on their residue number
  void ReorderAllResidues();

  /// \brief renumber residues of all chains
  ///
  /// \param start
  ///           Residues of every chain will be renumbered, whereas the first
  ///           residue gets the residue number start.
  ///
  /// \param keep_spacing
  ///           If set to false, residues will continously be renumbered ongoing from start.
  ///           Otherwise the spacings between the residues are kept.
  void RenumberAllResidues(int start, bool keep_spacing);


  void RenumberChain(ChainHandle chain, const ResNumList& new_numbers);

  /// \brief renumber residues of one chain
  ///
  /// \param chain
  ///           All residues of this chain will be renumbered according to the
  ///           parameters start and keep_spacing
  ///
  /// \param start
  ///           Residues of given chain will be renumbered, whereas the first
  ///           residue gets the residue number start.
  ///
  /// \param keep_spacing
  ///           If set to false, residues will continously be renumbered ongoing from start.
  ///           Otherwise the spacings between the residues are kept.
  void RenumberChain(const ChainHandle& chain, int start, bool keep_spacing);
    
  /// \brief Get edit mode of editor
  EditMode GetMode() const {return mode_;}
  
  /// \ brief return entity this editor works on
  EntityHandle GetEntity() const {return ent_;}

  /// \brief change the name of the atom to the new name  
  void RenameAtom(AtomHandle atom, const String& new_name);

  /// \brief Removes all residues and chains in the attached entity that don't 
  ///        contain any atoms
  void Prune(); 

protected:
  EditorBase(const EntityHandle& ent, EditMode mode);
  void UpdateTrace();
  EntityHandle     ent_;
  EditMode         mode_;
};
 
}} //ns

#endif
