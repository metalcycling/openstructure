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
#ifndef OST_RESIDUE_IMPL_HH
#define OST_RESIDUE_IMPL_HH

#include <boost/enable_shared_from_this.hpp>

#include <ost/geom/geom.hh>

#include <ost/mol/residue_prop.hh>
#include <ost/mol/entity_visitor_fw.hh>

#include <ost/mol/impl/residue_impl_fw.hh>
#include <ost/mol/impl/atom_impl_fw.hh>
#include <ost/mol/impl/chain_impl_fw.hh>
#include <ost/mol/impl/entity_impl_fw.hh>
#include <ost/mol/sec_structure.hh>
#include <ost/mol/impl/torsion_impl_fw.hh>
#include <ost/mol/impl/atom_group.hh>
#include <ost/mol/chem_class.hh>
#include <ost/mol/chem_type.hh>

#include <ost/generic_property.hh>
#include <ost/mol/property_id.hh>

namespace ost { namespace mol {  namespace impl {


/// \internal
/// \brief Residue implementation
class ResidueImpl: public GenericPropContainerImpl,
                   public boost::enable_shared_from_this<ResidueImpl>
                   
{
public:
  ResidueImpl(const EntityImplPtr& ent, const ChainImplPtr& ch,
              const ResNum& num, const ResidueKey& key);

  AtomImplPtr InsertAtom(const String& name, const geom::Vec3& pos, 
                         const String& ele);
  /// \brief insert new residue with exactly the same parameters as atom, but 
  ///     no bonds
  AtomImplPtr InsertAtom(const AtomImplPtr& atom);
  AtomImplPtr InsertAltAtom(const String& name, const String& alt_group,
                            const geom::Vec3& pos, const String& ele,
                            Real occ, Real b_factor);
  const ResNum& GetNumber() const {return num_;}
  void SetNumber(const ResNum& num) {num_=num;}

    const ResidueKey& GetKey() const {
    return key_;
  }

  void SetKey(const ResidueKey& key) {
    key_=key;
  }

  ChainImplPtr GetChain() const;

  AtomImplPtr GetCentralAtom() const;
  /*!
    explicitely set central atom

    if this is set, it will override the heuristic encoded
    in GetCentralAtom; pass an invalid ptr to deactivate again
  */
  void SetCentralAtom(const AtomImplPtr& a);

  geom::Vec3 GetCentralNormal() const;

  char GetOneLetterCode() const;

  void SetOneLetterCode(const char olc);

  ///\brief  Find atom by name.
  ///
  /// \todo   make the return type a AtomImplList, since there may be several
  ///         atoms with the given name.
  /// \return The first atom matchin the name, or an invalid AtomImplPtr if
  ///         there is no such atom.
  AtomImplPtr FindAtom(const String& aname) const;

  void Apply(EntityVisitor& v);

  /// \brief Get number of atoms of this residue
  int GetAtomCount() const;

  /// \brief Get number of bonds of this residue
  int GetBondCount() const;

  EntityImplPtr GetEntity() const;

  /// \brief Get phi torsion
  ///
  /// If no phi torsion has been assigned to this residue an invalid
  /// TorsionImplP will be returned.
  TorsionImplP GetPhiTorsion() const;

  /// \brief Retrieve omega torsion between this and the previous residue
  ///
  /// If no omega torsion has been assigned to this residue an invalid
  /// TorsionImplP will be returned.
  TorsionImplP GetOmegaTorsion() const;
  
  /// \brief Get psi torsion
  ///
  /// If no psi torsion has been assigned to this residue an invalid
  /// TorsionImplP will be returned.
  TorsionImplP GetPsiTorsion() const;

  ///\brief   Get secondary structure type. By default, the residue is in COIL
  ///         conformation.
  SecStructure GetSecStructure() const;

  /// \brief   Set secondary structure type
  void SetSecStructure(SecStructure ss);

  const AtomImplList& GetAtomList() const;

  AtomImplList& GetAtomList() {
    return atom_list_;
  }

  Real GetMass() const;
  geom::Vec3 GetCenterOfMass() const;
  geom::Vec3 GetCenterOfAtoms() const;
  
  geom::AlignedCuboid GetBounds() const;
  

  void DeleteAtom(const AtomImplPtr& atom);
  void DeleteAtoms(const String& atom_name);
  void DeleteAllAtoms();

  void AddTorsion(const TorsionImplP& torsion);

  const TorsionImplList& GetTorsionList() const {
    return torsion_list_;
  }
  /// \brief  Test if residue has alternative atoms with a certain name
  bool HasAltAtomGroup(const String& group) const;

  std::vector<String> GetAltAtomGroupNames() const;

  std::vector<String> GetAltAtomGroupNames(const AtomImplPtr& atom) const;
  bool HasAltAtoms() const;

  /// \brief retrieve average temperature factor of residue.
  Real GetAverageBFactor() const;
  
  bool SwitchAtomPos(const String& group);

  void AddAltAtomPos(const String& group,
                     const AtomImplPtr& atom,
                     const geom::Vec3& position,
                     Real occ, Real b_factor);
  geom::Vec3 GetAltAtomPos(const AtomImplPtr& atom, const String& group) const;
  Real GetAltAtomOcc(const AtomImplPtr& atom, const String& group) const;
  Real GetAltAtomBFactor(const AtomImplPtr& atom, const String& group) const;
  

  const String& GetCurrentAltGroupName() const {
    return curr_group_;
  }
  int GetIndex() const;

  /// \brief    Get fully qualified name of the residue, including chain name,
  ///           residue name and number
  /// \sa       AtomHandle::GetQualifiedName()
  String GetQualifiedName() const;

  /// \brief   Get name of residue. At the moment an alias for GetKey
  const String& GetName() const {
    return key_;
  }

  /// \brief   Get all torsions this residue is involved in.
  TorsionImplList& GetTorsionList()  {
    return torsion_list_;
  }

  void SetChemClass(ChemClass cc) {
    chem_class_=cc;
  }
  ChemClass GetChemClass() const {
    return chem_class_;
  }
  ChemType GetChemType() const {
    return chem_type_;
  }
  void SetChemType(ChemType ct) {
    chem_type_=ct;
  }

  TorsionImplP FindTorsion(const String& torsion_name) const;
  
  
  String GetStringProperty(Prop::ID prop_id) const;
                                                  
  Real GetFloatProperty(Prop::ID prop_id) const;
  
  int GetIntProperty(Prop::ID prop_id) const;
  
  void SetProtein(bool protein) { protein_=protein; }
  
  bool IsProtein() const { return protein_; }
  
  bool IsLigand() const { return ligand_; }
  void SetIsLigand(bool flag) { ligand_=flag; }

  void UpdateTransformedPos();

private:
  void AddAltAtom(const String& group, const AtomImplPtr& atom,
                  const geom::Vec3& position, Real occ, Real b_factor);
  void RemoveAltPositionsForAtom(const AtomImplPtr& atom);

  String                     curr_group_;
  AtomEntryGroups            alt_groups_;
  EntityImplW                ent_;
  ChainImplW                 chain_;
  ResNum                     num_;
  ResidueKey                 key_;
  AtomImplList               atom_list_;
  TorsionImplList            torsion_list_;
  SecStructure               sec_structure_;
  ChemClass                  chem_class_;
  ChemType                   chem_type_;
  char                       olc_;
  // whether the residue is part of the protein.
  // TODO: this should be fixed to be a enum'ed type aka
  // RESIDUE_TYPE type_;
  // where enum is one of protein, ligand, dna, lipid, etc
  bool                       protein_;
  bool                       ligand_;
  AtomImplPtr                central_atom_;
 
};

}}} // ns

#endif

