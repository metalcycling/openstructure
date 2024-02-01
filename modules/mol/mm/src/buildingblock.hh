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

#ifndef OST_MM_BUILDING_BLOCK_HH
#define OST_MM_BUILDING_BLOCK_HH

#include <vector>
#include <limits>
#include <set>

#include <boost/shared_ptr.hpp>

#include <ost/mol/mm/interaction.hh>
#include <ost/message.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/geom/vec3.hh>


namespace ost { namespace mol{ namespace mm{

class BuildingBlock;
typedef boost::shared_ptr<BuildingBlock> BuildingBlockPtr;

class BuildingBlock{
public:
  BuildingBlock() { }

  //copy constructor needs to be defined explicitely, since the forces
  //have to be newly created
  BuildingBlock(const BuildingBlock& block); 

  bool Match(const ost::mol::ResidueHandle& handle, bool match_connectivity, String& info) const;

  void Connect(ost::mol::ResidueHandle& handle, ost::mol::XCSEditor& ed);

  //getter functionality for all buildingblock members

  std::vector<String> GetAtoms() const { return atoms_; }

  std::vector<String> GetTypes() const { return types_; }

  std::vector<Real> GetCharges() const { return charges_; }

  std::vector<Real> GetMasses() const { return masses_; }

  String GetType(const String& name) const;

  Real GetCharge(const String& name) const;

  Real GetMass(const String& name) const;

  std::vector<InteractionPtr> GetBonds() const { return bonds_; }

  std::vector<InteractionPtr> GetAngles() const { return angles_; }

  std::vector<InteractionPtr> GetDihedrals() const { return dihedrals_; }

  std::vector<InteractionPtr> GetImpropers() const { return impropers_; }

  std::vector<InteractionPtr> GetCMaps() const { return cmaps_; }

  std::vector<InteractionPtr> GetExclusions() const { return exclusions_; }

  std::vector<InteractionPtr> GetConstraints() const { return constraints_;}

  //Add data to building block

  void AddAtom(const String& name, const String& type, Real charge, Real mass = std::numeric_limits<Real>::quiet_NaN());

  void AddBond(InteractionPtr p, bool replace_existing = false);

  void AddAngle(InteractionPtr p, bool replace_existing = false);

  void AddDihedral(InteractionPtr p, bool replace_existing = false);

  void AddImproper(InteractionPtr p, bool replace_existing = false);

  void AddExclusion(InteractionPtr p, bool replace_existing = false);

  void AddCMap(InteractionPtr p, bool replace_existing = false);

  void AddConstraint(InteractionPtr p, bool replace_existing = false);

  //modifiers

  //removes atom and all interactions associated to it
  void RemoveAtom(const String& name);

  //replaces atom in all interactions
  void ReplaceAtom(const String& name,const String& new_name,
                   const String& new_type, Real new_charge, 
                   Real new_mass = std::numeric_limits<Real>::quiet_NaN());


  //remove all interactions to previous or next residue
  void RemoveInteractionsToPrev();

  void RemoveInteractionsToNext();

  template <typename DS>
  void Serialize(DS& ds){
    int num_atoms = atoms_.size();
    ds & num_atoms;

    if(ds.IsSource()){
      atoms_ = std::vector<String>(num_atoms);
      types_ = std::vector<String>(num_atoms);
      charges_ = std::vector<Real>(num_atoms);
      masses_ = std::vector<Real>(num_atoms);
    }

    for(int i = 0; i < num_atoms; ++i){
      ds & atoms_[i];
      ds & types_[i];
      ds & charges_[i];
      ds & masses_[i];
    }

    int num_bonds = bonds_.size();
    int num_angles = angles_.size();
    int num_dihedrals = dihedrals_.size();
    int num_impropers = impropers_.size();
    int num_exclusions = exclusions_.size();
    int num_cmaps = cmaps_.size();
    int num_constraints = constraints_.size();

    ds & num_bonds;
    ds & num_angles;
    ds & num_dihedrals;
    ds & num_impropers;
    ds & num_exclusions;
    ds & num_cmaps;
    ds & num_constraints;

    for(int i = 0; i < num_bonds; ++i){
      int func_type;
      if(ds.IsSource()){
        ds & func_type;
        bonds_.push_back(InteractionPtr(new Interaction(FuncType(func_type))));
      }
      else{
        func_type = bonds_[i]->GetFuncType();
        ds & func_type;
      }
      ds & *(bonds_[i]);
    }

    for(int i = 0; i < num_angles; ++i){
      int func_type;
      if(ds.IsSource()){
        ds & func_type;
        angles_.push_back(InteractionPtr(new Interaction(FuncType(func_type))));
      }
      else{
        func_type = angles_[i]->GetFuncType();
        ds & func_type;
      }
      ds & *(angles_[i]);
    }

    for(int i = 0; i < num_dihedrals; ++i){
      int func_type;
      if(ds.IsSource()){
        ds & func_type;
        dihedrals_.push_back(InteractionPtr(new Interaction(FuncType(func_type))));
      }
      else{
        func_type = dihedrals_[i]->GetFuncType();
        ds & func_type;
      }
      ds & *(dihedrals_[i]);
    }

    for(int i = 0; i < num_impropers; ++i){
      int func_type;
      if(ds.IsSource()){
        ds & func_type;
        impropers_.push_back(InteractionPtr(new Interaction(FuncType(func_type))));
      }
      else{
        func_type = impropers_[i]->GetFuncType();
        ds & func_type;
      }
      ds & *(impropers_[i]);
    }

    for(int i = 0; i < num_exclusions; ++i){
      int func_type;
      if(ds.IsSource()){
        ds & func_type;
        exclusions_.push_back(InteractionPtr(new Interaction(FuncType(func_type))));
      }
      else{
        func_type = exclusions_[i]->GetFuncType();
        ds & func_type;
      }
      ds & *(exclusions_[i]);
    }

    for(int i = 0; i < num_cmaps; ++i){
      int func_type;
      if(ds.IsSource()){
        ds & func_type;
        cmaps_.push_back(InteractionPtr(new Interaction(FuncType(func_type))));
      }
      else{
        func_type = cmaps_[i]->GetFuncType();
        ds & func_type;
      }
      ds & *(cmaps_[i]);
    }

    for(int i = 0; i < num_constraints; ++i){
      int func_type;
      if(ds.IsSource()){
        ds & func_type;
        constraints_.push_back(InteractionPtr(new Interaction(FuncType(func_type))));
      }
      else{
        func_type = constraints_[i]->GetFuncType();
        ds & func_type;
      }
      ds & *(constraints_[i]);
    }
  }


private:

  int GetAtomIndex(const String& atom_name) const;
  void CheckInteractionToAdd(InteractionPtr p) const;

  std::vector<String> atoms_;
  std::vector<String> types_;
  std::vector<Real> charges_;
  std::vector<Real> masses_;
  std::vector<InteractionPtr> bonds_;
  std::vector<InteractionPtr> angles_;
  std::vector<InteractionPtr> dihedrals_;
  std::vector<InteractionPtr> impropers_;
  std::vector<InteractionPtr> exclusions_;
  std::vector<InteractionPtr> cmaps_;
  std::vector<InteractionPtr> constraints_;
};

}}}

#endif
