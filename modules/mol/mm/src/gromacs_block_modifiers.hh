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

#ifndef OST_MM_GROMCACS_BLOCK_MODIFIERS_HH
#define OST_MM_GROMCACS_BLOCK_MODIFIERS_HH


#include <boost/shared_ptr.hpp>

#include <ost/message.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/geom/vec3.hh>
#include <ost/mol/mm/block_modifiers.hh>
#include <ost/mol/mm/buildingblock.hh>


namespace ost{ namespace mol{ namespace mm{

class GromacsHydrogenConstructor;
class GromacsTerminiConstructor;
class GromacsBlockModifier;

typedef boost::shared_ptr<GromacsHydrogenConstructor> GromacsHydrogenConstructorPtr;
typedef boost::shared_ptr<GromacsTerminiConstructor> GromacsTerminiConstructorPtr;
typedef boost::shared_ptr<GromacsBlockModifier> GromacsBlockModifierPtr; 

class GromacsPositionRuleEvaluator{

public:
  static std::vector<geom::Vec3> EvaluatePosRule(int rule, int number,
                                                 const std::vector<geom::Vec3>& anchors);
};

class GromacsHydrogenConstructor : public HydrogenConstructor{

public:

  GromacsHydrogenConstructor() { }

  virtual void ApplyOnBuildingBlock(BuildingBlockPtr);

  virtual void ApplyOnResidue(ost::mol::ResidueHandle& res, ost::mol::XCSEditor& ed);

  void AddHydrogenRule(uint number, int method, 
                       const std::vector<String>& hydrogen_names, 
                       const std::vector<String>& anchors);

  // GET STUFF (for data extraction)
  uint GetNumHydrogenRules() const { return add_number_.size(); }
  // all getters for idx < GetNumHydrogenRules()
  int GetHydrogenRuleNumber(uint idx) const { return add_number_[idx]; }
  int GetHydrogenRuleMethod(uint idx) const { return methods_[idx]; }
  const std::vector<String>& 
  GetHydrogenRuleNames(uint idx) const { return hydrogen_names_[idx]; }
  const std::vector<String>& 
  GetHydrogenRuleAnchors(uint idx) const { return anchor_atom_names_[idx]; }

  virtual void OnSave(ost::io::BinaryDataSink& ds) { ds << *this; }

  virtual BlockModifierType GetBlockModifierType() { return GromacsBlockModifiers; }

  template <typename DS>
  void Serialize(DS& ds){
    int num_hydrogen_add_rules = add_number_.size();
    ds & num_hydrogen_add_rules;

    if(ds.IsSource()){
      add_number_ = std::vector<int>(num_hydrogen_add_rules);
      methods_ = std::vector<int>(num_hydrogen_add_rules);
      hydrogen_names_ = std::vector<std::vector<String> >(num_hydrogen_add_rules);
      anchor_atom_names_ = std::vector<std::vector<String> >(num_hydrogen_add_rules);
    }

    for(int i = 0; i < num_hydrogen_add_rules; ++i){
      ds & add_number_[i];
      ds & methods_[i];
      int num_hydrogen_names = hydrogen_names_[i].size();
      int num_anchor_names = anchor_atom_names_[i].size();
      ds & num_hydrogen_names;
      ds & num_anchor_names;
      if(ds.IsSource()){
        hydrogen_names_[i] = std::vector<String>(num_hydrogen_names);
        anchor_atom_names_[i] = std::vector<String>(num_anchor_names);
      }

      for(int j = 0; j < num_hydrogen_names; ++j){
        ds & hydrogen_names_[i][j];
      }
      for(int j = 0; j < num_anchor_names; ++j){
        ds & anchor_atom_names_[i][j];
      }
    }
  }

private:
  std::vector<int> add_number_;
  std::vector<int> methods_;
  std::vector<std::vector<String> > hydrogen_names_;
  std::vector<std::vector<String> > anchor_atom_names_;
};


class GromacsBlockModifier : public BlockModifier{

public:

  GromacsBlockModifier() { }

  virtual void ApplyOnBuildingBlock(BuildingBlockPtr p);

  virtual void ApplyOnResidue(ost::mol::ResidueHandle& res, ost::mol::XCSEditor& ed);

  void AddReplaceRule(const String& name, const String& new_name,
                                 const String& new_type, Real new_charge);

  void AddAddRule(int number, int method, 
                  const std::vector<String>& atom_names,
                  const std::vector<String>& anchors,
                  const String& type, Real charge);

  void AddBond(InteractionPtr p) { this->CheckInteractionToAdd(p,"BOND");
                                     bonds_.push_back(p); }

  void AddAngle(InteractionPtr p) { this->CheckInteractionToAdd(p,"ANGLE");
                                      angles_.push_back(p); }

  void AddDihedral(InteractionPtr p) { this->CheckInteractionToAdd(p,"DIHEDRAL");
                                         dihedrals_.push_back(p); }

  void AddImproper(InteractionPtr p) { this->CheckInteractionToAdd(p,"IMPROPER");
                                         impropers_.push_back(p); }

  void AddCMap(InteractionPtr p) { this->CheckInteractionToAdd(p,"CMAP");
                                     cmaps_.push_back(p); }

  void AddDeleteAtom(const String& atom_name) { delete_atom_names_.push_back(atom_name); }

  // GET STUFF (for data extraction)
  const std::vector<InteractionPtr>& GetBonds() const { return bonds_; }
  const std::vector<InteractionPtr>& GetAngles() const { return angles_; }
  const std::vector<InteractionPtr>& GetDihedrals() const { return dihedrals_; }
  const std::vector<InteractionPtr>& GetImpropers() const { return impropers_; }
  const std::vector<InteractionPtr>& GetCmaps() const { return cmaps_; }
  const std::vector<String>& GetDeleteAtoms() const { return delete_atom_names_; }
  // all idx for replace rules < GetNumReplaceRules()
  uint GetNumReplaceRules() const { return replace_old_atom_name_.size(); }
  const String& GetReplaceRuleName(uint idx) const {
    return replace_old_atom_name_[idx];
  }
  const String& GetReplaceRuleNewName(uint idx) const {
    return replace_new_atom_name_[idx];
  }
  const String& GetReplaceRuleNewType(uint idx) const {
    return replace_new_atom_type_[idx];
  }
  Real GetReplaceRuleNewCharge(uint idx) const {
    return replace_new_charge_[idx];
  }
  // all idx for add rules < GetNumAddRules()
  uint GetNumAddRules() const { return add_add_number_.size(); }
  int GetAddRuleNumber(uint idx) const { return add_add_number_[idx]; }
  int GetAddRuleMethod(uint idx) const { return add_methods_[idx]; }
  const std::vector<String>& 
  GetAddRuleNames(uint idx) const { return add_atom_names_[idx]; }
  const std::vector<String>& 
  GetAddRuleAnchors(uint idx) const { return add_anchor_atom_names_[idx]; }
  const String& GetAddRuleType(uint idx) const { return add_atom_types_[idx]; }
  Real GetAddRuleCharge(uint idx) const { return add_charges_[idx]; }

  virtual void OnSave(ost::io::BinaryDataSink& ds) { ds << *this; }

  virtual BlockModifierType GetBlockModifierType() { return GromacsBlockModifiers; }

  template <typename DS>
  void Serialize(DS& ds){

    int num_bonds = bonds_.size();
    int num_angles = angles_.size();
    int num_dihedrals = dihedrals_.size();
    int num_impropers = impropers_.size();
    int num_cmaps = cmaps_.size();

    ds & num_bonds;
    ds & num_angles;
    ds & num_dihedrals;
    ds & num_impropers;
    ds & num_cmaps;

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

    int num_replace_atoms = replace_old_atom_name_.size();
    ds & num_replace_atoms;

    if(ds.IsSource()){
      replace_old_atom_name_ = std::vector<String>(num_replace_atoms);
      replace_new_atom_name_ = std::vector<String>(num_replace_atoms);
      replace_new_atom_type_ = std::vector<String>(num_replace_atoms);
      replace_new_charge_ = std::vector<Real>(num_replace_atoms);
    }

    for(int i = 0; i < num_replace_atoms; ++i){
      ds & replace_old_atom_name_[i];
      ds & replace_new_atom_name_[i];
      ds & replace_new_atom_type_[i];
      ds & replace_new_charge_[i];
    }

    int num_add_atoms = add_add_number_.size();
    ds & num_add_atoms;

    if(ds.IsSource()){
      add_add_number_ = std::vector<int>(num_add_atoms);
      add_methods_ = std::vector<int>(num_add_atoms);
      add_atom_names_ = std::vector<std::vector<String> >(num_add_atoms);
      add_anchor_atom_names_ = std::vector<std::vector<String> >(num_add_atoms);
      add_atom_types_ = std::vector<String>(num_add_atoms);
      add_charges_ = std::vector<Real>(num_add_atoms);
    }

    for(int i = 0; i < num_add_atoms; ++i){
      ds & add_add_number_[i];
      ds & add_methods_[i];
      ds & add_atom_types_[i];
      ds & add_charges_[i];
      int num_add_atom_names = add_atom_names_[i].size();
      int num_add_anchor_atom_names = add_anchor_atom_names_[i].size();
      ds & num_add_atom_names;
      ds & num_add_anchor_atom_names;
      if(ds.IsSource()){
        add_atom_names_[i] = std::vector<String>(num_add_atom_names);
        add_anchor_atom_names_[i] = std::vector<String>(num_add_anchor_atom_names);
      }
      for(int j = 0; j < num_add_atom_names; ++j){
        ds & add_atom_names_[i][j];
      }
      for(int j = 0; j < num_add_anchor_atom_names; ++j){
        ds & add_anchor_atom_names_[i][j];
      }
    }

    int num_delete_atoms = delete_atom_names_.size();
    ds & num_delete_atoms;
    if(ds.IsSource()){
      delete_atom_names_ = std::vector<String>(num_delete_atoms);
    }
    for(int i = 0; i < num_delete_atoms; ++i){
      ds & delete_atom_names_[i];
    }
  }

private:

  void CheckInteractionToAdd(InteractionPtr p, const String& interaction_type) const;

  //Atoms, that will be added to the building block

  //Interactions, that will be added to the building block
  std::vector<InteractionPtr> bonds_;
  std::vector<InteractionPtr> angles_;
  std::vector<InteractionPtr> dihedrals_;
  std::vector<InteractionPtr> impropers_;
  std::vector<InteractionPtr> cmaps_;

  //paramters we need for the replace statement
  std::vector<String> replace_old_atom_name_;
  std::vector<String> replace_new_atom_name_;
  std::vector<String> replace_new_atom_type_;
  std::vector<Real> replace_new_charge_;

  //parameters we need for the add statement => gromacs dependent rules!
  std::vector<int> add_add_number_;
  std::vector<int> add_methods_;
  std::vector<std::vector<String> > add_atom_names_;
  std::vector<std::vector<String> > add_anchor_atom_names_;
  std::vector<String> add_atom_types_;
  std::vector<Real> add_charges_;

  //parameters we need for the delete statement
  std::vector<String> delete_atom_names_;
};




}}}

#endif
