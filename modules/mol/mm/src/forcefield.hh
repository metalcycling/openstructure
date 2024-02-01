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

#ifndef OST_MM_FORCE_FIELD_HH
#define OST_MM_FORCE_FIELD_HH

#include <vector>
#include <algorithm>


#include <boost/shared_ptr.hpp>
#include <boost/unordered_map.hpp>

#include <ost/message.hh>
#include <ost/mol/mm/buildingblock.hh>
#include <ost/mol/mm/interaction.hh>
#include <ost/mol/mm/block_modifiers.hh>
#include <ost/mol/mm/gromacs_block_modifiers.hh>

namespace ost { namespace mol{ namespace mm{

struct ResidueNames;
class Forcefield;

typedef boost::shared_ptr<ResidueNames> ResidueNamesPtr;
typedef boost::shared_ptr<ost::mol::mm::Forcefield> ForcefieldPtr;

struct ResidueNames{

  ResidueNames() { };

  ResidueNames(String a, String b, String c, String d):
               main(a),nter(b),cter(c),twoter(d) { }

  String main;
  String nter;
  String cter;
  String twoter;

  bool Contains(const String& name){
    return name == main || name == nter || name == cter || name == twoter;
  }

  template <typename DS>
  void Serialize(DS& ds){
    ds & main;
    ds & nter;
    ds & cter;
    ds & twoter;
  }
};


class Forcefield {
public:

  Forcefield(): gen_pairs_(true), fudge_LJ_(1.0),fudge_QQ_(1.0) { }

  static ForcefieldPtr Load(const String& filename);

  void Save(const String& filename);

  //Getter functions

  BuildingBlockPtr GetBuildingBlock(const String& name) const;

  BlockModifierPtr GetBlockModifier(const String& modifier_name) const;

  std::vector<String> GetBuildingBlockNames() const;

  String GetAtomType(const String& res_name, const String& atom_name) const;

  HydrogenConstructorPtr GetHydrogenConstructor(const String& name) const;

  BlockModifierPtr GetNTerModifier(const String& res_name, const String& ter_name = "") const;

  BlockModifierPtr GetCTerModifier(const String& res_name, const String& ter_name = "") const;

  InteractionPtr GetBond(const String& type1, 
                           const String& type2) const;

  InteractionPtr GetAngle(const String& type1, 
                            const String& type2,
                            const String& type3) const;

  std::vector<InteractionPtr> GetDihedrals(const String& type1,
                                             const String& type2,
                                             const String& type3,
                                             const String& type4) const;

  std::vector<InteractionPtr> GetImpropers(const String& type1,
                                             const String& type2,
                                             const String& type3,
                                             const String& type4) const;

  InteractionPtr GetCMap(const String& type1,
                           const String& type2,
                           const String& type3,
                           const String& type4,
                           const String& type5) const;

  InteractionPtr GetImplicitGenborn(const String& type1) const;

  InteractionPtr GetLJ(const String& type1, 
                             const String& type2, 
                             bool pair=false) const;

  InteractionPtr GetLJ(const String& type) const;

  InteractionPtr GetConstraint(const String& type1,
                                 const String& type2);

  Real GetMass(const String& type) const;

  Real GetFudgeLJ() const { return fudge_LJ_; }

  Real GetFudgeQQ() const { return fudge_QQ_; }

  //functions to add interactions and settings

  void AddBuildingBlock(const String& name, BuildingBlockPtr p) { building_blocks_[name] = p; }

  void AddBond(InteractionPtr p);

  void AddAngle(InteractionPtr p);

  void AddDihedral(InteractionPtr p);

  void AddImproper(InteractionPtr p);

  void AddCMap(InteractionPtr p);

  void AddImplicitGenborn(InteractionPtr p);

  void AddLJ(InteractionPtr p);

  void AddLJPair(InteractionPtr p);

  void AddConstraint(InteractionPtr p);

  void AddMass(const String& type, Real mass) { atom_masses_[type] = mass; }

  void SetFudgeLJ(Real f_lj) { fudge_LJ_ = f_lj; }

  void SetFudgeQQ(Real f_qq) { fudge_QQ_ = f_qq; }

  void SetGenPairs(bool gen_pairs) { gen_pairs_ = gen_pairs; }

  void AddResidueRenamingRule(const String& name,
                              const String& ff_main_name,
                              const String& ff_n_ter_name,
                              const String& ff_c_ter_name,
                              const String& ff_two_ter_name);

  void AddAtomRenamingRule(const String& res_name,
                           const String& old_atom_name,
                           const String& new_atom_name);

  String GetResidueRenamingMain(const String& name) const;

  String GetResidueRenamingNTer(const String& name) const;

  String GetResidueRenamingCTer(const String& name) const;

  String GetResidueRenamingTwoTer(const String& name) const;

  String GetAtomRenaming(const String& res_name, const String& atom_name) const;

  void AddHydrogenConstructor(const String& residue_name, HydrogenConstructorPtr p){
    hydrogen_constructors_[residue_name] = p;
  }

  void AddBlockModifier(const String& modifier_name,
                        BlockModifierPtr p) { block_modifiers_[modifier_name] = p; }

  void SetStandardCTer(const String& res_name, const String& ter_name) { standard_c_termini_[res_name] = ter_name; }

  void SetStandardNTer(const String& res_name, const String& ter_name) { standard_n_termini_[res_name] = ter_name; }


  //Renaming to the forcefield specific names (residues/atoms)
  void AssignFFSpecificNames(ost::mol::EntityHandle& handle, bool reverse=false) const;

  // get renaming rules (for data extraction)
  bool HasAtomRenamingRules(const String& res_name) const {
    return (atom_renaming_ff_specific_.find(res_name)
            != atom_renaming_ff_specific_.end());
  }
  typedef std::vector<std::pair<String,String> > AtomRenamingType;
  const AtomRenamingType& GetAtomRenamingRules(const String& res_name) const;

private:

  String AtomTypesToKeyword(std::vector<String>& types, bool allow_reordering = true) const;
  void CheckInteractionToAdd(InteractionPtr p, const String& interaction_type) const;

  //this is all nonbonded stuff
  bool gen_pairs_;
  Real fudge_LJ_;
  Real fudge_QQ_;

  boost::unordered_map<String, Real> atom_masses_;
  boost::unordered_map<String, BuildingBlockPtr> building_blocks_;
  boost::unordered_map<String, BlockModifierPtr> block_modifiers_;

  //the standard interactions
  boost::unordered_map<String,InteractionPtr> bonds_;
  boost::unordered_map<String,InteractionPtr> angles_;
  boost::unordered_map<String,InteractionPtr> lj_14_pairs_; 
  boost::unordered_map<String,InteractionPtr> constraints_;
  boost::unordered_map<String,InteractionPtr> cmaps_;
  boost::unordered_map<String,InteractionPtr> implicit_genborn_;
  boost::unordered_map<String,InteractionPtr> ljs_; 
  boost::unordered_map<String,std::vector<InteractionPtr> > dihedrals_;
  boost::unordered_map<String,std::vector<InteractionPtr> > improper_dihedrals_;

  boost::unordered_map<String, AtomRenamingType> atom_renaming_ff_specific_;
  boost::unordered_map<String, ResidueNamesPtr> res_renaming_ff_specific_;

  boost::unordered_map<String, HydrogenConstructorPtr> hydrogen_constructors_;
  boost::unordered_map<String, String> standard_n_termini_;
  boost::unordered_map<String, String> standard_c_termini_;
};

}}}//ns

#endif
