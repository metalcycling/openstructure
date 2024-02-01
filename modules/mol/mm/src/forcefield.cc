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

#include <ost/mol/mm/forcefield.hh>

namespace ost{ namespace mol{ namespace mm{

ForcefieldPtr Forcefield::Load(const String& filename){


  if (!boost::filesystem::exists(filename)) {
    std::stringstream ss;
    ss << "Could not open forcefield. File '"
       << filename << "' does not exist";
    throw ost::io::IOException(ss.str());
  }

  std::ifstream stream(filename.c_str(), std::ios_base::binary);
  ost::io::BinaryDataSource ds(stream);
  ForcefieldPtr ff_p(new Forcefield);

  bool gen_pairs;
  Real fudge_LJ;
  Real fudge_QQ;
  ds & gen_pairs;
  ds & fudge_LJ;
  ds & fudge_QQ;

  ff_p->SetGenPairs(gen_pairs);
  ff_p->SetFudgeLJ(fudge_LJ);
  ff_p->SetFudgeQQ(fudge_QQ);

  int num_atom_masses;
  ds & num_atom_masses;

  String atom_name;
  Real mass;
  for(int i = 0; i < num_atom_masses; ++i){
    ds & atom_name;
    ds & mass;
    ff_p->AddMass(atom_name,mass);
  }

  int num_building_blocks;
  ds & num_building_blocks;
  String block_name;
  for(int i = 0; i < num_building_blocks; ++i){
    ds & block_name;
    BuildingBlockPtr p(new BuildingBlock);
    ds & *p;
    ff_p->AddBuildingBlock(block_name,p);
  }
  
  int func_type;
  int num_bonds;
  ds & num_bonds;
  for(int i = 0; i < num_bonds; ++i){
    ds & func_type;
    InteractionPtr p(new Interaction(FuncType(func_type)));
    ds & *p;
    ff_p->AddBond(p);
  }

  int num_angles;
  ds & num_angles;
  for(int i = 0; i < num_angles; ++i){
    ds & func_type;
    InteractionPtr p(new Interaction(FuncType(func_type)));
    ds & *p;
    ff_p->AddAngle(p);
  }
  
  int num_lj_pairs;
  ds & num_lj_pairs;
  for(int i = 0; i < num_lj_pairs; ++i){
    ds & func_type;
    InteractionPtr p(new Interaction(FuncType(func_type)));
    ds & *p;
    ff_p->AddLJPair(p);
  }

  int num_constraints;
  ds & num_constraints;
  for(int i = 0; i < num_constraints; ++i){
    ds & func_type;
    InteractionPtr p(new Interaction(FuncType(func_type)));
    ds & *p;
    ff_p->AddConstraint(p);
  }
  
  int num_cmaps;
  ds & num_cmaps;
  for(int i = 0; i < num_cmaps; ++i){
    ds & func_type;
    InteractionPtr p(new Interaction(FuncType(func_type)));
    ds & *p;
    ff_p->AddCMap(p);
  }

  int num_genborn;
  ds & num_genborn;
  for(int i = 0; i < num_genborn; ++i){
    ds & func_type;
    InteractionPtr p(new Interaction(FuncType(func_type)));
    ds & *p;
    ff_p->AddImplicitGenborn(p);
  }

  int num_ljs;
  ds & num_ljs;
  for(int i = 0; i < num_ljs; ++i){
    ds & func_type;
    InteractionPtr p(new Interaction(FuncType(func_type)));
    ds & *p;
    ff_p->AddLJ(p);
  }


  int num_dihedrals;
  ds & num_dihedrals;
  int temp_num;
  for(int i = 0; i < num_dihedrals; ++i){
    ds & temp_num;
    for(int j = 0; j < temp_num; ++j){
      ds & func_type;
      InteractionPtr p(new Interaction(FuncType(func_type)));
      ds & *p;
      ff_p->AddDihedral(p);
    }
  }

  int num_impropers;
  ds & num_impropers;
  for(int i = 0; i < num_impropers; ++i){
    ds & temp_num;
    for(int j = 0; j < temp_num; ++j){
      ds & func_type;
      InteractionPtr p(new Interaction(FuncType(func_type)));
      ds & *p;
      ff_p->AddImproper(p);
    }
  }

  int num_atom_renamings;
  ds & num_atom_renamings;
  String res_name, old_name, new_name;
  for(int i = 0; i < num_atom_renamings; ++i){
    ds & res_name;
    ds & temp_num;
    for(int j = 0; j < temp_num; ++j){
      ds & old_name;
      ds & new_name;
      ff_p->AddAtomRenamingRule(res_name,old_name,new_name);
    }
  }

  int num_res_renaming;
  ds & num_res_renaming;
  for(int i = 0; i < num_res_renaming; ++i){
    ds & res_name;
    ResidueNamesPtr p(new ResidueNames);
    ds & *p;
    ff_p->AddResidueRenamingRule(res_name, p->main, p->nter, p->cter, p->twoter);
  }

  int num_hydrogen_constructors;
  int block_modifier_type;
  String block_modifier_name;
  ds & num_hydrogen_constructors;
  for(int i = 0; i < num_hydrogen_constructors; ++i){
    ds & block_modifier_type;
    ds & block_modifier_name;
    switch(block_modifier_type){
      case GromacsBlockModifiers:{
        GromacsHydrogenConstructorPtr p(new GromacsHydrogenConstructor);
        ds & *p;
        ff_p->AddHydrogenConstructor(block_modifier_name,p);
        break;
      }
      default:{
        throw ost::Error("Invalid Block Modifier type encountered when loading hydrogen constructor!");
      }

    }
  }

  int num_block_modifiers;
  ds & num_block_modifiers;
  String modifier_name;

  for(int i = 0; i < num_block_modifiers; ++i){
    ds & block_modifier_type;
    switch(block_modifier_type){
      case GromacsBlockModifiers:{
        GromacsBlockModifierPtr p(new GromacsBlockModifier);
        ds & modifier_name;
        ds & *p;
        ff_p->AddBlockModifier(modifier_name,p);
        break;
      }
      default:{
        throw ost::Error("Invalid Block Modifier type encountered when loading!");
      }
    }
    
  }

  int num_standard_n_termini;
  int num_standard_c_termini;

  ds & num_standard_n_termini;
  ds & num_standard_c_termini;

  for(int i = 0; i < num_standard_n_termini; ++i){
    ds & res_name;
    ds & modifier_name;
    ff_p->SetStandardNTer(res_name,modifier_name);
  }

  for(int i = 0; i < num_standard_c_termini; ++i){
    ds & res_name;
    ds & modifier_name;
    ff_p->SetStandardCTer(res_name,modifier_name);
  }

  return ff_p;
}

void Forcefield::Save(const String& filename){
  std::ofstream stream(filename.c_str(), std::ios_base::binary);
  io::BinaryDataSink ds(stream);


  ds & gen_pairs_;
  ds & fudge_LJ_;
  ds & fudge_QQ_;

  int num_atom_masses = atom_masses_.size();
  ds & num_atom_masses;
  for(boost::unordered_map<String,Real>::iterator i = atom_masses_.begin();
      i != atom_masses_.end(); ++i){
    ds & i->first;
    ds & i->second;
  }
  
  int num_building_blocks = building_blocks_.size();
  ds & num_building_blocks;
  for(boost::unordered_map<String,BuildingBlockPtr>::iterator i = building_blocks_.begin();
      i != building_blocks_.end(); ++i){
    ds & i->first;
    ds & *(i->second);
  }

  
  int func_type;
  int num_bonds = bonds_.size();
  ds & num_bonds;
  for(boost::unordered_map<String,InteractionPtr>::iterator i = bonds_.begin();
      i != bonds_.end(); ++i){
    func_type = int(i->second->GetFuncType());
    ds & func_type;
    ds & *(i->second);
  }
  
  int num_angles = angles_.size();
  ds & num_angles;
  for(boost::unordered_map<String,InteractionPtr>::iterator i = angles_.begin();
      i != angles_.end(); ++i){
    func_type = int(i->second->GetFuncType());
    ds & func_type;
    ds & *(i->second);
  }

  int num_lj_pairs = lj_14_pairs_.size();
  ds & num_lj_pairs;
  for(boost::unordered_map<String,InteractionPtr>::iterator i = lj_14_pairs_.begin();
      i != lj_14_pairs_.end(); ++i){
    func_type = int(i->second->GetFuncType());
    ds & func_type;
    ds & *(i->second);
  }
  
  int num_constraints = constraints_.size();
  ds & num_constraints;
  for(boost::unordered_map<String,InteractionPtr>::iterator i = constraints_.begin();
      i != constraints_.end(); ++i){
    func_type = int(i->second->GetFuncType());
    ds & func_type;
    ds & *(i->second);
  }

  int num_cmaps = cmaps_.size();
  ds & num_cmaps;
  for(boost::unordered_map<String,InteractionPtr>::iterator i = cmaps_.begin();
      i != cmaps_.end(); ++i){
    func_type = int(i->second->GetFuncType());
    ds & func_type;
    ds & *(i->second);
  }
  
  int num_genborn = implicit_genborn_.size();
  ds & num_genborn;
  for(boost::unordered_map<String,InteractionPtr>::iterator i = implicit_genborn_.begin();
      i != implicit_genborn_.end(); ++i){
    func_type = int(i->second->GetFuncType());
    ds & func_type;
    ds & *(i->second);
  }
  
  int num_ljs = ljs_.size();
  ds & num_ljs;
  for(boost::unordered_map<String,InteractionPtr>::iterator i = ljs_.begin();
      i != ljs_.end(); ++i){
    func_type = int(i->second->GetFuncType());
    ds & func_type;
    ds & *(i->second);
  }
  
  int num_dihedrals = dihedrals_.size();
  int vector_size;
  ds & num_dihedrals;
  for(boost::unordered_map<String,std::vector<InteractionPtr> >::iterator i = dihedrals_.begin();
      i != dihedrals_.end(); ++i){
    vector_size = i->second.size();
    ds & vector_size;
    for(std::vector<InteractionPtr>::iterator j = i->second.begin();
        j != i->second.end(); ++j){
      func_type = int((*j)->GetFuncType());
      ds & func_type;
      ds & **j;
    }
  }

  int num_impropers = improper_dihedrals_.size();
  ds & num_impropers;
  for(boost::unordered_map<String,std::vector<InteractionPtr> >::iterator i = improper_dihedrals_.begin();
      i != improper_dihedrals_.end(); ++i){
    vector_size = i->second.size();
    ds & vector_size;
    for(std::vector<InteractionPtr>::iterator j = i->second.begin();
        j != i->second.end(); ++j){
      func_type = int((*j)->GetFuncType());
      ds & func_type;
      ds & **j;
    }
  }
  

  int num_atom_renamings = atom_renaming_ff_specific_.size();
  ds & num_atom_renamings;
  for(boost::unordered_map<String,std::vector<std::pair<String,String> > >::iterator i = atom_renaming_ff_specific_.begin();
      i != atom_renaming_ff_specific_.end(); ++i){
    ds & i->first;
    vector_size = i->second.size();
    ds & vector_size;
    for(std::vector<std::pair<String,String> >::iterator j = i->second.begin();
        j != i->second.end(); ++j){
      ds & j->first;
      ds & j->second;
    }
  }

  int num_res_renaming = res_renaming_ff_specific_.size();
  ds & num_res_renaming;
  for(boost::unordered_map<String,ResidueNamesPtr>::iterator i = res_renaming_ff_specific_.begin();
      i != res_renaming_ff_specific_.end(); ++i){
    ds & i->first;
    ds & *(i->second);
  }

  int num_hydrogen_constructors = hydrogen_constructors_.size();
  ds & num_hydrogen_constructors;
  for(boost::unordered_map<String,HydrogenConstructorPtr>::iterator i = hydrogen_constructors_.begin();
      i != hydrogen_constructors_.end(); ++i){
    func_type = int(i->second->GetBlockModifierType());
    ds & func_type;
    ds & i->first;
    i->second->OnSave(ds);
  }

  int num_block_modifiers = block_modifiers_.size();
  ds & num_block_modifiers;
  for(boost::unordered_map<String,BlockModifierPtr>::iterator i = block_modifiers_.begin();
      i != block_modifiers_.end(); ++i){
    func_type = int(i->second->GetBlockModifierType());
    ds & func_type;
    ds & i->first;
    i->second->OnSave(ds);
  }

  int num_standard_n_termini = standard_c_termini_.size();
  int num_standard_c_termini = standard_n_termini_.size();
  ds & num_standard_n_termini;
  ds & num_standard_c_termini;

  for(boost::unordered_map<String, String>::iterator i = standard_n_termini_.begin();
      i != standard_n_termini_.end(); ++i){
    ds & i->first;
    ds & i->second;
  }

  for(boost::unordered_map<String, String>::iterator i = standard_c_termini_.begin();
      i != standard_c_termini_.end(); ++i){
    ds & i->first;
    ds & i->second;
  }

}

BuildingBlockPtr Forcefield::GetBuildingBlock(const String& name) const{
  boost::unordered_map<String,BuildingBlockPtr>::const_iterator i = building_blocks_.find(name);
  if(i == building_blocks_.end()){
    std::stringstream ss;
    ss << "Forcefield does not contain requested building block ";
    ss << "\"" << name << "\"";
    throw ost::Error(ss.str());
  }
  BuildingBlockPtr ff_p = i->second;
  //let's create a new building block to avoid messing around with the
  //force field building blocks
  BuildingBlockPtr p(new BuildingBlock(*ff_p));
  return p;
}

std::vector<String> Forcefield::GetBuildingBlockNames() const {

  std::vector<String> return_vec;
  for(boost::unordered_map<String,BuildingBlockPtr>::const_iterator i = 
      building_blocks_.begin(); i != building_blocks_.end(); ++i) {
    return_vec.push_back(i->first);
  }

  return return_vec;
} 

BlockModifierPtr Forcefield::GetBlockModifier(const String& modifier_name) const{
  boost::unordered_map<String,BlockModifierPtr>::const_iterator i = block_modifiers_.find(modifier_name);
  if(i == block_modifiers_.end()) return BlockModifierPtr();
  return i->second;
}

String Forcefield::GetAtomType(const String& res_name, const String& atom_name) const{
  BuildingBlockPtr block = this->GetBuildingBlock(res_name);
  if(!block){
    std::stringstream ss;
    ss << "There is no building block for " << res_name;
    ss << "in the forcefield!";
    throw ost::Error(ss.str());
  }
  return block->GetType(atom_name);
}

HydrogenConstructorPtr Forcefield::GetHydrogenConstructor(const String& name) const{
  boost::unordered_map<String,HydrogenConstructorPtr>::const_iterator i = hydrogen_constructors_.find(name);
  if(i != hydrogen_constructors_.end()) return i->second;
  return HydrogenConstructorPtr();
}

BlockModifierPtr Forcefield::GetNTerModifier(const String& res_name, const String& ter_name) const{
  String temp = ter_name;
  if(temp == ""){
    boost::unordered_map<String,String>::const_iterator it = standard_n_termini_.find(res_name);
    if(it != standard_n_termini_.end()) temp = it->second;
    else return BlockModifierPtr();
  }
  return this->GetBlockModifier(temp);
}

BlockModifierPtr Forcefield::GetCTerModifier(const String& res_name, const String& ter_name) const{
  String temp = ter_name;
  if(temp == ""){
    boost::unordered_map<String,String>::const_iterator it = standard_c_termini_.find(res_name);
    if(it != standard_c_termini_.end()) temp = it->second;
    else return BlockModifierPtr();
  }
  return this->GetBlockModifier(temp);
}

InteractionPtr Forcefield::GetBond(const String& type1, const String& type2) const{
  std::vector<String> types;
  types.push_back(type1);
  types.push_back(type2);
  boost::unordered_map<String,InteractionPtr>::const_iterator i = 
    bonds_.find(this->AtomTypesToKeyword(types));
  if(i != bonds_.end()) return i->second;
  throw ost::Error("Forcefield does not contain bond for requested atom types!");
}

InteractionPtr Forcefield::GetAngle(const String& type1, const String& type2,
                                      const String& type3) const{
  std::vector<String> types;
  types.push_back(type1);
  types.push_back(type2);
  types.push_back(type3);
  boost::unordered_map<String,InteractionPtr>::const_iterator i = 
    angles_.find(this->AtomTypesToKeyword(types));
  if(i != angles_.end()) return i->second;
  throw ost::Error("Forcefield does not contain angle for requested atom types!");
}

std::vector<InteractionPtr> Forcefield::GetDihedrals(const String& type1,
                                                       const String& type2,
                                                       const String& type3,
                                                       const String& type4) const{
  std::vector<String> types;
  types.push_back(type1);
  types.push_back(type2);
  types.push_back(type3);
  types.push_back(type4);

  boost::unordered_map<String,std::vector<InteractionPtr> >::const_iterator i = 
    dihedrals_.find(this->AtomTypesToKeyword(types));

  if(i != dihedrals_.end()) return i->second;

  //if we reach this point, we have to look for all dihedrals with wildcards...
  std::vector<InteractionPtr> dihedrals;
  std::vector<String> types_with_wildcards;
  
  //lets first check single wildcards
  for(int i = 0; i < 4; ++i){
    types_with_wildcards = types;
    types_with_wildcards[i] = "X";
    boost::unordered_map<String,std::vector<InteractionPtr> >::const_iterator j = 
      dihedrals_.find(this->AtomTypesToKeyword(types_with_wildcards));
    if(j != dihedrals_.end()) dihedrals.insert(dihedrals.end(),j->second.begin(),j->second.end());
  }

  //we check for two wildcards
  for(int i = 0; i < 4; ++i){
    for(int j = i+1; j < 4; ++j){
      types_with_wildcards = types;
      types_with_wildcards[i] = "X";
      types_with_wildcards[j] = "X";
      boost::unordered_map<String,std::vector<InteractionPtr> >::const_iterator k = 
        dihedrals_.find(this->AtomTypesToKeyword(types_with_wildcards));
      if(k != dihedrals_.end()) dihedrals.insert(dihedrals.end(),k->second.begin(),k->second.end());
    }
  }

  //as stupid as it may sound, we also look for dihedrals with 3 wildcards
  std::vector<String> all_wildcards;
  for(int i = 0; i < 4; ++i){
    all_wildcards.push_back("X");
  }
  for(int i = 0; i < 4; ++i){
    types_with_wildcards = all_wildcards;
    types_with_wildcards[i] = types[i];
    boost::unordered_map<String,std::vector<InteractionPtr> >::const_iterator j = 
      dihedrals_.find(this->AtomTypesToKeyword(types_with_wildcards));
    if(j != dihedrals_.end()) dihedrals.insert(dihedrals.end(),j->second.begin(),j->second.end());
  }  

  if(!dihedrals.empty()) return dihedrals;
  throw ost::Error("Forcefield does not contain any dihedrals for requested atom types!");
}

std::vector<InteractionPtr> Forcefield::GetImpropers(const String& type1,
                                                       const String& type2,
                                                       const String& type3,
                                                       const String& type4) const{
  std::vector<String> types;
  types.push_back(type1);
  types.push_back(type2);
  types.push_back(type3);
  types.push_back(type4);

  boost::unordered_map<String,std::vector<InteractionPtr> >::const_iterator i = 
    improper_dihedrals_.find(this->AtomTypesToKeyword(types));

  if(i != improper_dihedrals_.end()) return i->second;

  //if we reach this point, we have to look for all dihedrals with wildcards...
  std::vector<InteractionPtr> impropers;
  std::vector<String> types_with_wildcards;
  
  //lets first check single wildcards
  for(int i = 0; i < 4; ++i){
    types_with_wildcards = types;
    types_with_wildcards[i] = "X";
    boost::unordered_map<String,std::vector<InteractionPtr> >::const_iterator j = 
      improper_dihedrals_.find(this->AtomTypesToKeyword(types_with_wildcards));
    if(j != improper_dihedrals_.end()) impropers.insert(impropers.end(),j->second.begin(),j->second.end());
  }

  //we check for two wildcards
  for(int i = 0; i < 4; ++i){
    for(int j = i+1; j < 4; ++j){
      types_with_wildcards = types;
      types_with_wildcards[i] = "X";
      types_with_wildcards[j] = "X";
      boost::unordered_map<String,std::vector<InteractionPtr> >::const_iterator k = 
        improper_dihedrals_.find(this->AtomTypesToKeyword(types_with_wildcards));
      if(k != improper_dihedrals_.end()) impropers.insert(impropers.end(),k->second.begin(),k->second.end());
    }
  }

  //as stupid as it may sound, we also look for dihedrals with 3 wildcards
  std::vector<String> all_wildcards;
  for(int i = 0; i < 4; ++i){
    all_wildcards.push_back("X");
  }
  for(int i = 0; i < 4; ++i){
    types_with_wildcards = all_wildcards;
    types_with_wildcards[i] = types[i];
    boost::unordered_map<String,std::vector<InteractionPtr> >::const_iterator j = 
      improper_dihedrals_.find(this->AtomTypesToKeyword(types_with_wildcards));
    if(j != improper_dihedrals_.end()) impropers.insert(impropers.end(),j->second.begin(),j->second.end());
  }  

  if(!impropers.empty()) return impropers;

  throw ost::Error("Forcefield does not contain any impropers for requested atom types!");
}

InteractionPtr Forcefield::GetCMap(const String& type1, const String& type2,
                                     const String& type3, const String& type4, 
                                     const String& type5) const{

  std::vector<String> types;
  types.push_back(type1);
  types.push_back(type2);
  types.push_back(type3);
  types.push_back(type4);
  types.push_back(type5);
  boost::unordered_map<String,InteractionPtr>::const_iterator i = 
    cmaps_.find(this->AtomTypesToKeyword(types,false));
  if(i != cmaps_.end()) return i->second;

  throw ost::Error("Forcefield does not contain cmap for requested atom types!");
}

InteractionPtr Forcefield::GetImplicitGenborn(const String& type) const{

  boost::unordered_map<String,InteractionPtr>::const_iterator i = 
    implicit_genborn_.find(type);
  if(i != implicit_genborn_.end()) return i->second;
  throw ost::Error("Forcefield does not contain implicit genborn parameters for requested atom type!");
}

InteractionPtr Forcefield::GetLJ(const String& type) const{

  boost::unordered_map<String,InteractionPtr>::const_iterator i = 
    ljs_.find(type);

  if(i != ljs_.end()){
    return i->second;
  }
  throw ost::Error("Forcefield does not contain lennard jones parameters for requested atom types!");
}

InteractionPtr Forcefield::GetLJ(const String& type1, 
                                       const String& type2,
                                       bool pair) const{
  std::vector<String> types;
  types.push_back(type1);
  types.push_back(type2);

  if(pair){
    boost::unordered_map<String,InteractionPtr>::const_iterator i = 
      lj_14_pairs_.find(this->AtomTypesToKeyword(types));
    if(i == lj_14_pairs_.end()){
      if(!gen_pairs_){
        std::stringstream ss;
        ss << "was not able to find lj parameters between atoms of type ";
        ss << type1 << " and " << type2 << " with gen_pairs ";
        ss << "set to no!";
        throw ost::Error(ss.str());
      }
    }
    else{
      return i->second;
    }
  }
  
  InteractionPtr temp1, temp2;
  temp1 = this->GetLJ(type1);
  temp2 = this->GetLJ(type2);
  if(temp1 && temp2){
    InteractionPtr return_ptr(new Interaction(LJPair));
    std::vector<Real> param1, param2, param;
    param1 = temp1->GetParam();
    param2 = temp2->GetParam();
    param.push_back(0.5*(param1[0]+param2[0]));
    param.push_back(sqrt(param1[1]*param2[1]));
    return_ptr->SetTypes(types);
    return_ptr->SetParam(param);
    return return_ptr;
  }
  
  throw ost::Error("Forcefield does not contain lennard jones parameters for requested atom types!");

}

InteractionPtr Forcefield::GetConstraint(const String& type1,
                                           const String& type2){
  std::vector<String> types;
  types.push_back(type1);
  types.push_back(type2);
  boost::unordered_map<String,InteractionPtr>::const_iterator i = 
    constraints_.find(this->AtomTypesToKeyword(types));
  if(i != constraints_.end()) return i->second;
  throw ost::Error("Forcefield does not contain constraint for requested atom types!");
}

Real Forcefield::GetMass(const String& atom_type) const{
  boost::unordered_map<String,Real>::const_iterator it = atom_masses_.find(atom_type);
  if(it == atom_masses_.end()){
    std::stringstream ss;
    ss << "There is no mass defined for atom of type " << atom_type;
    throw ost::Error(ss.str());
  }
  return it->second;
}

void Forcefield::AddBond(InteractionPtr p){
  this->CheckInteractionToAdd(p,"BOND");
  std::vector<String> types = p->GetTypes();
  bonds_[this->AtomTypesToKeyword(types)] = p;
}

void Forcefield::AddAngle(InteractionPtr p){
  this->CheckInteractionToAdd(p,"ANGLE");
  std::vector<String> types = p->GetTypes();
  angles_[this->AtomTypesToKeyword(types)] = p;
}

void Forcefield::AddDihedral(InteractionPtr p){
  this->CheckInteractionToAdd(p,"DIHEDRAL");
  std::vector<String> types = p->GetTypes();
  String keyword = this->AtomTypesToKeyword(types);
  if(dihedrals_.find(keyword) == dihedrals_.end()){
    dihedrals_[keyword] = std::vector<InteractionPtr>();
  }
  dihedrals_[keyword].push_back(p);
}

void Forcefield::AddImproper(InteractionPtr p){
  this->CheckInteractionToAdd(p,"IMPROPER");
  std::vector<String> types = p->GetTypes();
  String keyword = this->AtomTypesToKeyword(types);
  if(improper_dihedrals_.find(keyword) == improper_dihedrals_.end()){
    improper_dihedrals_[keyword] = std::vector<InteractionPtr>();
  }
  improper_dihedrals_[keyword].push_back(p);
}

void Forcefield::AddCMap(InteractionPtr p){
  this->CheckInteractionToAdd(p,"CMAP");
  std::vector<String> types = p->GetTypes();
  cmaps_[this->AtomTypesToKeyword(types,false)] = p;
}

void Forcefield::AddImplicitGenborn(InteractionPtr p){
  this->CheckInteractionToAdd(p,"IMPLICIT_GENBORN");
  std::vector<String> types = p->GetTypes();
  implicit_genborn_[this->AtomTypesToKeyword(types)] = p;
}

void Forcefield::AddLJ(InteractionPtr p){
  this->CheckInteractionToAdd(p,"LJ");
  std::vector<String> types = p->GetTypes();
  ljs_[types[0]] = p;
}

void Forcefield::AddLJPair(InteractionPtr p){
  this->CheckInteractionToAdd(p,"LJ_PAIR");
  std::vector<String> types = p->GetTypes();
  lj_14_pairs_[this->AtomTypesToKeyword(types)] = p;
}

void Forcefield::AddConstraint(InteractionPtr p){
  this->CheckInteractionToAdd(p,"CONSTRAINT");
  std::vector<String> types = p->GetTypes();
  constraints_[this->AtomTypesToKeyword(types)] = p;
}

void Forcefield::AddResidueRenamingRule(const String& name,
                                        const String& ff_main_name,
                                        const String& ff_n_ter_name,
                                        const String& ff_c_ter_name,
                                        const String& ff_two_ter_name){
  ResidueNamesPtr p(new ResidueNames(ff_main_name,
                                     ff_n_ter_name,
                                     ff_c_ter_name,
                                     ff_two_ter_name));
  res_renaming_ff_specific_[name] = p;
}

void Forcefield::AddAtomRenamingRule(const String& res_name,
                                     const String& old_atom_name,
                                     const String& new_atom_name){
  if(atom_renaming_ff_specific_.find(res_name) == atom_renaming_ff_specific_.end()){
    std::vector<std::pair<String,String> > new_vec;
    atom_renaming_ff_specific_[res_name] = new_vec;
  }

  atom_renaming_ff_specific_[res_name].push_back(std::make_pair(old_atom_name,new_atom_name));
}

String Forcefield::GetResidueRenamingMain(const String& name) const{
  boost::unordered_map<String, ResidueNamesPtr>::const_iterator i = res_renaming_ff_specific_.find(name);
  if(i == res_renaming_ff_specific_.end()){
    return name;
  }
  return i->second->main;  
}

String Forcefield::GetResidueRenamingNTer(const String& name) const{
  boost::unordered_map<String, ResidueNamesPtr>::const_iterator i = res_renaming_ff_specific_.find(name);
  if(i == res_renaming_ff_specific_.end()){
    return name;
  }
  return i->second->nter;  
}


String Forcefield::GetResidueRenamingCTer(const String& name) const{
  boost::unordered_map<String, ResidueNamesPtr>::const_iterator i = res_renaming_ff_specific_.find(name);
  if(i == res_renaming_ff_specific_.end()){
    return name;
  }
  return i->second->cter;  
}

String Forcefield::GetResidueRenamingTwoTer(const String& name) const{
  boost::unordered_map<String, ResidueNamesPtr>::const_iterator i = res_renaming_ff_specific_.find(name);
  if(i == res_renaming_ff_specific_.end()){
    return name;
  }
  return i->second->twoter;  
}

String Forcefield::GetAtomRenaming(const String& res_name,
                                   const String& atom_name) const {
  boost::unordered_map<String,AtomRenamingType>::const_iterator i;
  i = atom_renaming_ff_specific_.find(res_name);
  if(i == atom_renaming_ff_specific_.end()) return atom_name;
  for(std::vector<std::pair<String,String> >::const_iterator j = i->second.begin(), e = i->second.end();
      j != e; ++j){
    if(atom_name == j->first) return j->second; 
  } 
  return atom_name;
}

const Forcefield::AtomRenamingType&
Forcefield::GetAtomRenamingRules(const String& res_name) const {
  boost::unordered_map<String,AtomRenamingType>::const_iterator i;
  i = atom_renaming_ff_specific_.find(res_name);
  if (i == atom_renaming_ff_specific_.end()) {
    throw ost::Error(String("No atom renaming rules for ") + res_name);
  }
  return i->second;
}

void Forcefield::AssignFFSpecificNames(ost::mol::EntityHandle& handle, bool reverse) const{

  //the whole FF specific renaming procedure assumes a standard naming,
  //in our case the gromacs defined standard naming

  ost::mol::XCSEditor ed = handle.EditXCS();
  ost::mol::ResidueHandleList res_list = handle.GetResidueList();
  ost::mol::AtomHandleList atom_list;
  ost::mol::AtomHandle atom;
  ost::mol::ResidueHandle residue;

  if(!reverse){

    ost::mol::ResidueHandle next,prev;
    bool n_ter,c_ter; 
    //ost::mol::AtomHandle peptide_n,peptide_c,nucleotide_p,nucleotide_o; 

    for(ost::mol::ResidueHandleList::iterator i = res_list.begin(); 
        i != res_list.end(); ++i){

      if(res_renaming_ff_specific_.find(i->GetName()) == res_renaming_ff_specific_.end()){
        continue;
      }

      n_ter = false;
      c_ter = false;

      //peptide_n = i->FindAtom("N");
      //nucleotide_p = i->FindAtom("P");

      prev = i->GetPrev();
      next = i->GetNext();

      if(!prev.IsValid()){
        n_ter = true;
      }
      //else{
        //peptide_c = prev.FindAtom("C");
        //nucleotide_o = prev.FindAtom("O3");

        //if(!ost::mol::BondExists(peptide_n,peptide_c) && !ost::mol::BondExists(nucleotide_p,nucleotide_o)){
          //n_ter = true;
        //}
      //}

      if(!next.IsValid()){
        c_ter = true;
      }
      //else{
        //peptide_n = next.FindAtom("N");
        //peptide_c = i->FindAtom("C");
        //nucleotide_p = next.FindAtom("P");
        //nucleotide_o = i->FindAtom("O3");

        //if(!ost::mol::BondExists(peptide_n,peptide_c) && !ost::mol::BondExists(nucleotide_p,nucleotide_o)){
          //c_ter = true;
        //}
      //}
      if(c_ter && n_ter){
        ed.RenameResidue(*i,res_renaming_ff_specific_.find(i->GetName())->second->twoter);
        continue;
      }
      if(n_ter){
        ed.RenameResidue(*i,res_renaming_ff_specific_.find(i->GetName())->second->nter);
        continue;
      }
      if(c_ter){
        ed.RenameResidue(*i,res_renaming_ff_specific_.find(i->GetName())->second->cter);
        continue;
      }
      ed.RenameResidue(*i,res_renaming_ff_specific_.find(i->GetName())->second->main);
    }

    for(ost::mol::ResidueHandleList::iterator i = res_list.begin();i!=res_list.end();++i){
      boost::unordered_map<String,std::vector<std::pair<String,String> > >::const_iterator map_it = atom_renaming_ff_specific_.find(i->GetName());
      if(map_it!=atom_renaming_ff_specific_.end()){
        for(std::vector<std::pair<String,String> >::const_iterator j = map_it->second.begin();
            j!=map_it->second.end();++j){
          atom = i->FindAtom(j->first);
          if(atom.IsValid()) ed.RenameAtom(atom,j->second); 
        }
      }
    }  
    return;
  }

  //it seems, that we have to do the reverse ff renaming => ff specific to gromacs standard

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();i!=res_list.end();++i){
    boost::unordered_map<String,std::vector<std::pair<String,String> > >::const_iterator map_it = atom_renaming_ff_specific_.find(i->GetName());
    if(map_it!=atom_renaming_ff_specific_.end()){
      for(std::vector<std::pair<String,String> >::const_iterator j = map_it->second.begin();
          j!=map_it->second.end();++j){
        atom = i->FindAtom(j->second);
        if(atom.IsValid()) ed.RenameAtom(atom,j->first); 
      }
    }
  }  

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin(); i != res_list.end(); ++i){
    for(boost::unordered_map<String,ResidueNamesPtr>::const_iterator j = res_renaming_ff_specific_.begin();
        j != res_renaming_ff_specific_.end(); ++j){
      if(j->second->Contains(i->GetName())){
        ed.RenameResidue(*i,j->first);
        break;
      }
    }
  }
}

String Forcefield::AtomTypesToKeyword(std::vector<String>& types, bool allow_reordering) const{
  bool revert = false;
  int num_types = types.size();

  int comparisons = num_types/2; //integer division => truncate
  for(int i = 0; i < comparisons; ++i){
    //note, that in case of equality the next iteration gets entered
    if(types[i] < types[num_types - i - 1]) break;
    else if(types[i] > types[num_types - i - 1]){
      revert = true;
      break;
    }
  }

  if(revert && allow_reordering) std::reverse(types.begin(),types.end());
  std::stringstream ss;
  std::vector<String>::iterator i = types.begin();
  ss << *i++;
  for (; i != types.end(); ++i) {
      ss << "_" << *i;
  }
  return ss.str();
}

void Forcefield::CheckInteractionToAdd(InteractionPtr p, const String& interaction_type) const{
  bool valid = true;
  if(p->GetTypes().empty()) valid = false;
  if(!p->IsParametrized()) valid = false;
  if(!valid){
    throw ost::Error("Expect interaction to have types and parameters properly set!");
  }
  if(interaction_type == "BOND"){
    if(p->GetFuncType() != HarmonicBond){
      throw ost::Error("Invalid function type encountered when adding bond!");
    } return;
  }
  if(interaction_type == "ANGLE"){
    if(p->GetFuncType() != HarmonicAngle && p->GetFuncType() != UreyBradleyAngle){
      throw ost::Error("Invalid function type encountered when adding angle!");  
    } return;
  }
  if(interaction_type == "DIHEDRAL"){
    if(p->GetFuncType() != PeriodicDihedral){       
      throw ost::Error("Invalid function type encountered when adding dihedral!");  
    } return;
  }
  if(interaction_type == "IMPROPER"){
    if(p->GetFuncType() != PeriodicImproper && p->GetFuncType() != HarmonicImproper){       
      throw ost::Error("Invalid function type encountered when adding improper!"); 
    } return;
  }
  if(interaction_type == "CMAP"){
    if(p->GetFuncType() != CMap){       
      throw ost::Error("Invalid function type encountered when adding cmap!");  
    } return;
  }
  if(interaction_type == "IMPLICIT_GENBORN"){
    if(p->GetFuncType() != GBSA){       
      throw ost::Error("Invalid function type encountered when adding implicit genborn!");  
    } return;
  }
  if(interaction_type == "LJ"){
    if(p->GetFuncType() != LJ){       
      throw ost::Error("Invalid function type encountered when adding lj!");  
    } return;
  }
  if(interaction_type == "LJ_PAIR"){
    if(p->GetFuncType() != LJPair){       
      throw ost::Error("Invalid function type encountered when adding lj pair!");  
    } return;
  }
  if(interaction_type == "CONSTRAINT"){
    if(p->GetFuncType() != DistanceConstraint){       
      throw ost::Error("Invalid function type encountered when adding constraint!");  
    } return;
  }

  std::stringstream ss;
  ss << "Cannot validate added interaction of type "<< interaction_type << "!";  
  throw ost::Error(ss.str());

}

}}}//ns
