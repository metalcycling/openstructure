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

#include <ost/mol/mm/buildingblock.hh>

namespace ost{ namespace mol{ namespace mm{

BuildingBlock::BuildingBlock(const BuildingBlock& block){

  atoms_ = block.GetAtoms();
  types_ = block.GetTypes();
  charges_ = block.GetCharges();
  masses_ = block.GetMasses();

  //We have to force to create new interactions with new!!!
  //There would be to danger to mess around with an interaction,
  //that is associated with another building block.

  std::vector<InteractionPtr> bonds = block.GetBonds();
  std::vector<InteractionPtr> angles = block.GetAngles();
  std::vector<InteractionPtr> dihedrals = block.GetDihedrals();
  std::vector<InteractionPtr> impropers = block.GetImpropers();
  std::vector<InteractionPtr> exclusions = block.GetExclusions();
  std::vector<InteractionPtr> cmaps = block.GetCMaps();
  std::vector<InteractionPtr> constraints = block.GetConstraints();

  for(std::vector<InteractionPtr>::iterator i = bonds.begin();
      i != bonds.end(); ++i){
    bonds_.push_back(InteractionPtr(new Interaction(**i)));
  }
  for(std::vector<InteractionPtr>::iterator i = angles.begin();
      i != angles.end(); ++i){
    angles_.push_back(InteractionPtr(new Interaction(**i)));
  }
  for(std::vector<InteractionPtr>::iterator i = dihedrals.begin();
      i != dihedrals.end(); ++i){
    dihedrals_.push_back(InteractionPtr(new Interaction(**i)));
  }
  for(std::vector<InteractionPtr>::iterator i = impropers.begin();
      i != impropers.end(); ++i){
    impropers_.push_back(InteractionPtr(new Interaction(**i)));
  }
  for(std::vector<InteractionPtr>::iterator i = exclusions.begin();
      i != exclusions.end(); ++i){
    exclusions_.push_back(InteractionPtr(new Interaction(**i)));
  }
  for(std::vector<InteractionPtr>::iterator i = cmaps.begin();
      i != cmaps.end(); ++i){
    cmaps_.push_back(InteractionPtr(new Interaction(**i)));
  }
  for(std::vector<InteractionPtr>::iterator i = constraints.begin();
      i != constraints.end(); ++i){
    constraints_.push_back(InteractionPtr(new Interaction(**i)));
  }
}

int BuildingBlock::GetAtomIndex(const String& atom_name) const{

  uint ind = std::find(atoms_.begin(),atoms_.end(),atom_name) - atoms_.begin();
  
  if(ind<atoms_.size()){
    return ind;
  }
  return -1;
}

bool BuildingBlock::Match(const ost::mol::ResidueHandle& handle, bool match_connectivity, String& info) const{

  //check for size
  if(uint(handle.GetAtomCount()) != atoms_.size()){
    std::stringstream ss;
    ss << "Expected "<<atoms_.size()<<" atoms, got "<<handle.GetAtomCount();
    info = ss.str();
    return false;
  }

  //check for presence of particular atoms
  for(std::vector<String>::const_iterator i = atoms_.begin(); i!=atoms_.end(); ++i){
    if(!handle.FindAtom(*i).IsValid()){
      std::stringstream ss;
      ss << "Buildingblock has atom of name " << *i;
      ss << ", which is not present in residue";
      info = ss.str();
      return false;
    }
  }

  if(!match_connectivity) return true;

  //check connectivity by extracting all unique bonds from handle,
  //also bonds to other residues.
  std::set<std::pair<String,String> > raw_bonds_handle;
  std::vector<InteractionPtr> bonds_handle;
  ost::mol::AtomHandleList atom_list = handle.GetAtomList();
  ost::mol::ResidueHandle prev = handle.GetPrev();
  ost::mol::ResidueHandle next = handle.GetNext();
  ost::mol::AtomHandle atom_one;
  ost::mol::AtomHandle atom_two;
  String name_one;
  String name_two;
  std::vector<String> names;
  names.push_back("");
  names.push_back("");
  ost::mol::BondHandleList bond_list;

  for(ost::mol::AtomHandleList::iterator i = atom_list.begin();
      i != atom_list.end(); ++i){
    bond_list = i->GetBondList();
    for(ost::mol::BondHandleList::iterator j = bond_list.begin(); 
        j != bond_list.end(); ++j){
      atom_one = j->GetFirst();
      atom_two = j->GetSecond();
      name_one = atom_one.GetName();
      name_two = atom_two.GetName();

      if(atom_one.GetResidue() == prev) name_one = "-"+name_one;
      else if(atom_one.GetResidue() == next) name_one = "+"+name_one;
      else if(atom_one.GetResidue() != handle) continue; //in the building block
                                                         //there can only be atoms
                                                         //to previous or subsequent
                                                         //residues. So we don't check
                                                         //cases as they occur for example
                                                         //in disulfid bridges
      if(atom_two.GetResidue() == prev) name_two = "-"+name_two;
      else if(atom_two.GetResidue() == next) name_two = "+"+name_two;
      else if(atom_two.GetResidue() != handle) continue; //dito

      if(atom_one.GetIndex()<atom_two.GetIndex()) raw_bonds_handle.insert(std::make_pair(name_one,name_two));
      else raw_bonds_handle.insert(std::make_pair(name_two,name_one)); 

    }
  }

  for(std::set<std::pair<String,String> >::iterator i = raw_bonds_handle.begin();
      i != raw_bonds_handle.end(); ++i){
    InteractionPtr p(new Interaction(HarmonicBond));
    names[0] = i->first;
    names[1] = i->second;
    p->SetNames(names);
    bonds_handle.push_back(p);
  }

  bool found;
  //let's first check for bonds
  for(std::vector<InteractionPtr>::const_iterator i = bonds_.begin();
      i!=bonds_.end(); ++i){
    found = false;
    for(std::vector<InteractionPtr>::iterator j = bonds_handle.begin();
        j != bonds_handle.end(); ++j){
      if((*j)->MatchNames((*i)->GetNames())){
        found = true;
        std::vector<String> names = (*j)->GetNames();
        raw_bonds_handle.erase(std::make_pair(names[0],names[1]));
        break;
      }
    }
    if(!found){
      std::stringstream ss;
      std::vector<String> names_for_info = (*i)->GetNames();
      ss << "Buildingblock defines bond between \"" <<names_for_info[0];
      ss << "\" and \"" << names_for_info[1] << "\". Those atoms are not connected in your residue";
      info = ss.str();
      return false;
    }
  }
  //If there is a distance constraint matching a remaining raw bond, we also
  //remove it from raw bonds
  for(std::vector<InteractionPtr>::const_iterator i = constraints_.begin();
      i != constraints_.end(); ++i){
    for(std::vector<InteractionPtr>::iterator j = bonds_handle.begin();
        j != bonds_handle.end(); ++j){
      if((*j)->MatchNames((*i)->GetNames())){
        std::vector<String> names = (*j)->GetNames();
        raw_bonds_handle.erase(std::make_pair(names[0],names[1]));
        break;
      }
    }
  }
  //we finally check whether there is an additional internal bond in the
  //residue. Note, that a remaining bond to another residue (starting with +/-)
  //is not a criterium for failing the match procedure
  for(std::set<std::pair<String,String> >::iterator i = raw_bonds_handle.begin();
      i != raw_bonds_handle.end(); ++i){
    if( i->first.at(0) == '+' || i->first.at(0) == '-') continue;
    if( i->second.at(0) == '+' || i->second.at(0) == '-') continue;
    std::stringstream ss;
    ss << "Your residue has a bond between atoms of name \"" << i->first << "\" and \"";
    ss << i->second << "\", which is not defined in the buildingblock";
    info = ss.str(); 
    return false;
  }
  return true;
}

void BuildingBlock::Connect(ost::mol::ResidueHandle& handle, ost::mol::XCSEditor& ed){

  if(!handle.IsValid()) throw ost::Error("Invalid Residue encountered in connecting procedure!");

  ost::mol::ResidueHandle prev = handle.GetPrev();
  ost::mol::ResidueHandle next = handle.GetNext();

  //Let's build the connectivity
  ost::mol::AtomHandle atom1,atom2;
  std::vector<String> names;
  String name1, name2;


  for(std::vector<InteractionPtr>::iterator i = bonds_.begin();
      i != bonds_.end(); ++i){

    names = (*i)->GetNames();
    name1 = names[0];
    name2 = names[1];

    if(name1[0] == '-'){
      if(!prev.IsValid()) throw ost::Error("Cannot access previous residue in connecting procedure!");
      atom1 = prev.FindAtom(name1.substr(1));
    }
    else if(name1[0] == '+'){
      if(!next.IsValid()) throw ost::Error("Cannot access next residue in connecting procedure!");
      atom1 = next.FindAtom(name1.substr(1));
    }
    else{
      atom1 = handle.FindAtom(name1);
    }

    if(name2[0] == '-'){
      if(!prev.IsValid()) throw ost::Error("Cannot access previous residue in connecting procedure!");
      atom2 = prev.FindAtom(name2.substr(1));
    }
    else if(name2[0] == '+'){
      if(!next.IsValid()) throw ost::Error("Cannot access next residue in connecting procedure!");
      atom2 = next.FindAtom(name2.substr(1));
    }
    else{
      atom2 = handle.FindAtom(name2);
    }

    if(!atom1.IsValid() || !atom2.IsValid()) throw ost::Error("Could not find required atom in connecting procedure!");
   
    ed.Connect(atom1,atom2);
  }
}

String BuildingBlock::GetType(const String& name) const{
  int index = this->GetAtomIndex(name);
  if(index!=-1) return types_[index];
  std::stringstream ss;
  ss << "Could not find requested atom of name \""<<name;
  ss << "\" in buildingblock!";
  throw ost::Error(ss.str());
}

Real BuildingBlock::GetCharge(const String& name) const{
  int index = this->GetAtomIndex(name);
  if(index!=-1) return charges_[index];
  std::stringstream ss;
  ss << "Could not find requested atom of name \""<<name;
  ss << "\" in buildingblock!";
  throw ost::Error(ss.str());
}

Real BuildingBlock::GetMass(const String& name) const{
  int index = this->GetAtomIndex(name);
  if(index!=-1) return masses_[index];
  std::stringstream ss;
  ss << "Could not find requested atom of name \""<<name;
  ss << "\" in buildingblock!";
  throw ost::Error(ss.str());
}


void BuildingBlock::AddBond(InteractionPtr p, bool replace_existing){
  this->CheckInteractionToAdd(p);
  if(replace_existing){
    std::vector<String> names = p->GetNames();
    for(uint i = 0; i < bonds_.size(); ++i){
      if(bonds_[i]->MatchNames(names)){
        bonds_.erase(bonds_.begin() + i);
        --i;
      }
    }
  }
  bonds_.push_back(p);
}

void BuildingBlock::AddAngle(InteractionPtr p, bool replace_existing){
  this->CheckInteractionToAdd(p);
  if(replace_existing){
    std::vector<String> names = p->GetNames();
    for(uint i = 0; i < angles_.size(); ++i){
      if(angles_[i]->MatchNames(names)){
        angles_.erase(angles_.begin() + i);
        --i;
      }
    }
  }
  angles_.push_back(p);
}

void BuildingBlock::AddDihedral(InteractionPtr p, bool replace_existing){
  this->CheckInteractionToAdd(p);
  if(replace_existing){
    std::vector<String> names = p->GetNames();
    for(uint i = 0; i < dihedrals_.size(); ++i){
      if(dihedrals_[i]->MatchNames(names)){
        dihedrals_.erase(dihedrals_.begin() + i);
        --i;
      }
    }
  }
  dihedrals_.push_back(p);
}

void BuildingBlock::AddImproper(InteractionPtr p, bool replace_existing){
  this->CheckInteractionToAdd(p);
  if(replace_existing){
    std::vector<String> names = p->GetNames();
    for(uint i = 0; i < impropers_.size(); ++i){
      if(impropers_[i]->MatchNames(names)){
        impropers_.erase(impropers_.begin() + i);
        --i;
      }
    }
  }
  impropers_.push_back(p);
}

void BuildingBlock::AddExclusion(InteractionPtr p, bool replace_existing){
  this->CheckInteractionToAdd(p);
  if(replace_existing){
    std::vector<String> names = p->GetNames();
    for(uint i = 0; i < exclusions_.size(); ++i){
      if(exclusions_[i]->MatchNames(names)){
        exclusions_.erase(exclusions_.begin() + i);
        --i;
      }
    }
  }
  exclusions_.push_back(p);
}

void BuildingBlock::AddCMap(InteractionPtr p, bool replace_existing){
  this->CheckInteractionToAdd(p);
  if(replace_existing){
    std::vector<String> names = p->GetNames();
    for(uint i = 0; i < cmaps_.size(); ++i){
      if(cmaps_[i]->MatchNames(names)){
        cmaps_.erase(cmaps_.begin() + i);
        --i;
      }
    }
  }
  cmaps_.push_back(p);
}

void BuildingBlock::AddConstraint(InteractionPtr p, bool replace_existing){
  this->CheckInteractionToAdd(p);
  if(replace_existing){
    std::vector<String> names = p->GetNames();
    for(uint i = 0; i < constraints_.size(); ++i){
      if(constraints_[i]->MatchNames(names)){
        constraints_.erase(constraints_.begin() + i);
        --i;
      }
    }
  }
  constraints_.push_back(p);
}

void BuildingBlock::RemoveAtom(const String& name){
  int index = this->GetAtomIndex(name);
  if(index != -1){
    atoms_.erase(atoms_.begin()+index);
    types_.erase(types_.begin()+index);
    charges_.erase(charges_.begin()+index);
    masses_.erase(masses_.begin()+index);
  }

  //All interactions, that contain this atom get deleted as well
  //note, that we're looking at all interactions independently whether tere
  //really is such an atom in the building block
  for(uint i = 0; i < bonds_.size(); ++i){
    if(bonds_[i]->HasName(name)){
      bonds_.erase(bonds_.begin()+i);
      --i;
    }
  }
  for(uint i = 0; i < angles_.size(); ++i){
    if(angles_[i]->HasName(name)){
      angles_.erase(angles_.begin()+i);
      --i;
    }
  }
  for(uint i = 0; i < dihedrals_.size(); ++i){
    if(dihedrals_[i]->HasName(name)){
      dihedrals_.erase(dihedrals_.begin()+i);
      --i;
    }
  }
  for(uint i = 0; i < impropers_.size(); ++i){
    if(impropers_[i]->HasName(name)){
      impropers_.erase(impropers_.begin()+i);
      --i;
    }
  }
  for(uint i = 0; i < exclusions_.size(); ++i){
    if(exclusions_[i]->HasName(name)){
      exclusions_.erase(exclusions_.begin()+i);
      --i;
    }
  }
  for(uint i = 0; i < cmaps_.size(); ++i){
    if(cmaps_[i]->HasName(name)){
      cmaps_.erase(cmaps_.begin()+i);
      --i;
    }
  }
  for(uint i = 0; i < constraints_.size(); ++i){
    if(constraints_[i]->HasName(name)){
      constraints_.erase(constraints_.begin()+i);
      --i;
    }
  }
}

void BuildingBlock::ReplaceAtom(const String& name,const String& new_name,
                                const String& new_type, Real new_charge, Real new_mass){

  int index = this->GetAtomIndex(name);
  if(index != -1){
    atoms_[index] = new_name;
    types_[index] = new_type;
    charges_[index] = new_charge;
    masses_[index] = new_mass;

    //we don't want to search every interaction associated with this
    //atom, so we call the ReplaceAtom method for every interaction
    //nothing will happen, if the atom of interest is not part of
    //that interaction.

    for(std::vector<InteractionPtr>::iterator i = bonds_.begin();
        i!=bonds_.end(); ++i){
      (*i)->ReplaceAtom(name, new_name, new_type);
    }

    for(std::vector<InteractionPtr>::iterator i = angles_.begin();
        i!=angles_.end(); ++i){
      (*i)->ReplaceAtom(name, new_name, new_type);
    }

    for(std::vector<InteractionPtr>::iterator i = dihedrals_.begin();
        i!=dihedrals_.end(); ++i){
      (*i)->ReplaceAtom(name, new_name, new_type);
    }

    for(std::vector<InteractionPtr>::iterator i = impropers_.begin();
        i!=impropers_.end(); ++i){
      (*i)->ReplaceAtom(name, new_name, new_type);
    }

    for(std::vector<InteractionPtr>::iterator i = exclusions_.begin();
        i!=exclusions_.end(); ++i){
      (*i)->ReplaceAtom(name, new_name, new_type);
    }

    for(std::vector<InteractionPtr>::iterator i = cmaps_.begin();
        i!=cmaps_.end(); ++i){
      (*i)->ReplaceAtom(name, new_name, new_type);
    }
  }
}

void BuildingBlock::AddAtom(const String& name, const String& type, Real charge, Real mass){
  if (std::find(atoms_.begin(),atoms_.end(),name) != atoms_.end()) return;
  atoms_.push_back(name);
  types_.push_back(type);
  charges_.push_back(charge);
  masses_.push_back(mass);
}

void BuildingBlock::RemoveInteractionsToPrev(){
  std::vector<String> names;
  bool has_prev_atom;

  for(uint i = 0; i < bonds_.size(); ++i){
    has_prev_atom = false;
    names = bonds_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '-'){
        has_prev_atom = true;
        break;
      }
    }
    if(has_prev_atom){
      bonds_.erase(bonds_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < angles_.size(); ++i){
    has_prev_atom = false;
    names = angles_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '-'){
        has_prev_atom = true;
        break;
      }
    }
    if(has_prev_atom){
      angles_.erase(angles_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < dihedrals_.size(); ++i){
    has_prev_atom = false;
    names = dihedrals_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '-'){
        has_prev_atom = true;
        break;
      }
    }
    if(has_prev_atom){
      dihedrals_.erase(dihedrals_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < impropers_.size(); ++i){
    has_prev_atom = false;
    names = impropers_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '-'){
        has_prev_atom = true;
        break;
      }
    }
    if(has_prev_atom){
      impropers_.erase(impropers_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < exclusions_.size(); ++i){
    has_prev_atom = false;
    names = exclusions_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '-'){
        has_prev_atom = true;
        break;
      }
    }
    if(has_prev_atom){
      exclusions_.erase(exclusions_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < cmaps_.size(); ++i){
    has_prev_atom = false;
    names = cmaps_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '-'){
        has_prev_atom = true;
        break;
      }
    }
    if(has_prev_atom){
      cmaps_.erase(cmaps_.begin()+i);
      --i;
    }
  }
}

void BuildingBlock::RemoveInteractionsToNext(){
  std::vector<String> names;
  bool has_next_atom;

  for(uint i = 0; i < bonds_.size(); ++i){
    has_next_atom = false;
    names = bonds_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '+'){
        has_next_atom = true;
        break;
      }
    }
    if(has_next_atom){
      bonds_.erase(bonds_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < angles_.size(); ++i){
    has_next_atom = false;
    names = angles_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '+'){
        has_next_atom = true;
        break;
      }
    }
    if(has_next_atom){
      angles_.erase(angles_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < dihedrals_.size(); ++i){
    has_next_atom = false;
    names = dihedrals_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '+'){
        has_next_atom = true;
        break;
      }
    }
    if(has_next_atom){
      dihedrals_.erase(dihedrals_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < impropers_.size(); ++i){
    has_next_atom = false;
    names = impropers_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '+'){
        has_next_atom = true;
        break;
      }
    }
    if(has_next_atom){
      impropers_.erase(impropers_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < exclusions_.size(); ++i){
    has_next_atom = false;
    names = exclusions_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '+'){
        has_next_atom = true;
        break;
      }
    }
    if(has_next_atom){
      exclusions_.erase(exclusions_.begin()+i);
      --i;
    }
  }

  for(uint i = 0; i < cmaps_.size(); ++i){
    has_next_atom = false;
    names = cmaps_[i]->GetNames();
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      if((*j)[0] == '+'){
        has_next_atom = true;
        break;
      }
    }
    if(has_next_atom){
      cmaps_.erase(cmaps_.begin()+i);
      --i;
    }
  }
}

void BuildingBlock::CheckInteractionToAdd(InteractionPtr p) const{
  if(p->GetNames().empty()){
    throw ost::Error("Expect interaction to have names properly set!");
  }
}

}}}//ns
