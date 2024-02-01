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

#include <ost/mol/mm/interaction.hh>

namespace ost { namespace mol{ namespace mm{

Interaction::Interaction(FuncType func_type): func_type_(func_type), 
                                              set_parameters_(false), 
                                              has_type_wildcard_(false),
                                              has_name_wildcard_(false){ }


void Interaction::SetTypes(std::vector<String> types){
  if(!this->CheckSetNamesTypes(types)){
    throw ost::Error("Tried to set invalid number of types to interaction!");
  }
  atom_types_ = types;
  if(std::find(atom_types_.begin(),atom_types_.end(),"X") != atom_types_.end()){
    has_type_wildcard_ = true;
  }
  else has_type_wildcard_ = false;
}

void Interaction::SetNames(std::vector<String> names){
  if(!this->CheckSetNamesTypes(names)){
    throw ost::Error("Tried to set invalid number of names to interaction!");
  }
  atom_names_ = names;
  if(std::find(atom_names_.begin(),atom_names_.end(),"X") != atom_names_.end()){
    has_name_wildcard_ = true;
  }
  else has_name_wildcard_ = false;
}

void Interaction::SetParam(std::vector<Real>& parameters){
  if(!this->CheckSetParam(parameters)){
    throw ost::Error("Tried to set invalid number of parameters to interaction!");
  }  
  parameters_ = parameters;
  set_parameters_ = true;
}

ost::mol::AtomHandleList Interaction::GetAtoms(const ost::mol::ResidueHandle& res) const{

  ost::mol::AtomHandleList return_list;
  ost::mol::ResidueHandle prev = res.GetPrev();
  ost::mol::ResidueHandle next = res.GetNext();
  String atom_name;
  ost::mol::AtomHandle atom;

  for(std::vector<String>::const_iterator i = atom_names_.begin();
      i != atom_names_.end(); ++i){
    atom_name = (*i);
    if(atom_name[0] == '+'){
      atom_name = atom_name.substr(1);
      if(!next.IsValid()){
        throw ost::Error("Access to nonexistent next residue requested!");
      }
      atom = next.FindAtom(atom_name);
      if(atom.IsValid()){
        return_list.push_back(atom);
        continue;
      }
      throw ost::Error("Could not find required atom "+*i+"!");
    }
    if(atom_name[0] == '-'){
      atom_name = atom_name.substr(1);
      if(!prev.IsValid()){
        throw ost::Error("Access to nonexistent previous residue requested!");
      }
      atom = prev.FindAtom(atom_name);
      if(atom.IsValid()){
        return_list.push_back(atom);
        continue;
      }
      throw ost::Error("Could not find required atom "+*i+"!");
    }
    atom = res.FindAtom(atom_name);
    if(atom.IsValid()){
      return_list.push_back(atom);
      continue;
    }
    throw ost::Error("Could not find required atom "+*i+"!");
  }

  return return_list;
}

bool Interaction::MatchTypes(const std::vector<String>& atom_types) const {

  if(atom_types_.size() == 0) return false;
  if(atom_types.size() != atom_types_.size()) return false;

  bool match = true;
  bool reverse_match = true;

  uint i,j;

  for(i = 0, j = atom_types_.size()-1; i < atom_types_.size(); ++i,--j){
    if((atom_types[i] != atom_types_[i]) && (atom_types_[i] != "X")) match = false;
    if((atom_types[i] != atom_types_[j]) && (atom_types_[j] != "X")) reverse_match = false;
  }

  return match || reverse_match;
}

bool Interaction::MatchNames(const std::vector<String>& atom_names) const {

  if(atom_names_.size() == 0) return false;
  if(atom_names.size() != atom_names_.size()) return false;

  bool match = true;
  bool r_match = true;

  uint i,j;

  for(i = 0, j = atom_names_.size()-1; i < atom_names_.size(); ++i,--j){
    if((atom_names[i] != atom_names_[i]) && (atom_names_[i] != "X")) match = false;
    if((atom_names[i] != atom_names_[j]) && (atom_names_[j] != "X")) r_match = false;
  }

  return match || r_match;
}

bool Interaction::ReplaceAtom(const String& name, const String& new_name, const String& new_type){

  for(uint i = 0; i < atom_names_.size(); ++i){
    if(atom_names_[i] == name){
      atom_names_[i] = new_name;
      if(i < atom_types_.size()) atom_types_[i] = new_type;
      //it could be, that wildcard stuff has changed
      has_type_wildcard_ = std::find(atom_types_.begin(),atom_types_.end(),"X") != atom_types_.end();
      has_name_wildcard_ = std::find(atom_names_.begin(),atom_names_.end(),"X") != atom_names_.end();
      return true;
    }
  }
  return false;
}

bool Interaction::HasName(const String& name) const{
  return std::find(atom_names_.begin(),atom_names_.end(),name) != atom_names_.end();
}

bool Interaction::HasType(const String& type) const{
  return std::find(atom_types_.begin(),atom_types_.end(),type) != atom_types_.end();
}

bool Interaction::CheckSetNamesTypes(std::vector<String>& types){

  switch(func_type_){
    case HarmonicBond: return types.size() == 2;
    case HarmonicAngle: return types.size() == 3;
    case UreyBradleyAngle: return types.size() == 3;
    case PeriodicDihedral: return types.size() == 4;
    case PeriodicImproper: return types.size() == 4;
    case HarmonicImproper: return types.size() == 4;
    case CMap: return types.size() == 5;
    case LJ: return types.size() == 1;
    case LJPair: return types.size() == 2;
    case GBSA: return types.size() == 1;
    case DistanceConstraint: return types.size() == 2;
    case Exclusion: return types.size() == 2;
    case HarmonicPositionRestraint: return types.size() == 1;
    case HarmonicDistanceRestraint: return types.size() == 2;
    default: return false;
  }
}

bool Interaction::CheckSetParam(std::vector<Real>& param){

  switch(func_type_){
    case HarmonicBond: return param.size() == 2;
    case HarmonicAngle: return param.size() == 2;
    case UreyBradleyAngle: return param.size() == 4;
    case PeriodicDihedral: return param.size() == 3;
    case PeriodicImproper: return param.size() == 3;
    case HarmonicImproper: return param.size() == 2;
    case CMap:{
      int num_param = param.size();
      if(num_param < 1) return false;
      int x = param[0];
      return num_param == (x*x + 1);
    }
    case LJ: return param.size() == 2;
    case LJPair: return param.size() == 2;
    case GBSA: return param.size() == 2;
    case DistanceConstraint: return param.size() == 1;
    case Exclusion: return param.empty();
    case HarmonicPositionRestraint: return param.size() == 7;
    case HarmonicDistanceRestraint: return param.size() == 2;
    default: return false;
  }
  
}

}}} //ns