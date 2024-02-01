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

#include <ost/mol/mm/gromacs_block_modifiers.hh>

namespace ost{ namespace mol{ namespace mm{

std::vector<geom::Vec3> GromacsPositionRuleEvaluator::EvaluatePosRule(
                                      int rule, int number,
                                      const std::vector<geom::Vec3>& anchors) {

  std::vector<geom::Vec3> positions;

  switch(rule){

    case 1:{
      geom::Vec3 i_j = geom::Normalize(anchors[0]-anchors[1]);
      geom::Vec3 i_k = geom::Normalize(anchors[0]-anchors[2]);
      geom::Vec3 temp = geom::Normalize(i_j+i_k);
      positions.push_back(anchors[0]+temp);
      break;
    }

    case 2:{
      Real x = 0.333806859234;
      Real y = 0.942641491092;

      geom::Vec3 i = anchors[0];
      geom::Vec3 j = anchors[1];
      geom::Vec3 k = anchors[2];

      geom::Vec3 j_i = geom::Normalize(i-j);
      geom::Vec3 q = geom::Normalize(geom::Cross(geom::Cross(k-j,j_i),j_i));

      positions.push_back(i+x*j_i+y*q);
      break;
    }

    case 3:{
      Real x = 0.5;
      Real y = 0.866025403784;

      geom::Vec3 i = anchors[0];
      geom::Vec3 j = anchors[1];
      geom::Vec3 k = anchors[2];
      
      geom::Vec3 j_i = geom::Normalize(i-j);
      geom::Vec3 q = geom::Normalize(geom::Cross(geom::Cross(k-j,j_i),j_i));

      positions.push_back(i+x*j_i-y*q);
      positions.push_back(i+x*j_i+y*q);
      break;
    }

    case 4:{
      Real x = 0.333313247568;
      Real y = 0.942816142732;
      Real z = 0.866025403784;


      geom::Vec3 i = anchors[0];
      geom::Vec3 j = anchors[1];
      geom::Vec3 k = anchors[2];
      
      geom::Vec3 j_i = geom::Normalize(i-j);
      geom::Vec3 q = geom::Normalize(geom::Cross(geom::Cross(k-j,j_i),j_i));

      geom::Vec3 n1 = i+x*j_i+y*q;

      geom::Vec3 temp = geom::Normalize(geom::Cross(j_i,n1-i));

      geom::Vec3 n2 = i+x*j_i-0.5*q+z*temp;
      geom::Vec3 n3 = i+x*j_i-0.5*q-z*temp;

      positions.push_back(n1);
      positions.push_back(n2);
      if(number==3) positions.push_back(n3);
      break;
    }

    case 5:{
      geom::Vec3 i = anchors[0];
      geom::Vec3 i_j = geom::Normalize(i-anchors[1]);
      geom::Vec3 i_k = geom::Normalize(i-anchors[2]);
      geom::Vec3 i_l = geom::Normalize(i-anchors[3]);
      geom::Vec3 temp = geom::Normalize(i_j+i_k+i_l);
      positions.push_back(i+temp);
      break;
    }

    case 6:{
      Real x = 0.81649043092;
      Real y = 0.577358966516;

      geom::Vec3 i = anchors[0];
      geom::Vec3 i_j = geom::Normalize(i-anchors[1]);
      geom::Vec3 i_k = geom::Normalize(i-anchors[2]);
      geom::Vec3 temp_one = geom::Normalize(i_j+i_k);
      geom::Vec3 temp_two = geom::Normalize(geom::Cross(i_j,i_k));
      positions.push_back(i+y*temp_one+x*temp_two);
      positions.push_back(i+y*temp_one-x*temp_two);
      break;
    }

    case 7:{
      //not implemented yet, do nothing
      break;
    }

    case 8:{
      Real x = 1.36*0.522498564716;
      Real y = 1.36*0.852640164354;

      geom::Vec3 i = anchors[0];
      geom::Vec3 j = anchors[1];
      geom::Vec3 k = anchors[2];
      
      geom::Vec3 j_i = geom::Normalize(i-j);
      geom::Vec3 q = geom::Normalize(geom::Cross(geom::Cross(k-j,j_i),j_i));

      positions.push_back(i+x*j_i-y*q);
      positions.push_back(i+x*j_i+y*q);
      break;
    }

    case 9:{
      Real x1 = 1.23*0.51503807491;
      Real x2 = 1.25*0.422618261741;
      Real y1 = 1.23*0.857167300702;
      Real y2 = 1.25*0.906307787037;

      geom::Vec3 i = anchors[0];
      geom::Vec3 j = anchors[1];
      geom::Vec3 k = anchors[2];
      
      geom::Vec3 j_i = geom::Normalize(i-j);
      geom::Vec3 q = geom::Normalize(geom::Cross(geom::Cross(k-j,j_i),j_i));

      positions.push_back(i+x1*j_i-y1*q);
      positions.push_back(i+x2*j_i+y2*q);

      //we have to construct a third position hanging at the second one...

      Real x = 0.333806859234;
      Real y = 0.942641491092;

      i = positions[1];
      j = anchors[0];
      k = anchors[1];

      j_i = geom::Normalize(i-j);
      q = geom::Normalize(geom::Cross(geom::Cross(k-j,j_i),j_i));

      positions.push_back(i+x*j_i+y*q);
      break;
    }

    case 10:{
      //not implemented yet, do nothing
      break;
    }

    case 11:{
      //not implemented yet, do nothing
      break;
    }

    default:{
      break;
    }
  }
  return positions;
}


void GromacsHydrogenConstructor::ApplyOnBuildingBlock(BuildingBlockPtr p) {
  //we assume, that the hydrogens are already correct in the building block...
}

void GromacsHydrogenConstructor::ApplyOnResidue(ost::mol::ResidueHandle& res,
                                                ost::mol::XCSEditor& ed) {

  if(!res.IsValid()) return;

  int number;
  int method;
  String prefix;

  std::vector<String> anchor_names;
  ost::mol::AtomHandleList anchor_atoms;

  std::vector<String> hydrogen_names;
  std::vector<geom::Vec3> hydrogen_positions;
  std::vector<geom::Vec3> anchor_positions;

  for(uint a=0;a<add_number_.size();++a){

    hydrogen_names = hydrogen_names_[a];

    number = add_number_[a];
    method = methods_[a];
    anchor_names = anchor_atom_names_[a];

    hydrogen_positions.clear();
    anchor_positions.clear();
    anchor_atoms.clear();

    String atom_name;
    ost::mol::ResidueHandle temp_res;
    ost::mol::AtomHandle atom;
    bool found_anchors = true;

    for (std::vector<String>::iterator it = anchor_names.begin();
         it!=anchor_names.end(); ++it) {

      temp_res = res;
      atom_name = *it;

      if(atom_name[0] == '-'){
        atom_name.erase(0,1);
        temp_res = res.GetPrev();
        if(!temp_res.IsValid()){
          //happens in case of termini... In this case the adding definitions
          //are taken from the termini database...
          found_anchors = false;
          break;
        }
      }
      if(atom_name[0] == '+'){
        atom_name.erase(0,1);
        temp_res = res.GetNext();
        if(!temp_res.IsValid()){
          //happens in case of termini... In this case the adding definitions
          //are taken from the termini database...
          found_anchors = false;
          break;
        }
      }

      atom = temp_res.FindAtom(atom_name);
      if(!atom.IsValid()){
        std::stringstream ss;
        ss << "cannot find required anchor atom " << atom_name << " in residue ";
        ss << res.GetQualifiedName() << " to add hydrogen!";
        throw ost::Error(ss.str());
      }
      
      anchor_atoms.push_back(atom);
      anchor_positions.push_back(atom.GetPos());
    }
    
    if(!found_anchors){
      continue;
    }
    
    hydrogen_positions
     = GromacsPositionRuleEvaluator::EvaluatePosRule(method, number,
                                                     anchor_positions);

    for(int b=0;b<number;++b){
      //only add hydrogen if not already present!
      atom = res.FindAtom(hydrogen_names[b]);
      if (!atom.IsValid()) {
        atom = ed.InsertAtom(res, hydrogen_names[b], hydrogen_positions[b], "H");
      }
      else {
        ed.SetAtomPos(atom, hydrogen_positions[b]);
      }
    }
  }
}


void GromacsHydrogenConstructor::AddHydrogenRule(
                                uint number, int method,
                                const std::vector<String>& hydrogen_names,
                                const std::vector<String>& anchors) {

  if(number != hydrogen_names.size()){
    std::stringstream ss;
    ss << "Tried to add hydrogen rule to hydrogen constructor with add number ";
    ss << number << ", but "<< hydrogen_names.size() << "hydrogen names";
  }

  hydrogen_names_.push_back(hydrogen_names);
  add_number_.push_back(number);
  methods_.push_back(method);
  anchor_atom_names_.push_back(anchors);
}

void GromacsBlockModifier::AddReplaceRule(const String& name, const String& new_name,
                                               const String& new_type, Real new_charge){
  replace_old_atom_name_.push_back(name);
  replace_new_atom_name_.push_back(new_name);
  replace_new_atom_type_.push_back(new_type);
  replace_new_charge_.push_back(new_charge);
}

void GromacsBlockModifier::AddAddRule(int number, int method, 
                                           const std::vector<String>& atom_names,
                                           const std::vector<String>& anchors,
                                           const String& type, Real charge){
  add_add_number_.push_back(number);
  add_methods_.push_back(method);
  add_atom_names_.push_back(atom_names);
  add_anchor_atom_names_.push_back(anchors);
  add_atom_types_.push_back(type);
  add_charges_.push_back(charge);
}

void GromacsBlockModifier::ApplyOnBuildingBlock(BuildingBlockPtr p){

  //let's resolve the replace statement
  for(uint i = 0; i < replace_old_atom_name_.size(); ++i){
    p->ReplaceAtom(replace_old_atom_name_[i], replace_new_atom_name_[i],
                   replace_new_atom_type_[i], replace_new_charge_[i]);

  }
  //let's resolve the add statement
  for(uint i = 0; i < add_add_number_.size(); ++i){
    String type = add_atom_types_[i];
    Real charge = add_charges_[i];
    std::vector<String> anchor_atoms = add_anchor_atom_names_[i];
    std::vector<String> names = add_atom_names_[i];

    //let's add the stuff to the buildingblock
    for(std::vector<String>::iterator j = names.begin();
        j != names.end(); ++j){
      p->AddAtom(*j,type,charge);
      //the bond, that connects the newly added atom to its anchor is not necessarily
      //defined. Let's add it here. Note, that nothing happens, if this specific bond
      //is already part of the building block.
      InteractionPtr b_p(new Interaction(p->GetBonds()[0]->GetFuncType()));
      std::vector<String> bond_names;
      std::vector<String> bond_types;
      bond_names.push_back(anchor_atoms[0]);
      bond_names.push_back(*j);
      bond_types.push_back(p->GetType(anchor_atoms[0]));
      bond_types.push_back(type);
      b_p->SetNames(bond_names);
      b_p->SetTypes(bond_types);
      p->AddBond(b_p);
    }
  }

  //let's resolve the delete statement
  for(std::vector<String>::iterator i = delete_atom_names_.begin();
      i != delete_atom_names_.end(); ++i){
    p->RemoveAtom(*i);
  }
 
  //we simply add the new interactions to the building block
  for(std::vector<InteractionPtr>::iterator i = bonds_.begin();
      i != bonds_.end(); ++i){
    p->AddBond(*i, true);
  }

  for(std::vector<InteractionPtr>::iterator i = angles_.begin();
      i != angles_.end(); ++i){
    p->AddAngle(*i, true);
  }

  for(std::vector<InteractionPtr>::iterator i = dihedrals_.begin();
      i != dihedrals_.end(); ++i){
    p->AddDihedral(*i, true);
  }

  for(std::vector<InteractionPtr>::iterator i = impropers_.begin();
      i != impropers_.end(); ++i){
    p->AddImproper(*i, true);
  }

  for(std::vector<InteractionPtr>::iterator i = cmaps_.begin();
      i != cmaps_.end(); ++i){
    p->AddCMap(*i, true);
  } 
}

void GromacsBlockModifier::ApplyOnResidue(ost::mol::ResidueHandle& res,
                                          ost::mol::XCSEditor& ed) {

  //let's resolve the replace statement
  for(uint i = 0; i < replace_old_atom_name_.size(); ++i){
    ost::mol::AtomHandle at = res.FindAtom(replace_old_atom_name_[i]);
    if(at.IsValid()){
      ed.RenameAtom(at,replace_new_atom_name_[i]);
    }
  }


  //let's resolve the add statement
  for(uint i = 0; i < add_add_number_.size(); ++i){
    int number = add_add_number_[i];
    int method = add_methods_[i];
    std::vector<String> anchor_names = add_anchor_atom_names_[i];
    //String type = add_atom_types_[i];
    std::vector<geom::Vec3> anchor_positions;
    ost::mol::AtomHandleList anchor_atoms;
    std::vector<String> names = add_atom_names_[i];

    for(std::vector<String>::iterator j = anchor_names.begin();
        j != anchor_names.end(); ++j){
      ost::mol::AtomHandle at = res.FindAtom(*j);
      if(!at.IsValid()){
        std::stringstream ss;
        ss << "Cannot find anchor atom required to construct termini";
        ss << " in provided residue!";
        throw ost::Error(ss.str());
      }
      anchor_positions.push_back(at.GetPos());
      anchor_atoms.push_back(at);
    }
    std::vector<geom::Vec3> positions
     = GromacsPositionRuleEvaluator::EvaluatePosRule(method, number,
                                                     anchor_positions);

    String ele = "H";
    if(method == 8 || method == 9){
      ele = "O";
    }
    for(int j=0; j<number; ++j){
      ost::mol::AtomHandle at = res.FindAtom(names[j]);
      if(!at.IsValid()) at = ed.InsertAtom(res,names[j],positions[j],ele);
      else ed.SetAtomPos(at,positions[j]);
      ed.Connect(at, anchor_atoms[0]);
    }
  }
  //let's resolve the delete statement
  for(std::vector<String>::iterator i = delete_atom_names_.begin();
      i != delete_atom_names_.end(); ++i){
    ost::mol::AtomHandle at = res.FindAtom(*i);
    if(at.IsValid()){
      ed.DeleteAtom(at);
    }
  }
}


void GromacsBlockModifier::CheckInteractionToAdd(InteractionPtr p, const String& interaction_type) const{
  if(p->GetNames().empty()) throw ost::Error("Expect interaction to have names properly set!");

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

  std::stringstream ss;
  ss << "Cannot validate added interaction of type "<< interaction_type << "!";  
  throw ost::Error(ss.str());

}

}}}//ns
