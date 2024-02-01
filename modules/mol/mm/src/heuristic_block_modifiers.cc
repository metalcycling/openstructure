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

#include <ost/mol/mm/heuristic_block_modifiers.hh>
#include <ost/geom/geom.hh>
#include <math.h>

namespace ost{ namespace mol{ namespace mm{


void HeuristicHydrogenConstructor::ApplyOnBuildingBlock(BuildingBlockPtr p){
  //we assume, that the hydrogens are already correct in the building block...
}

HeuristicHydrogenConstructor::HeuristicHydrogenConstructor(BuildingBlockPtr block){
  std::vector<String> atom_names = block->GetAtoms();
  std::vector<InteractionPtr> bonds = block->GetBonds();

  for(std::vector<String>::iterator aname = atom_names.begin(); 
      aname != atom_names.end(); ++aname){
    if((*aname)[0] == 'H') continue; //it's a hydrogen

    std::vector<String> antecedent_names, hydrogen_names;
    String aname_other;

    for (std::vector<InteractionPtr>::iterator i = bonds.begin() ; 
         i != bonds.end(); ++i){

      std::vector<String> bond_names=(*i)->GetNames();

      if (*aname!=bond_names[0] && *aname!=bond_names[1]) continue;//bond has nothing to do
                                                                   //with current atom
      if (*aname==bond_names[0]) aname_other=bond_names[1];
      else aname_other=bond_names[0]; 

      //check whether it's a hydrogen
      if (aname_other[0] == 'H') {
        hydrogen_names.push_back(aname_other);
      }
      else{
        antecedent_names.push_back(aname_other);
      }
    }
    if(!hydrogen_names.empty()){
      if(antecedent_names.empty()){
        std::stringstream ss;
        ss << "Can only reconstruct hydrogens when at least one ";
        ss << "antecedent atom is present!";
        throw ost::Error(ss.str());
      }
      if(antecedent_names.size() + hydrogen_names.size() > 4){
        std::stringstream ss;
        ss << "Antecedent and hydrogen atoms sum up to more than four ";
        ss << "around anchor atom "<<*aname<<". This is not possible with ";
        ss << "the used block modifier implementation!";
        throw ost::Error(ss.str());
      }
      anchor_atom_names_.push_back(*aname);
      antecedent_names_.push_back(antecedent_names);
      hydrogen_names_.push_back(hydrogen_names);
    }
  } 
}

void HeuristicHydrogenConstructor::ApplyOnResidue(ost::mol::ResidueHandle& res, ost::mol::XCSEditor& edi){

  if(!res.IsValid()) return;

  for(uint i = 0; i < anchor_atom_names_.size(); ++i){

    ost::mol::AtomHandle anchor_atom = res.FindAtom(anchor_atom_names_[i]);
    if(!anchor_atom.IsValid()){
      std::stringstream ss;
      ss << "Anchor atom of name " << anchor_atom_names_[i];
      ss << " cannot be found in residue "<<res.GetQualifiedName()<<"!";
      throw ost::Error(ss.str());
    }

    geom::Vec3List antecedent_positions;
    for(std::vector<String>::iterator j = antecedent_names_[i].begin();
        j != antecedent_names_[i].end(); ++j){
      ost::mol::AtomHandle at = res.FindAtom(*j);
      if(!at.IsValid()){
        std::stringstream ss;
        ss << "Cannot find required antecedent atom of name "<<*j;
        ss << " in residue "<<res.GetQualifiedName()<<"!";
        throw ost::Error(ss.str());
      }
      antecedent_positions.push_back(at.GetPos());
    }

    Real bond_length=1.08;
    String ele(1,'H');
    geom::Vec3 antecedents_cm=antecedent_positions.GetCenter();
    geom::Vec3 pa=anchor_atom.GetPos();
    geom::Vec3 v=geom::Normalize(pa-antecedents_cm);
    geom::Vec3 n=geom::OrthogonalVector(v),direction=geom::Vec3(),hpos=geom::Vec3();
    geom::Mat3 R=geom::Mat3();
    Real bond_angle;
    AtomHandle ha=mol::AtomHandle();

    switch (hydrogen_names_[i].size()){
      case 1:{  
        if (antecedent_names_[i].size()==1){
          bond_angle=0.4*M_PI;
        }
        else {
          bond_angle=0.0;
        }
        R=geom::AxisRotation(n,bond_angle);
        direction=R*v;
        hpos=pa+bond_length*direction;
        ha = res.FindAtom(hydrogen_names_[i][0]);
        if(!ha.IsValid()) ha=edi.InsertAtom(res,hydrogen_names_[i][0],hpos,ele);
        else edi.SetAtomPos(ha,hpos);
        break;
      }
      case 2:{
        if (antecedent_names_[i].size()==1){
          bond_angle=M_PI/3.;
        }
        else {
          bond_angle = 0.4*M_PI;
          n=geom::Normalize(antecedent_positions[0]-antecedent_positions[1]);
        }

        R=geom::AxisRotation(n,bond_angle);
        direction=R*v;
        hpos=pa+bond_length*direction;
        ha = res.FindAtom(hydrogen_names_[i][0]);
        if(!ha.IsValid()) ha=edi.InsertAtom(res,hydrogen_names_[i][0],hpos,ele);
        else edi.SetAtomPos(ha,hpos);

        R=geom::AxisRotation(v,M_PI);
        direction=R*direction;
        hpos=pa+bond_length*direction;
        ha = res.FindAtom(hydrogen_names_[i][1]);
        if(!ha.IsValid()) ha=edi.InsertAtom(res,hydrogen_names_[i][1],hpos,ele);
        else edi.SetAtomPos(ha,hpos);

        break;
      }
      case 3:{
        bond_angle=0.4*M_PI;

        R=geom::AxisRotation(n,bond_angle);
        direction=R*v;
        hpos=pa+bond_length*direction;
        ha = res.FindAtom(hydrogen_names_[i][0]);
        if(!ha.IsValid()) ha=edi.InsertAtom(res,hydrogen_names_[i][0],hpos,ele);
        else edi.SetAtomPos(ha,hpos);

        R=geom::AxisRotation(v,2.*M_PI/3.);
        direction=R*direction;
        hpos=pa+bond_length*direction;
        ha = res.FindAtom(hydrogen_names_[i][1]);
        if(!ha.IsValid()) ha=edi.InsertAtom(res,hydrogen_names_[i][1],hpos,ele);
        else edi.SetAtomPos(ha,hpos);

        direction=R*direction;
        hpos=pa+bond_length*direction;
        ha = res.FindAtom(hydrogen_names_[i][2]);
        if(!ha.IsValid()) ha=edi.InsertAtom(res,hydrogen_names_[i][2],hpos,ele);
        else edi.SetAtomPos(ha,hpos);
        break;
      }
    }
  }
}

void HeuristicBlockModifier::ApplyOnBuildingBlock(BuildingBlockPtr p){

}

void HeuristicBlockModifier::ApplyOnResidue(ost::mol::ResidueHandle& res, ost::mol::XCSEditor& ed){

}


}}}//ns


















