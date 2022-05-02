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

#include <ost/mol/mm/topology_creator.hh>
#include <ost/mol/mm/modeller.hh>

namespace{

int GetAtomIndex(uint residue_index, 
                 ost::mol::ResidueHandleList& res_list, 
                 std::map<long,int>& atom_indices,
                 const String& atom_name){

  if(atom_name[0] == '+'){
    ost::mol::ResidueHandle actual_res = res_list[residue_index];
    ost::mol::ResidueHandle next_residue = actual_res.GetNext();
    if(!next_residue){
      throw ost::Error("Could not find required atom when settings up topology!");
    }
    String next_atom_name = atom_name.substr(1);
    ost::mol::AtomHandle at = next_residue.FindAtom(next_atom_name);
    if(!at.IsValid()){
      throw ost::Error("Could not find required atom when settings up topology!");      
    }
    return atom_indices[at.GetHashCode()];
  }

  if(atom_name[0] == '-'){
    ost::mol::ResidueHandle actual_res = res_list[residue_index];
    ost::mol::ResidueHandle prev_residue = actual_res.GetPrev();
    if(!prev_residue){
      throw ost::Error("Could not find required atom when settings up topology!");
    }
    String prev_atom_name = atom_name.substr(1);
    ost::mol::AtomHandle at = prev_residue.FindAtom(prev_atom_name);
    if(!at.IsValid()){
      throw ost::Error("Could not find required atom when settings up topology!");      
    }
    return atom_indices[at.GetHashCode()];
  }

  ost::mol::AtomHandle at = res_list[residue_index].FindAtom(atom_name);
  if(!at.IsValid()){
    throw ost::Error("Could not find required atom when settings up topology!");      
  }
  return atom_indices[at.GetHashCode()];
}






}

namespace ost{ namespace mol{ namespace mm{

TopologyPtr TopologyCreator::Create(ost::mol::EntityHandle& ent, 
                                    const SettingsPtr settings){


  ost::mol::ResidueHandleList res_list = ent.GetResidueList();
  ost::mol::XCSEditor ed = ent.EditXCS(ost::mol::BUFFERED_EDIT);

  //rename all the stuff to the gromacs naming...
  Modeller::AssignGromacsNaming(ent);

  //even if not reconnected yet, it gets assumed, that
  //peptide or nucleotide bonds are correctly set in the input entity
  //to allow the identification and tagging of terminal residues
  //note, that only nucleotides or aminoacids are allowed to build termini...
  //this should be enough for most needs

  ost::mol::ResidueHandle next,prev;
  bool n_ter,c_ter; 
  ost::mol::AtomHandle peptide_n,peptide_c,nucleotide_p,nucleotide_o; 

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin(); 
      i != res_list.end(); ++i){

    n_ter = false;
    c_ter = false;

    peptide_n = i->FindAtom("N");
    nucleotide_p = i->FindAtom("P");

    if(!peptide_n.IsValid() && !nucleotide_p.IsValid()) continue;
    //in this case the residue is neither a peptide nor a nucleotide.
    //There won't be a termini anyway...

    prev = i->GetPrev();
    next = i->GetNext();

    if(!prev.IsValid()){
      n_ter = true;
    }
    else{
      peptide_c = prev.FindAtom("C");
      nucleotide_o = prev.FindAtom("O3");

      if(!ost::mol::BondExists(peptide_n,peptide_c) && !ost::mol::BondExists(nucleotide_p,nucleotide_o)){
        n_ter = true;
      }
    }

    if(!next.IsValid()){
      c_ter = true;
    }
    else{
      peptide_n = next.FindAtom("N");
      peptide_c = i->FindAtom("C");
      nucleotide_p = next.FindAtom("P");
      nucleotide_o = i->FindAtom("O3");

      if(!ost::mol::BondExists(peptide_n,peptide_c) && !ost::mol::BondExists(nucleotide_p,nucleotide_o)){
        c_ter = true;
      }
    }

    if(n_ter) i->SetBoolProp("n_ter",true);
    if(c_ter) i->SetBoolProp("c_ter",true);
  }


  ForcefieldPtr ff = settings->forcefield;
  if(!ff){
    throw ost::Error("You have to assign a valid forcefield to the settings object when creating a topology!");
  }
  //The full parametrization procedure relies on the naming, that is force field specific
  
  ost::mol::BondHandleList bond_list = ent.GetBondList();
  ed.DeleteBonds(bond_list);
  ed.UpdateICS();

  if(settings->generate_disulfid_bonds){
    Modeller::GenerateDisulfidBonds(ent);
  }

  ff->AssignFFSpecificNames(ent);

  //we first generate the building blocks and directly construct hydrogens/termini and stuff
  std::vector<BuildingBlockPtr> building_blocks;

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin(); 
      i != res_list.end(); ++i){
    BuildingBlockPtr block = ff->GetBuildingBlock(i->GetName());
    HydrogenConstructorPtr hc = ff->GetHydrogenConstructor(i->GetName());
    if(hc){
      hc->ApplyOnBuildingBlock(block);
      hc->ApplyOnResidue(*i,ed);
    }
    //check for n terminus
    if(i->HasProp("n_ter")){
      String exception_name = "";
      if(settings->termini_exceptions->HasException(*i)){
        exception_name = settings->termini_exceptions->GetException(*i);
      }
      BlockModifierPtr bm = ff->GetNTerModifier(i->GetName(), exception_name);
      if(bm){
        bm->ApplyOnBuildingBlock(block);
        bm->ApplyOnResidue(*i,ed);
        block->RemoveInteractionsToPrev();
      }
    }
    if(i->HasProp("c_ter")){
      String exception_name = "";
      if(settings->termini_exceptions->HasException(*i)){
        exception_name = settings->termini_exceptions->GetException(*i);
      }
      BlockModifierPtr bm = ff->GetCTerModifier(i->GetName(), exception_name);
      if(bm){
        bm->ApplyOnBuildingBlock(block);
        bm->ApplyOnResidue(*i,ed);
        block->RemoveInteractionsToNext();
      }
    }
    block->Connect(*i,ed);
    building_blocks.push_back(block);
  }

  //The editor won't be needed anymore, let's update the internal coordinate system
  ed.UpdateICS();


  //check, whether all building blocks match now with the previously applied modifications
  String match_fail_info;
  for(unsigned int i = 0; i < building_blocks.size(); ++i){
    if(!building_blocks[i]->Match(res_list[i], false, match_fail_info)){
      std::stringstream ss;
      ss << "Residue "<< res_list[i].GetQualifiedName() << " does not match the force field ";
      ss << "definition. "<<match_fail_info;
      throw ost::Error(ss.str());
    }
  }

  std::vector<Real> initial_masses(ent.GetAtomCount());

  TopologyPtr top = TopologyPtr(new Topology(initial_masses));

  ost::mol::AtomHandleList atom_list = ent.GetAtomList();

  std::map<long,int> atom_indices;
  std::map<long,int> residue_indices;
  int actual_index = 0;

  for(ost::mol::AtomHandleList::const_iterator i = atom_list.begin(), e = atom_list.end(); 
      i != e; ++i){
    atom_indices[i->GetHashCode()] = actual_index;
    ++actual_index;
  }

  actual_index = 0;
  for(ost::mol::ResidueHandleList::const_iterator i = res_list.begin(), e = res_list.end(); 
      i != e; ++i){
    residue_indices[i->GetHashCode()] = actual_index;
    ++actual_index;
  }


  //let's extract atom specific informations
  //extract information about single atoms
  std::vector<std::vector<uint> > bound_to;
  std::vector<String> atom_types;
  std::vector<Real> atom_charges;
  std::vector<String> atom_names;
  std::vector<Real> atom_masses;
  std::vector<String> residue_names_of_atoms;
  ost::mol::AtomHandleList bonded_atoms;
  std::vector<uint> temp;
  int residue_index;

  //extract masses, types, charges and bondpartners

  for(ost::mol::AtomHandleList::iterator i = atom_list.begin();
      i!=atom_list.end();++i){
    residue_index = residue_indices[i->GetResidue().GetHashCode()];
    residue_names_of_atoms.push_back(i->GetResidue().GetName());
    atom_types.push_back(building_blocks[residue_index]->GetType(i->GetName()));
    atom_charges.push_back(building_blocks[residue_index]->GetCharge(i->GetName()));
    atom_names.push_back(i->GetName());
    atom_masses.push_back(ff->GetMass(atom_types.back()));
    bonded_atoms = i->GetBondPartners();
    temp.clear();
    for(ost::mol::AtomHandleList::iterator j = bonded_atoms.begin();
        j!=bonded_atoms.end();++j){
      temp.push_back(atom_indices[j->GetHashCode()]);
    }
    bound_to.push_back(temp);
  }

  //let's set the proper masses
  top->SetMasses(atom_masses);


  //we simply set all charges to zero if we want to kill the electrostatics
  if(settings->kill_electrostatics){
    for(std::vector<Real>::iterator i = atom_charges.begin();
        i != atom_charges.end(); ++i){
      *i = 0;
    }
  }

  top->SetCharges(atom_charges);

  //the interaction read from the structure will be stored in here
  std::set<Index<2> > bonds;
  std::set<Index<3> > angles;
  std::set<Index<4> > dihedrals;
  std::set<Index<4> > impropers;
  std::set<Index<5> > cmaps;
  std::set<Index<2> > exclusions;
  std::set<Index<2> > pairs_14;
  std::set<Index<2> > distance_constraints;
  //per atom parameters
  std::vector<InteractionPtr> ljs;
  std::vector<InteractionPtr> gbsa;
  //extract all bonded interactions directly from the structure itself,
  //except the impropers, cmaps and exclusions. Those two get read from the building
  //blocks.
  //Note, that the bonds build the basis to build angles and the angles 
  //build the basis to build dihedrals... so we have to extract bonds 
  //even when only dihedrals have to be added.
  //Another special case is when nonbonded terms have to be added.
  //They often exclude interactions originating from closely bonded
  //atoms. We will use the bonds, angles and dihedrals later on to
  //define these closely bonded atoms, so if nonbonded terms have to
  //added, we also need to get the bonds, angles and dihedrals.
  if(settings->add_bonds || settings->add_angles || 
     settings->add_dihedrals || settings->add_nonbonded){
    //build all bonds
    for(uint i=0 ; i<bound_to.size() ;++i){
      for(uint j = 0; j<bound_to[i].size(); ++j){
        bonds.insert(Index<2>(std::min(i,bound_to[i][j]),std::max(i,bound_to[i][j])));
      }
    }
  }
  if(settings->add_angles || settings->add_dihedrals || 
     settings->add_nonbonded){
    //build all angles
    for(std::set<Index<2> >::iterator i = bonds.begin();
        i!=bonds.end();++i){
      uint atom_1 = (*i)[0];
      uint atom_2 = (*i)[1];
      
      for(std::vector<uint>::iterator j = bound_to[atom_1].begin();
          j!=bound_to[atom_1].end(); ++j){
        if((*j) != atom_2){
          angles.insert(Index<3>(std::min((*j),atom_2),atom_1,std::max((*j),atom_2)));
        }
      }
      for(std::vector<uint>::iterator j = bound_to[atom_2].begin();
          j!=bound_to[atom_2].end(); ++j){
        if((*j) != atom_1){
          angles.insert(Index<3>(std::min((*j),atom_1),atom_2,std::max((*j),atom_1)));
        }
      }
    }
  }
  if(settings->add_dihedrals || settings->add_nonbonded){
    //build all dihedrals
    for(std::set<Index<3> >::iterator i = angles.begin();
        i!=angles.end();++i){
      uint atom_1 = (*i)[0];
      uint atom_2 = (*i)[1];
      uint atom_3 = (*i)[2];
      for(std::vector<uint>::iterator j = bound_to[atom_1].begin();
          j!=bound_to[atom_1].end();++j){
        if((*j)!=atom_2){
          if((*j)<atom_3) dihedrals.insert(Index<4>(*j,atom_1,atom_2,atom_3));
          else dihedrals.insert(Index<4>(atom_3,atom_2,atom_1,*j)); 
        }
      }
      for(std::vector<uint>::iterator j = bound_to[atom_3].begin();
          j!=bound_to[atom_3].end();++j){
        if((*j)!=atom_2){
          if((*j)<atom_1) dihedrals.insert(Index<4>(*j,atom_3,atom_2,atom_1));
          else dihedrals.insert(Index<4>(atom_1,atom_2,atom_3,*j));
        }
      }
    }
  }
  //impropers exclusions and cmaps get read from residuetemplates!
  std::vector<InteractionPtr> interaction_list;
  std::vector<String> interaction_atom_names;

  if(settings->add_impropers){
    for(uint i = 0; i < building_blocks.size(); ++i){
      interaction_list = building_blocks[i]->GetImpropers();
      for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        interaction_atom_names = (*j)->GetNames();
        Index<4> improper_index;
        for(uint k = 0; k < 4; ++k){
          improper_index[k] = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[k]);
        }
        impropers.insert(improper_index);
      }
    }
  }

  if(settings->add_cmaps){
    for(uint  i = 0; i < building_blocks.size(); ++i){
      interaction_list = building_blocks[i]->GetCMaps();
      for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        interaction_atom_names = (*j)->GetNames();
        Index<5> cmap_index;
        for(uint k = 0; k < 5; ++k){
          cmap_index[k] = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[k]);
        }
        cmaps.insert(cmap_index);
      }
    }
  }
  
  if(settings->add_nonbonded && settings->add_exclusions){
    for(uint i = 0; i < building_blocks.size(); ++i){
      interaction_list = building_blocks[i]->GetExclusions();
      for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        interaction_atom_names = (*j)->GetNames();
        uint one = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[0]);
        uint two = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[1]);
        exclusions.insert(Index<2>(std::min(one,two),std::max(one,two)));
      }
    }
  }
  if(settings->add_nonbonded){

    //find exclusions and excpetions based on connectivity
    //add exclusions using previously extracted bonds => 1,2 interactions
    for(std::set<Index<2> >::iterator i = bonds.begin();
        i!=bonds.end(); ++i){
      exclusions.insert(*i);
    }
    //add exclusions using previously extracted angles => 1,3 interactions
    for(std::set<Index<3> >::iterator i = angles.begin();
        i!=angles.end(); ++i){
      exclusions.insert(Index<2>((*i)[0],(*i)[2]));
    }
    //add exclusions using previously extracted dihedrals => 1,4 interactions
    for(std::set<Index<4> >::iterator i = dihedrals.begin();
        i!=dihedrals.end(); ++i){
      pairs_14.insert(Index<2>((*i)[0],(*i)[3]));
    }

    //remove all exclusions from the 1,4 interactions
    for(std::set<Index<2> >::iterator i = exclusions.begin();
        i != exclusions.end(); ++i){
      pairs_14.erase(*i);
    }
  }

  //we will need this goddamn variable in a lot of places...
  std::vector<Real> parameters;


  //add constraints from building blocks
  //note, that the constraints are in some cases parametrized.
  //If there is no constraint distance specified, the actual distance
  //between the corresponding atoms is taken.
  //if there are several constraints on the same atom pair,
  //only the first is added to avoid contradicting constraints
  for(uint i = 0; i < building_blocks.size(); ++i){
    interaction_list = building_blocks[i]->GetConstraints();
    for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
        j != interaction_list.end(); ++j){
      interaction_atom_names = (*j)->GetNames();
      uint one = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[0]);
      uint two = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[1]);
      Index<2> index(std::min(one,two),std::max(one,two));
      //we don't want contradicting constraints!
      if(distance_constraints.find(index) != distance_constraints.end()) continue;
      distance_constraints.insert(index);
      Real distance;
      if((*j)->IsParametrized()){
        parameters = (*j)->GetParam();
        distance = parameters[0];
      }
      else distance = geom::Distance(atom_list[one].GetPos(),atom_list[two].GetPos())/10;
      top->AddDistanceConstraint(one,two,distance);
    }
  }

  std::set<String> water_res_names;
  water_res_names.insert("SOL");
  water_res_names.insert("HOH");
  water_res_names.insert("TIP3");

  //add distance constrains given the corresponding settings
  if(settings->constrain_bonds || settings->constrain_hbonds || settings->rigid_water){
    for(std::set<Index<2> >::iterator i = bonds.begin();
        i != bonds.end(); ++i){
      if(settings->constrain_bonds){
        if(distance_constraints.find(*i) != distance_constraints.end()) continue;
        distance_constraints.insert(*i);
        Real distance;
        if(settings->ideal_bond_length_constraints){
          InteractionPtr bond_ptr = ff->GetBond(atom_types[(*i)[0]],atom_types[(*i)[1]]);
          parameters = bond_ptr->GetParam();
          distance = parameters[0];
        }
        else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[1]].GetPos())/10;
        top->AddDistanceConstraint((*i)[0],(*i)[1],distance);
        continue;
      }
      if(settings->constrain_hbonds){
        if(atom_names[(*i)[0]][0] == 'H' || atom_names[(*i)[1]][0] == 'H'){
          if(distance_constraints.find(*i) != distance_constraints.end()) continue;
          distance_constraints.insert(*i);
          Real distance;
          if(settings->ideal_bond_length_constraints){
            InteractionPtr bond_ptr = ff->GetBond(atom_types[(*i)[0]],atom_types[(*i)[1]]);
            parameters = bond_ptr->GetParam();
            distance = parameters[0];            
          }
          else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[1]].GetPos())/10;
          top->AddDistanceConstraint((*i)[0],(*i)[1],distance);
          continue;
        }
      }
      if(settings->rigid_water){
        if(water_res_names.find(residue_names_of_atoms[(*i)[0]]) != water_res_names.end() && 
           water_res_names.find(residue_names_of_atoms[(*i)[1]]) != water_res_names.end()){
          if(distance_constraints.find(*i) != distance_constraints.end()) continue;
          distance_constraints.insert(*i);
          Real distance;
          if(settings->ideal_bond_length_constraints) distance = 0.09572; //OH bond length in water in CHARMM forcefield...
          else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[1]].GetPos())/10;
          top->AddDistanceConstraint((*i)[0],(*i)[1],distance);
          continue;        
        }
      }
    }
  }

  //set angle constraints given the corresponding settings
  std::set<Index<3> > constrained_angles;
  if(settings->rigid_water){
    for(std::set<Index<3> >::iterator i = angles.begin();
        i != angles.end(); ++i){
      if(water_res_names.find(residue_names_of_atoms[(*i)[0]]) != water_res_names.end() && 
         water_res_names.find(residue_names_of_atoms[(*i)[1]]) != water_res_names.end() &&
         water_res_names.find(residue_names_of_atoms[(*i)[2]]) != water_res_names.end()){
        //we only have to add the H-H distance constant, the O-H distance is already
        //constrained above
        Real distance;
        if(settings->ideal_bond_length_constraints) distance = 0.15139; //HH distance taken from CHARMM      
        else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[2]].GetPos())/10;
        top->AddDistanceConstraint((*i)[0],(*i)[2],distance);
        constrained_angles.insert(*i); 
      }
    }
  }
  
  //remove all constraints from the bonds and angles
  for(std::set<Index<2> >::iterator i = distance_constraints.begin();
      i != distance_constraints.end(); ++i){
    bonds.erase(*i);
  }

  for(std::set<Index<3> >::iterator i = constrained_angles.begin(); 
      i!= constrained_angles.end(); ++i){
    angles.erase(*i);
  }


  //the force definitions from the forcefield can be overwritten by parameters
  //defined in the residue building blocks.
  //We therefore check which forces have been defined in the building blocks
  //and don't search them any more in the forcefield forces.
  std::vector<String> types;
  if(settings->add_bonds){
    //handle bonds
    for(uint i = 0; i < res_list.size(); ++i){
      interaction_list = building_blocks[i]->GetBonds();
      for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          interaction_atom_names = (*j)->GetNames();
          uint one = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[0]);
          uint two = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[1]); 
          Index<2> bond_index(std::min(one,two),std::max(one,two));
          bonds.erase(bond_index);
          //There are only harmonic bonds supported
          parameters = (*j)->GetParam();
          if(parameters[1] != 0.0) {
            top->AddHarmonicBond(one,two,parameters[0],parameters[1]);
          }
        }
      }
    }
    //add pointers of bond definitions, that are not already parametrized
    //in the building blocks
    for(std::set<Index<2> >::iterator i = bonds.begin();
        i!=bonds.end(); ++i){    

      InteractionPtr bond_ptr;
      try{
        bond_ptr = ff->GetBond(atom_types[(*i)[0]],atom_types[(*i)[1]]);
      } catch(ost::Error& e){
        if(settings->strict_interactions){
          std::stringstream ss;
          ss << "Failed to parametrize bond between: " << atom_list[(*i)[0]];
          ss << " and "<<atom_list[(*i)[1]] << ". ";
          ss << "To ignore these things, set strict_interactions in the settings";
          ss << " object to False (not recommended).";
          throw ost::Error(ss.str());
        }
        continue;  // ignore it and continue with next interaction
      }

      parameters = bond_ptr->GetParam();
      if(parameters[1] != 0.0) {
        top->AddHarmonicBond((*i)[0],(*i)[1],parameters[0],parameters[1]);
      }
    }
  }

  if(settings->add_angles){
    //handle angles
    for(uint i = 0; i < res_list.size(); ++i){
      interaction_list = building_blocks[i]->GetAngles();
      for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          interaction_atom_names = (*j)->GetNames();
          uint one = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[0]);
          uint two = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[1]); 
          uint three = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[2]);
          Index<3> angle_index(std::min(one,three),two,std::max(one,three));
          if(angles.find(angle_index) == angles.end()){
            std::stringstream ss;
            ss << "Building block for residue " << res_list[i].GetQualifiedName();
            ss << " defines angle, that doesn't exist!";
            throw ost::Error(ss.str());
          }
          angles.erase(angle_index);
          parameters = (*j)->GetParam();
          switch((*j)->GetFuncType()){
            case HarmonicAngle:{
              if(parameters[1] != 0.0) {
                top->AddHarmonicAngle(angle_index[0],angle_index[1],angle_index[2],
                                      parameters[0], parameters[1]);
              }
              break;
            }
            case UreyBradleyAngle:{
              if(parameters[1] != 0.0 || parameters[3] != 0.0) {
                top->AddUreyBradleyAngle(angle_index[0],angle_index[1],angle_index[2],
                                          parameters[0], parameters[1], parameters[2],
                                          parameters[3]);
              }
              break;
            }
            default:{
              throw ost::Error("Observed unknown function type for angle interaction!");
            }
          }
        }
      }
    }
    //add pointers of angle definitions, that are not already parametrized
    //in the building blocks
    for(std::set<Index<3> >::iterator i = angles.begin();
        i!=angles.end(); ++i){

      InteractionPtr angle_ptr;
      try{
        angle_ptr = ff->GetAngle(atom_types[(*i)[0]],
                                 atom_types[(*i)[1]],
                                 atom_types[(*i)[2]]);
      } catch(ost::Error& e){
        if(settings->strict_interactions){
          std::stringstream ss;
          ss << "Failed to parametrize angle between: " << atom_list[(*i)[0]];
          ss << ", "<<atom_list[(*i)[1]] << " and " << atom_list[(*i)[2]] << ".";
          ss << "To ignore these things, set strict_interactions in the settings";
          ss << " object to False (not recommended).";
          throw ost::Error(ss.str());
        }
        continue;  // ignore it and continue with next interaction
      }

      parameters = angle_ptr->GetParam();
      switch(angle_ptr->GetFuncType()){
        case HarmonicAngle:{
          if(parameters[1] != 0.0) {
            top->AddHarmonicAngle((*i)[0],(*i)[1],(*i)[2],
                                   parameters[0],parameters[1]);
          }
          break;
        }
        case UreyBradleyAngle:{
          if(parameters[1] != 0.0 || parameters[3] != 0.0) {
            top->AddUreyBradleyAngle((*i)[0],(*i)[1],(*i)[2],
                                      parameters[0],parameters[1],
                                      parameters[2],parameters[3]);
          }
          break;
        }
        default:{
          throw ost::Error("Observed invalid function type for angle interaction!");
        }
      }
    }
  }

  if(settings->add_dihedrals){
    //handle dihedrals
    std::set<Index<4> > dihedrals_to_delete;
    for(uint i = 0; i < res_list.size(); ++i){
      interaction_list = building_blocks[i]->GetDihedrals();
      for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          interaction_atom_names = (*j)->GetNames();
          uint one = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[0]);
          uint two = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[1]); 
          uint three = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[2]);
          uint four = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[3]);
          Index<4> dihedral_index;
          if(one<four){
            dihedral_index[0] = one;
            dihedral_index[1] = two;
            dihedral_index[2] = three;
            dihedral_index[3] = four;
          }
          else{
            dihedral_index[0] = four;
            dihedral_index[1] = three;
            dihedral_index[2] = two;
            dihedral_index[3] = one;
          }
          if(dihedrals.find(dihedral_index) == dihedrals.end()){
            std::stringstream ss;
            ss << "Building block for residue " << res_list[i].GetQualifiedName();
            ss << "defines dihedral, that doesn't exist!";
            throw ost::Error(ss.str());
          }
          
          dihedrals_to_delete.insert(dihedral_index);
          //only periodic dihedrals are supported... 
          parameters = (*j)->GetParam();
          if(parameters[2] != 0.0) {
            top->AddPeriodicDihedral(one,two,three,four,
                                     parameters[0],parameters[1],parameters[2]);
          }
        }
      }
    }
    
    for(std::set<Index<4> >::iterator i = dihedrals_to_delete.begin();
        i != dihedrals_to_delete.end(); ++i){
      dihedrals.erase(*i);
    }

    //add pointers of dihedrals definitions, that are not already parametrized
    //in the building blocks
    for(std::set<Index<4> >::iterator i = dihedrals.begin();
        i!=dihedrals.end(); ++i){

      std::vector<InteractionPtr> dihedral_ptr;

      try{
        dihedral_ptr = ff->GetDihedrals(atom_types[(*i)[0]],
                                        atom_types[(*i)[1]],
                                        atom_types[(*i)[2]],
                                        atom_types[(*i)[3]]);
      } catch(ost::Error& e){
        if(settings->strict_interactions){
          std::stringstream ss;
          ss << "Failed to parametrize dihedral between: " << atom_list[(*i)[0]];
          ss << ", "<<atom_list[(*i)[1]] << ", " << atom_list[(*i)[2]] << " and ";
          ss << atom_list[(*i)[3]] << ".";
          ss << "To ignore these things, set strict_interactions in the settings";
          ss << " object to False (not recommended).";
          throw ost::Error(ss.str());
        }
        continue;  // ignore it and continue with next interaction
      }

      for(std::vector<InteractionPtr>::iterator j = dihedral_ptr.begin();
          j != dihedral_ptr.end(); ++j){
        parameters = (*j)->GetParam();
        if(parameters[2] != 0.0) {
          top->AddPeriodicDihedral((*i)[0],(*i)[1],(*i)[2],(*i)[3],
                                    parameters[0],parameters[1],parameters[2]);
        }
      }
    }
  }

  if(settings->add_impropers){
    //handle impropers
    std::set<Index<4> > impropers_to_delete;
    for(uint i = 0; i < building_blocks.size(); ++i){
      interaction_list = building_blocks[i]->GetImpropers();
      for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          //we do not have to care about the ordering, since we extracted the
          //impropers from the building block => they will have the same ordering 
          //anyway  (or at least hopefully ;) 
          interaction_atom_names = (*j)->GetNames();
          uint one = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[0]);
          uint two = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[1]); 
          uint three = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[2]);
          uint four = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[3]);
          Index<4> improper_index(one,two,three,four);
          impropers_to_delete.insert(improper_index);
          parameters = (*j)->GetParam();
          switch((*j)->GetFuncType()){
            case PeriodicImproper:{
              if(parameters[2] != 0.0) {
                top->AddPeriodicImproper(one,two,three,four,
                                         parameters[0],parameters[1],parameters[2]);
              }
              break;
            }
            case HarmonicImproper:{
              if(parameters[1] != 0.0) {
                top->AddHarmonicImproper(one,two,three,four,
                                         parameters[0],parameters[1]);
              }
              break;
            }
            default:{
              throw ost::Error("Observed invalid function type when adding improper interaction!");
            }
          }
        }
      }
    }
    
    for(std::set<Index<4> >::iterator i = impropers_to_delete.begin();
        i != impropers_to_delete.end(); ++i){
      impropers.erase(*i);
    }
    
    //add pointers of improper definitions, that are not already parametrized
    //in the building blocks

    for(std::set<Index<4> >::iterator i = impropers.begin();
        i!=impropers.end(); ++i){

      std::vector<InteractionPtr> improper_ptr;

      try{
        improper_ptr = ff->GetImpropers(atom_types[(*i)[0]],
                                        atom_types[(*i)[1]],
                                        atom_types[(*i)[2]],
                                        atom_types[(*i)[3]]);
      } catch(ost::Error& e){
        if(settings->strict_interactions){
          std::stringstream ss;
          ss << "Failed to parametrize improper between: " << atom_list[(*i)[0]];
          ss << ", "<<atom_list[(*i)[1]] << ", " << atom_list[(*i)[2]] << " and ";
          ss << atom_list[(*i)[3]] << ".";
          ss << "To ignore these things, set strict_interactions in the settings";
          ss << " object to False (not recommended).";
          throw ost::Error(ss.str());
        }
        continue;  // ignore it and continue with next interaction
      }

      for(std::vector<InteractionPtr>::iterator j = improper_ptr.begin();
          j != improper_ptr.end(); ++j){
        parameters = (*j)->GetParam();
        switch((*j)->GetFuncType()){
          case PeriodicImproper:{
            if(parameters[2] != 0.0) {
              top->AddPeriodicImproper((*i)[0],(*i)[1],(*i)[2],(*i)[3],
                                        parameters[0],parameters[1],parameters[2]);
            }
            break;
          }
          case HarmonicImproper:{
            if(parameters[1] != 0.0) {
              top->AddHarmonicImproper((*i)[0],(*i)[1],(*i)[2],(*i)[3],
                                        parameters[0],parameters[1]);
            }
            break;
          }
          default:{
            throw ost::Error("Observed invalid function type when adding improper interaction!");
          }
        }
      }
    }
  }

  if(settings->add_cmaps){
    //handle cmaps
    for(uint i = 0; i < res_list.size(); ++i){
      interaction_list = building_blocks[i]->GetCMaps();
      for(std::vector<InteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          //we do not have to care about the ordering, since we extracted the
          //cmaps from the building block => they will have the same ordering 
          //anyway
          interaction_atom_names = (*j)->GetNames();
          uint one = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[0]);
          uint two = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[1]); 
          uint three = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[2]);
          uint four = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[3]);
          uint five = GetAtomIndex(i,res_list,atom_indices,interaction_atom_names[4]);
          Index<5> cmap_index(one,two,three,four,five);
          cmaps.erase(cmap_index);
          parameters = (*j)->GetParam();
          int dimension = int(parameters[0]);
          parameters.erase(parameters.begin());
          top->AddCMap(cmap_index[0],cmap_index[1],cmap_index[2],cmap_index[3],
                       cmap_index[4],dimension,parameters);
        }
      }
    }

    //add pointers of cmap definitions, that are not already parametrized
    //in the building blocks
    for(std::set<Index<5> >::iterator i = cmaps.begin();
        i!=cmaps.end(); ++i){

      InteractionPtr cmap_ptr;

      try{
        cmap_ptr = ff->GetCMap(atom_types[(*i)[0]],
                               atom_types[(*i)[1]],
                               atom_types[(*i)[2]],
                               atom_types[(*i)[3]],
                               atom_types[(*i)[4]]);
      } catch(ost::Error& e){
        if(settings->strict_interactions){
          std::stringstream ss;
          ss << "Failed to parametrize cmap between: " << atom_list[(*i)[0]];
          ss << ", "<<atom_list[(*i)[1]] << ", " << atom_list[(*i)[2]] << ", ";
          ss << atom_list[(*i)[3]] << " and " << atom_list[(*i)[4]] << ".";
          ss << "To ignore these things, set strict_interactions in the settings";
          ss << " object to False (not recommended).";
          throw ost::Error(ss.str());
        }
        continue;  // ignore it and continue with next interaction
      }
      parameters = cmap_ptr->GetParam();
      int dimension = int(parameters[0]);
      parameters.erase(parameters.begin());
      top->AddCMap((*i)[0],(*i)[1],(*i)[2],(*i)[3],(*i)[4],dimension,parameters);
    }
  }

  if(settings->add_nonbonded){

    InteractionPtr lj_interaction;
    std::vector<Real> sigmas;
    std::vector<Real> epsilons;

    for(std::vector<String>::iterator i = atom_types.begin();
        i != atom_types.end(); ++i){
      lj_interaction = ff->GetLJ(*i);
      if(!lj_interaction){
        std::stringstream ss;
        ss << "Failed to find LJ parametrization for type \""<<*i<<"\"";
        ss << " in given forcefield!";
        throw ost::Error(ss.str());
      }
      parameters = lj_interaction->GetParam();
      sigmas.push_back(parameters[0]);
      epsilons.push_back(parameters[1]);
    }
    top->SetSigmas(sigmas);
    top->SetEpsilons(epsilons);
    //take care of the 1-4 interactions
    for(std::set<Index<2> >::iterator i = pairs_14.begin();
        i != pairs_14.end(); ++i){
      //we can be sure, that the parameters are valid, since we found all single
      //parameters. The forcefield can therefore build the new parametrization with
      //the combination rules... Except: gen_pairs_ is set to no, but then an error
      // is thrown anyway
      lj_interaction = ff->GetLJ(atom_types[(*i)[0]], 
                                 atom_types[(*i)[1]],
                                 true);
      parameters = lj_interaction->GetParam();
      top->AddLJPair((*i)[0],(*i)[1],parameters[0],parameters[1]);
    }
    //take care of the exclusions
    for(std::set<Index<2> >::iterator i = exclusions.begin();
        i != exclusions.end(); ++i){
      top->AddExclusion((*i)[0],(*i)[1]);
    }
  }
  if(settings->add_gbsa){
    std::vector<Real> radii;
    std::vector<Real> scaling;
    for(uint i = 0; i < atom_list.size(); ++i){
      InteractionPtr gbsa_ptr = ff->GetImplicitGenborn(atom_types[i]);
      if(!gbsa_ptr){
        std::stringstream ss;
        ss << "Structure contains atom of type " << types[i] << " which is not";
        ss << "parametrized in the forcefield!"<<std::endl;
        throw ost::Error(ss.str());
      }
      parameters = gbsa_ptr->GetParam();
      radii.push_back(parameters[0]);
      scaling.push_back(parameters[1]);
    }
    top->SetGBSARadii(radii);
    top->SetOBCScalings(scaling);
  }

  top->SetFudgeQQ(ff->GetFudgeQQ());
  top->SetFudgeLJ(ff->GetFudgeLJ());

  if(!settings->keep_ff_specific_naming){
    //let's rename to the pdb standard!
    //first, we rename to the gromacs standars...
    ff->AssignFFSpecificNames(ent,true);
    //still call the assign gromacs names function...
    Modeller::AssignGromacsNaming(ent);
    //let's finally do the renaming to PDB standard
    Modeller::AssignPDBNaming(ent);
  }


  return top;

}

}}}//ns
