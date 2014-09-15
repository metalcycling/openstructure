#include <ost/mol/mm/topology_creator.hh>
#include <ost/mol/mm/mm_modeller.hh>

namespace ost{ namespace mol{ namespace mm{

TopologyPtr TopologyCreator::Create(const ost::mol::EntityHandle& handle, 
                                    const MMSettingsPtr settings){

  if(settings->constrain_hangles == true && settings->constrain_hbonds == false){
    throw ost::Error("If hangles is true, hbonds must also be true in settings object!");
  }

  //make sure the input gets not modified
  ost::mol::EntityHandle ent = handle.Copy();
  ost::mol::ResidueHandleList res_list = ent.GetResidueList();
  ost::mol::XCSEditor ed = ent.EditXCS(ost::mol::BUFFERED_EDIT);


  //rename all the stuff to the gromacs naming...
  MMModeller::AssignGromacsNaming(ent);


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
    MMModeller::GenerateDisulfidBonds(ent);
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

  TopologyPtr top = TopologyPtr(new Topology(ent,initial_masses));
  ent = top->GetEntity();

  //note, that we have to get the residue list again, since there is a new entity handle
  //created when initializing the topology
  res_list = ent.GetResidueList();
  ost::mol::AtomHandleList atom_list = ent.GetAtomList();

  //let's extract atom specific informations
  //extract information about single atoms
  std::vector<std::vector<uint> > bound_to;
  std::vector<String> atom_types;
  std::vector<Real> atom_charges;
  std::vector<String> atom_names;
  std::vector<Real> atom_masses;
  std::vector<String> atom_elements;
  std::vector<String> residue_names_of_atoms;
  ost::mol::AtomHandleList bonded_atoms;
  std::vector<uint> temp;
  //extract masses, types, charges and bondpartners

  ost::mol::ResidueHandle temp_res;

  for(ost::mol::AtomHandleList::iterator i = atom_list.begin();
      i!=atom_list.end();++i){
    temp_res = i->GetResidue();
    residue_names_of_atoms.push_back(temp_res.GetName());
    atom_elements.push_back(i->GetElement());
    atom_types.push_back(building_blocks[top->GetResidueIndex(temp_res)]->GetType(i->GetName()));
    atom_charges.push_back(building_blocks[top->GetResidueIndex(temp_res)]->GetCharge(i->GetName()));
    //atom_masses.push_back(building_blocks[top->GetResidueIndex(temp_res)]->GetMass(i->GetName()));
    //if( atom_masses.back() != atom_masses.back() ){
    //  atom_masses.back() = ff->GetMass(atom_types.back());  
    //}
    atom_masses.push_back(ff->GetMass(atom_types.back()));
    bonded_atoms = i->GetBondPartners();
    temp.clear();
    for(ost::mol::AtomHandleList::iterator j = bonded_atoms.begin();
        j!=bonded_atoms.end();++j){
      temp.push_back(top->GetAtomIndex(*j));
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
  std::vector<MMInteractionPtr> ljs;
  std::vector<MMInteractionPtr> gbsa;
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
  std::vector<MMInteractionPtr> interaction_list;
  std::vector<String> interaction_atom_names;

  if(settings->add_impropers){
    for(uint i = 0; i < building_blocks.size(); ++i){
      interaction_list = building_blocks[i]->GetImpropers();
      for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        interaction_atom_names = (*j)->GetNames();
        Index<4> improper_index;
        for(uint k = 0; k < 4; ++k){
          improper_index[k] = top->GetAtomIndex(i,interaction_atom_names[k]);
        }
        impropers.insert(improper_index);
      }
    }
  }

  if(settings->add_cmaps){
    for(uint  i = 0; i < building_blocks.size(); ++i){
      interaction_list = building_blocks[i]->GetCMaps();
      for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        interaction_atom_names = (*j)->GetNames();
        Index<5> cmap_index;
        for(uint k = 0; k < 5; ++k){
          cmap_index[k] = top->GetAtomIndex(i,interaction_atom_names[k]);
        }
        cmaps.insert(cmap_index);
      }
    }
  }
  
  if(settings->add_nonbonded && settings->add_exclusions){
    for(uint i = 0; i < building_blocks.size(); ++i){
      interaction_list = building_blocks[i]->GetExclusions();
      for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        interaction_atom_names = (*j)->GetNames();
        uint one = top->GetAtomIndex(i,interaction_atom_names[0]);
        uint two = top->GetAtomIndex(i,interaction_atom_names[1]);
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
    for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
        j != interaction_list.end(); ++j){
      interaction_atom_names = (*j)->GetNames();
      uint one = top->GetAtomIndex(i, interaction_atom_names[0]);
      uint two = top->GetAtomIndex(i, interaction_atom_names[1]);
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

  //add distance constrains given the corresponding settings
  if(settings->constrain_bonds || settings->constrain_hbonds || settings->rigid_water){
    for(std::set<Index<2> >::iterator i = bonds.begin();
        i != bonds.end(); ++i){
      if(settings->constrain_bonds){
        if(distance_constraints.find(*i) != distance_constraints.end()) continue;
        distance_constraints.insert(*i);
        Real distance;
        if(settings->ideal_bond_length_constraints){
          MMInteractionPtr bond_ptr = ff->GetBond(atom_types[(*i)[0]],atom_types[(*i)[1]]);
          parameters = bond_ptr->GetParam();
          distance = parameters[0];
        }
        else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[1]].GetPos())/10;
        top->AddDistanceConstraint((*i)[0],(*i)[1],distance);
        continue;
      }
      if(settings->constrain_hbonds){
        //if(atom_elements[(*i)[0]] == "H" || atom_elements[(*i)[1]] == "H"){
        if(atom_masses[(*i)[0]] < 1.1 || atom_masses[(*i)[1]] < 1.1){
          if(distance_constraints.find(*i) != distance_constraints.end()) continue;
          distance_constraints.insert(*i);
          Real distance;
          if(settings->ideal_bond_length_constraints){
            MMInteractionPtr bond_ptr = ff->GetBond(atom_types[(*i)[0]],atom_types[(*i)[1]]);
            parameters = bond_ptr->GetParam();
            distance = parameters[0];            
          }
          else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[1]].GetPos())/10;
          top->AddDistanceConstraint((*i)[0],(*i)[1],distance);
          continue;
        }
      }
      if(settings->rigid_water){
        if(residue_names_of_atoms[(*i)[0]] == "SOL" || residue_names_of_atoms[(*i)[1]] == "SOL"){
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
      if(residue_names_of_atoms[(*i)[0]] == "SOL" || residue_names_of_atoms[(*i)[0]] == "SOL" ||
         residue_names_of_atoms[(*i)[2]] == "SOL"){
        //even for ideal_bond_length_constraints, we calculate the ideal angle every time...
        //could be replaced...
        Real distance;
        if(settings->ideal_bond_length_constraints) distance = 0.15139; //HH distance taken from CHARMM      
        else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[2]].GetPos())/10;
        top->AddDistanceConstraint((*i)[0],(*i)[2],distance);
        constrained_angles.insert(*i); 
      }
    }
  }

  if(settings->constrain_hangles){
    for(std::set<Index<3> >::iterator i = angles.begin();
        i != angles.end(); ++i){
      if(atom_masses[(*i)[0]] < 1.1 && atom_masses[(*i)[2]] < 1.1){
        //two hydrogens...
        if(constrained_angles.find(*i) != constrained_angles.end()) continue;
        Real distance;
        if(settings->ideal_bond_length_constraints){
          MMInteractionPtr bond_one = ff->GetBond(atom_types[(*i)[0]],atom_types[(*i)[1]]);
          MMInteractionPtr bond_two = ff->GetBond(atom_types[(*i)[1]],atom_types[(*i)[2]]);
          MMInteractionPtr angle = ff->GetAngle(atom_types[(*i)[0]],atom_types[(*i)[1]],atom_types[(*i)[2]]);
          std::vector<Real> parameters;
          Real l1,l2,a;
          parameters = bond_one->GetParam();
          l1 = parameters[0];
          parameters = bond_two->GetParam();
          l2 = parameters[0];
          parameters = angle->GetParam();
          a = parameters[0];
          distance = sqrt(l1*l1+l2*l2-2*l1*l2*cos(a));
        }      
        else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[2]].GetPos())/10;
        top->AddDistanceConstraint((*i)[0],(*i)[2],distance);
        constrained_angles.insert(*i); 
        continue;
      }

      if(atom_masses[(*i)[1]] > 15.0 && atom_masses[(*i)[1]] < 17.0){
        //central atom is an oxygen
        if(atom_masses[(*i)[0]] < 1.1 || atom_masses[(*i)[2]] < 1.1){
          //a hydrogen is attached to the oxygen...
          if(constrained_angles.find(*i) != constrained_angles.end()) continue;
          Real distance;
          if(settings->ideal_bond_length_constraints){
            MMInteractionPtr bond_one = ff->GetBond(atom_types[(*i)[0]],atom_types[(*i)[1]]);
            MMInteractionPtr bond_two = ff->GetBond(atom_types[(*i)[1]],atom_types[(*i)[2]]);
            MMInteractionPtr angle = ff->GetAngle(atom_types[(*i)[0]],atom_types[(*i)[1]],atom_types[(*i)[2]]);
            std::vector<Real> parameters;
            Real l1,l2,a;
            parameters = bond_one->GetParam();
            l1 = parameters[0];
            parameters = bond_two->GetParam();
            l2 = parameters[0];
            parameters = angle->GetParam();
            a = parameters[0];
            distance = sqrt(l1*l1+l2*l2-2*l1*l2*cos(a));
          }      
          else distance = geom::Distance(atom_list[(*i)[0]].GetPos(),atom_list[(*i)[2]].GetPos())/10;
          top->AddDistanceConstraint((*i)[0],(*i)[2],distance);
          constrained_angles.insert(*i); 
          continue;
        }
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
    for(uint i = 0; i < top->GetNumResidues(); ++i){
      interaction_list = building_blocks[i]->GetBonds();
      for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          interaction_atom_names = (*j)->GetNames();
          uint one = top->GetAtomIndex(i, interaction_atom_names[0]);
          uint two = top->GetAtomIndex(i, interaction_atom_names[1]); 
          Index<2> bond_index(std::min(one,two),std::max(one,two));
          bonds.erase(bond_index);
          //There are only harmonic bonds supported
          parameters = (*j)->GetParam();
          top->AddHarmonicBond(one,two,parameters[0],parameters[1]);
        }
      }
    }
    //add pointers of bond definitions, that are not already parametrized
    //in the building blocks
    for(std::set<Index<2> >::iterator i = bonds.begin();
        i!=bonds.end(); ++i){    

      MMInteractionPtr bond_ptr;
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
      top->AddHarmonicBond((*i)[0],(*i)[1],parameters[0],parameters[1]);
    }
  }

  if(settings->add_angles){
    //handle angles
    for(uint i = 0; i < top->GetNumResidues(); ++i){
      interaction_list = building_blocks[i]->GetAngles();
      for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          interaction_atom_names = (*j)->GetNames();
          uint one = top->GetAtomIndex(i, interaction_atom_names[0]);
          uint two = top->GetAtomIndex(i, interaction_atom_names[1]); 
          uint three = top->GetAtomIndex(i, interaction_atom_names[2]);
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
              top->AddHarmonicAngle(angle_index[0],angle_index[1],angle_index[2],
                                     parameters[0], parameters[1]);
              break;
            }
            case UreyBradleyAngle:{
              top->AddUreyBradleyAngle(angle_index[0],angle_index[1],angle_index[2],
                                        parameters[0], parameters[1], parameters[2],
                                        parameters[3]);
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

      MMInteractionPtr angle_ptr;
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
          top->AddHarmonicAngle((*i)[0],(*i)[1],(*i)[2],
                                 parameters[0],parameters[1]);
          break;
        }
        case UreyBradleyAngle:{
          top->AddUreyBradleyAngle((*i)[0],(*i)[1],(*i)[2],
                                    parameters[0],parameters[1],
                                    parameters[2],parameters[3]);
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
    for(uint i = 0; i < top->GetNumResidues(); ++i){
      interaction_list = building_blocks[i]->GetDihedrals();
      for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          interaction_atom_names = (*j)->GetNames();
          uint one = top->GetAtomIndex(i, interaction_atom_names[0]);
          uint two = top->GetAtomIndex(i, interaction_atom_names[1]); 
          uint three = top->GetAtomIndex(i, interaction_atom_names[2]);
          uint four = top->GetAtomIndex(i, interaction_atom_names[3]);
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
          dihedrals.erase(dihedral_index);
          //only periodic dihedrals are supported... 
          parameters = (*j)->GetParam();
          top->AddPeriodicDihedral(one,two,three,four,
                                   parameters[0],parameters[1],parameters[2]);
        }
      }
    }

    //add pointers of dihedrals definitions, that are not already parametrized
    //in the building blocks
    for(std::set<Index<4> >::iterator i = dihedrals.begin();
        i!=dihedrals.end(); ++i){

      std::vector<MMInteractionPtr> dihedral_ptr;

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

      for(std::vector<MMInteractionPtr>::iterator j = dihedral_ptr.begin();
          j != dihedral_ptr.end(); ++j){
        parameters = (*j)->GetParam();
        top->AddPeriodicDihedral((*i)[0],(*i)[1],(*i)[2],(*i)[3],
                                  parameters[0],parameters[1],parameters[2]);
      }
    }
  }

  if(settings->add_impropers){
    //handle impropers
    for(uint i = 0; i < building_blocks.size(); ++i){
      interaction_list = building_blocks[i]->GetImpropers();
      for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          //we do not have to care about the ordering, since we extracted the
          //impropers from the building block => they will have the same ordering 
          //anyway  (or at least hopefully ;) 
          interaction_atom_names = (*j)->GetNames();
          uint one = top->GetAtomIndex(i, interaction_atom_names[0]);
          uint two = top->GetAtomIndex(i, interaction_atom_names[1]); 
          uint three = top->GetAtomIndex(i, interaction_atom_names[2]);
          uint four = top->GetAtomIndex(i, interaction_atom_names[3]);
          Index<4> improper_index(one,two,three,four);
          impropers.erase(improper_index);
          parameters = (*j)->GetParam();
          switch((*j)->GetFuncType()){
            case PeriodicImproper:{
              top->AddPeriodicImproper(one,two,three,four,
                                       parameters[0],parameters[1],parameters[2]);
              break;
            }
            case HarmonicImproper:{
              top->AddHarmonicImproper(one,two,three,four,
                                       parameters[0],parameters[1]);
              break;
            }
            default:{
              throw ost::Error("Observed invalid function type when adding improper interaction!");
            }
          }
        }
      }
    }
    //add pointers of improper definitions, that are not already parametrized
    //in the building blocks

    for(std::set<Index<4> >::iterator i = impropers.begin();
        i!=impropers.end(); ++i){

      std::vector<MMInteractionPtr> improper_ptr;

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

      for(std::vector<MMInteractionPtr>::iterator j = improper_ptr.begin();
          j != improper_ptr.end(); ++j){
        parameters = (*j)->GetParam();
        switch((*j)->GetFuncType()){
          case PeriodicImproper:{
            top->AddPeriodicImproper((*i)[0],(*i)[1],(*i)[2],(*i)[3],
                                      parameters[0],parameters[1],parameters[2]);
            break;
          }
          case HarmonicImproper:{
            top->AddHarmonicImproper((*i)[0],(*i)[1],(*i)[2],(*i)[3],
                                      parameters[0],parameters[1]);
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
    for(uint i = 0; i < top->GetNumResidues(); ++i){
      interaction_list = building_blocks[i]->GetCMaps();
      for(std::vector<MMInteractionPtr>::iterator j = interaction_list.begin();
          j != interaction_list.end(); ++j){
        if((*j)->IsParametrized()){
          //we do not have to care about the ordering, since we extracted the
          //cmaps from the building block => they will have the same ordering 
          //anyway
          interaction_atom_names = (*j)->GetNames();
          uint one = top->GetAtomIndex(i, interaction_atom_names[0]);
          uint two = top->GetAtomIndex(i, interaction_atom_names[1]); 
          uint three = top->GetAtomIndex(i, interaction_atom_names[2]);
          uint four = top->GetAtomIndex(i, interaction_atom_names[3]);
          uint five = top->GetAtomIndex(i, interaction_atom_names[4]);
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

      MMInteractionPtr cmap_ptr;

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

    MMInteractionPtr lj_interaction;
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
    for(uint i = 0; i < top->GetNumAtoms(); ++i){
      MMInteractionPtr gbsa_ptr = ff->GetImplicitGenborn(atom_types[i]);
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


  return top;

}

}}}//ns
