#include <ost/mol/mm/modeller.hh>


namespace ost{ namespace mol{ namespace mm{

void Modeller::GenerateDisulfidBonds(ost::mol::EntityHandle& handle){

  ost::mol::ResidueHandleList res_list = handle.GetResidueList();
  ost::mol::XCSEditor ed = handle.EditXCS();

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();
      i != res_list.end(); ++i){
    if(i->GetName() == "CYS" || i->GetName() == "CYX" || i->GetName() == "CYS2" || i->GetName() == "CYM"){
      ost::mol::AtomHandle s = i->FindAtom("SG");
      if(s.IsValid()){
        //quick check, whether max bond partners is already reached
        if(s.GetBondCount() >= 2) continue;
        //check, whether SG is already bound to a SG
        ost::mol::AtomHandleList bonded_atoms = s.GetBondPartners();
        bool already_bound = false;
        for(ost::mol::AtomHandleList::iterator bound_it = bonded_atoms.begin();
            bound_it != bonded_atoms.end(); ++bound_it){
          if(bound_it->GetName() == "SG"){
            if(bound_it->GetResidue().GetName() == "CYS" || bound_it->GetResidue().GetName() == "CYX" ||
               bound_it->GetResidue().GetName() == "CYS2" || bound_it->GetResidue().GetName() == "CYM"){
              already_bound = true;
              break;
            }
          }
        }
        if(already_bound) continue;
        ost::mol::AtomHandleList in_reach = handle.FindWithin(s.GetPos(),2.5);
        for(ost::mol::AtomHandleList::iterator j = in_reach.begin();
            j != in_reach.end(); ++j){
          if(j->GetName() == "SG" && j->GetResidue() != *i){
            if(j->GetBondCount() >= 2) continue;
            if(!ost::mol::BondExists(s,*j)){
              ed.Connect(s,*j);
            }
            ed.RenameResidue(*i,"CYS2");
            ed.RenameResidue(j->GetResidue(),"CYS2");
          }
        }
      }
    }
  }
}

void Modeller::GenerateCYSHEMEBonds(ost::mol::EntityHandle& handle){
  
  ost::mol::ResidueHandleList res_list = handle.GetResidueList();
  ost::mol::XCSEditor ed = handle.EditXCS();

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();
      i != res_list.end(); ++i){
    if(i->GetName() == "HEM" || i->GetName() == "HEME"){
      ost::mol::AtomHandle fe = i->FindAtom("FE");
      if(fe.IsValid()){
        ost::mol::AtomHandleList in_reach = handle.FindWithin(fe.GetPos(), 3.0);
        for(ost::mol::AtomHandleList::iterator j = in_reach.begin();
            j != in_reach.end(); ++j){
          if( j->GetName() == "SG"){
            if(!ost::mol::BondExists(fe,*j)){
              ed.Connect(fe,*j);
            }
            ed.RenameResidue(*i,"HEME");
            ed.RenameResidue(j->GetResidue(),"CYS2");
          }
        }
      }
      ost::mol::AtomHandle cab = i->FindAtom("CAB");
      if(cab.IsValid()){
        ost::mol::AtomHandleList in_reach = handle.FindWithin(cab.GetPos(), 0.23);
        for(ost::mol::AtomHandleList::iterator j = in_reach.begin();
            j != in_reach.end(); ++j){
          if(j->GetName() == "SG"){
            if(!BondExists(cab,*j)){
              ed.Connect(cab,*j);
            }
            ed.RenameResidue(*i,"HEME");
            ed.RenameResidue(j->GetResidue(),"CYS2");
          }
        }
      }
      ost::mol::AtomHandle cac = i->FindAtom("CAC");
      if(cac.IsValid()){
        ost::mol::AtomHandleList in_reach = handle.FindWithin(cac.GetPos(),0.23);
        for(ost::mol::AtomHandleList::iterator j = in_reach.begin();
            j != in_reach.end(); ++j){
          if(j->GetName() == "SG"){
            if(!BondExists(cac,*j)){
              ed.Connect(cac,*j);
            }
            ed.RenameResidue(*i,"HEME");
            ed.RenameResidue(j->GetResidue(),"CYS2");
          }
        }
      }
    }
  }
}

void Modeller::GenerateHISHEMEBonds(ost::mol::EntityHandle& handle){

  ost::mol::ResidueHandleList res_list = handle.GetResidueList();
  ost::mol::XCSEditor ed = handle.EditXCS();

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();
      i != res_list.end(); ++i){
    if(i->GetName() == "HEM" || i->GetName() == "HEME"){
      ost::mol::AtomHandle fe = i->FindAtom("FE");
      if(fe.IsValid()){
        ost::mol::AtomHandleList in_reach = handle.FindWithin(fe.GetPos(),2.5);
        for(ost::mol::AtomHandleList::iterator j = in_reach.begin();
            j != in_reach.end(); ++j){
          if(j->GetName() == "NE2" && j->GetResidue().GetName() == "HIS"){
            if(!ost::mol::BondExists(fe,*j)){
              ed.Connect(fe,*j);
            }
            ed.RenameResidue(*i,"HEME");
            ed.RenameResidue(j->GetResidue(),"HIS1");
          }
        }
      }
    }
  }
}

void Modeller::GenerateMETHEMEBonds(ost::mol::EntityHandle& handle){

  ost::mol::ResidueHandleList res_list = handle.GetResidueList();
  ost::mol::XCSEditor ed = handle.EditXCS();

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();
      i != res_list.end(); ++i){
    if(i->GetName() == "HEM" || i->GetName() == "HEME"){
      ost::mol::AtomHandle fe = i->FindAtom("FE");
      if(fe.IsValid()){
        ost::mol::AtomHandleList in_reach = handle.FindWithin(fe.GetPos(),2.9);
        for(ost::mol::AtomHandleList::iterator j = in_reach.begin();
            j != in_reach.end(); ++j){
          if(j->GetName() == "SD" && j->GetResidue().GetName() == "MET"){
            if(!ost::mol::BondExists(fe,*j)){
              ed.Connect(fe,*j);
            }
            ed.RenameResidue(*i,"HEME");
          }
        }
      }
    }
  }
}

void Modeller::LowerPrecision(ost::mol::EntityHandle& handle){

  ost::mol::AtomHandleList atom_list = handle.GetAtomList();
  ost::mol::XCSEditor ed = handle.EditXCS(ost::mol::BUFFERED_EDIT);
  geom::Vec3 pos;
  geom::Vec3 truncated_pos;

  for(ost::mol::AtomHandleList::iterator i = atom_list.begin();
      i != atom_list.end(); ++i){
    pos = i->GetPos();
    truncated_pos[0] = Real(round(pos[0]*100))/100;
    truncated_pos[1] = Real(round(pos[1]*100))/100;
    truncated_pos[2] = Real(round(pos[2]*100))/100;
    ed.SetAtomPos(*i,truncated_pos);
  }
  ed.UpdateICS();
}

void Modeller::AssignPDBNaming(ost::mol::EntityHandle& handle){

  std::map<String,String> gromacs_to_pdb;
  gromacs_to_pdb["ARGN"] = "ARG";
  gromacs_to_pdb["ASPH"] = "ASP";
  gromacs_to_pdb["CYS2"] = "CYS";
  gromacs_to_pdb["GLUH"] = "GLU";
  gromacs_to_pdb["HISD"] = "HIS";
  gromacs_to_pdb["HISE"] = "HIS";
  gromacs_to_pdb["HISH"] = "HIS";
  gromacs_to_pdb["HIS1"] = "HIS";
  gromacs_to_pdb["LYSN"] = "LYS";
  gromacs_to_pdb["HEME"] = "HEM";

  ost::mol::XCSEditor ed = handle.EditXCS();
  ost::mol::ResidueHandleList res_list = handle.GetResidueList();
  String res_name;
  ost::mol::AtomHandle atom;

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin(); 
      i != res_list.end(); ++i){
    res_name = i->GetName();
    if(gromacs_to_pdb.find(res_name) != gromacs_to_pdb.end()){
      ed.RenameResidue(*i,gromacs_to_pdb[res_name]);
    }
    //this is nasty hardcoded stuff
    //in case of ILE gromacs seems to define wrong standard naming
    //we have to keep this wrong standard naming to remain consistent
    //with the namings in the forcefields.
    //to avoid errors in loading a pdb structure with these namings
    //we have to hardcode a correct conversion
    if(res_name == "ILE"){
      atom = i->FindAtom("CD");
      if(atom.IsValid()){
        ed.RenameAtom(atom,"CD1");
      }
      atom = i->FindAtom("HD1");
      if(atom.IsValid()){
        ed.RenameAtom(atom,"HD11");
      }      
      atom = i->FindAtom("HD2");
      if(atom.IsValid()){
        ed.RenameAtom(atom,"HD12");
      }
      atom = i->FindAtom("HD3");
      if(atom.IsValid()){
        ed.RenameAtom(atom,"HD13");
      }
    }
  }
}

void Modeller::AssignGromacsNaming(ost::mol::EntityHandle& handle){

  ost::mol::XCSEditor ed = handle.EditXCS();

  //Let's first rename as much as possible to the gromacs standard
  std::map<String,String> residue_renaming;

  //a proper solution should be found for that (assign proper protonation)
  residue_renaming["HIS"] = "HISE";
  //this is awful....


  residue_renaming["HID"] = "HISD";
  residue_renaming["HIE"] = "HISE";
  residue_renaming["HSP"] = "HISH";
  residue_renaming["LSN"] = "LYSN";
  residue_renaming["ASPP"] = "ASPH";
  residue_renaming["GLUP"] = "GLUH";
  residue_renaming["HOH"] = "SOL";
  residue_renaming["TIP3"] = "SOL";
  residue_renaming["HOH"] = "SOL";
  residue_renaming["HEME"] = "HEM";
  //following stuff are required to rename the charmm names that are different,
  //as it for example comes from CHARMM GUI
  residue_renaming["POT"] = "K";
  residue_renaming["CLA"] = "CL";
  residue_renaming["SOD"] = "NA";
  residue_renaming["CAL"] = "CA";
  residue_renaming["CES"] = "CS";
  residue_renaming["ZN2"] = "ZN";

  ost::mol::ResidueHandleList res_list = handle.GetResidueList();
  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();
      i != res_list.end(); ++i){
    if(residue_renaming.find(i->GetName()) != residue_renaming.end()){
      ed.RenameResidue(*i,residue_renaming[i->GetName()]);
    }
  }

  //this data is taken from the gromacs xlateat.dat file
  std::map<String, std::map<String, String> > renaming_to_standard;
  renaming_to_standard["ILE"] = std::map<String,String>();
  renaming_to_standard["HOH"] = std::map<String,String>();
  renaming_to_standard["HO4"] = std::map<String,String>();
  renaming_to_standard["HO5"] = std::map<String,String>();
  renaming_to_standard["HEME"] = std::map<String,String>();
  renaming_to_standard["protein"] = std::map<String,String>();

  renaming_to_standard["ILE"]["CD1"] = "CD";
  renaming_to_standard["ILE"]["HD11"] = "HD1";
  renaming_to_standard["ILE"]["HD12"] = "HD2";
  renaming_to_standard["ILE"]["HD13"] = "HD3";

  renaming_to_standard["SOL"]["O"] = "OW";
  renaming_to_standard["SOL"]["OW1"] = "OW";
  renaming_to_standard["SOL"]["OH2"] = "OW"; //specific for CHARMM-GUI
  renaming_to_standard["SOL"]["H1"] = "HW1"; //specific for CHARMM-GUI
  renaming_to_standard["SOL"]["H2"] = "HW2"; //specific for CHARMM-GUI

  renaming_to_standard["HEM"]["N_A"] = "NA";
  renaming_to_standard["HEM"]["N_B"] = "NB";
  renaming_to_standard["HEM"]["N_C"] = "NC";
  renaming_to_standard["HEM"]["N_D"] = "ND";

  renaming_to_standard["protein"]["O1"] = "O";
  renaming_to_standard["protein"]["O2"] = "OXT";
  renaming_to_standard["protein"]["OT1"] = "O";
  renaming_to_standard["protein"]["OT2"] = "OXT";
  renaming_to_standard["protein"]["OT"] = "OXT";
  renaming_to_standard["protein"]["O'"] = "O";
  renaming_to_standard["protein"]["O''"] = "OXT";
  renaming_to_standard["protein"]["OC1"] = "O";
  renaming_to_standard["protein"]["OC2"] = "OXT";
  renaming_to_standard["protein"]["HN"] = "H";
  renaming_to_standard["protein"]["HT1"] = "H1";
  renaming_to_standard["protein"]["HT2"] = "H2";
  renaming_to_standard["protein"]["HT3"] = "H3";

  //again CHARMM specific stuff
  renaming_to_standard["K"]["POT"] = "K";
  renaming_to_standard["CL"]["CLA"] = "CL";
  renaming_to_standard["NA"]["SOD"] = "NA";
  renaming_to_standard["CA"]["CAL"] = "CA";
  renaming_to_standard["CS"]["CES"] = "CS";



  //let's apply the explicit residue names definition of aboves weird data structure
  ost::mol::AtomHandleList renaming_atom_list;
  String res_name;
  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();
      i != res_list.end(); ++i){
    res_name = i->GetName();
    if(renaming_to_standard.find(res_name) == renaming_to_standard.end()) continue;
    renaming_atom_list = i->GetAtomList();
    for(ost::mol::AtomHandleList::iterator j = renaming_atom_list.begin();
        j != renaming_atom_list.end(); ++j){
      if(renaming_to_standard[res_name].find(j->GetName()) == renaming_to_standard[res_name].end()) continue;
      ed.RenameAtom(*j,renaming_to_standard[res_name][j->GetName()]);
    }
  }

  //entries with of name protein have to be applied on all peptide residues
  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();
      i != res_list.end(); ++i){
    if(!i->IsPeptideLinking()) continue;
    renaming_atom_list = i->GetAtomList();
    for(ost::mol::AtomHandleList::iterator j = renaming_atom_list.begin();
        j != renaming_atom_list.end(); ++j){
      if(renaming_to_standard["protein"].find(j->GetName()) == renaming_to_standard["protein"].end()) continue;
      ed.RenameAtom(*j,renaming_to_standard["protein"][j->GetName()]);
    }
  }
}

}}}
