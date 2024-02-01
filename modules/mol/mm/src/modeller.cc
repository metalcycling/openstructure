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

#include <ost/mol/mm/modeller.hh>

namespace{

bool IsCysteine(const ost::mol::ResidueHandle& res){
  String res_name = res.GetName();
  if(res_name[0] == 'C'){
    if(res_name == "CYS" || res_name == "CYX" || 
      res_name == "CYS2" || res_name == "CYM"){
      return true;
    }
  }
  return false;
}

bool BoundToSG(const ost::mol::AtomHandle& atom){
  ost::mol::AtomHandleList bound_atoms = atom.GetBondPartners();
  for(ost::mol::AtomHandleList::iterator i = bound_atoms.begin();
      i != bound_atoms.end(); ++i){
    if(i->GetResidue() == atom.GetResidue()) continue;
    if(i->GetName() != "SG") continue;
    if(IsCysteine(i->GetResidue())) return true;
  }
  return false;
}

}


namespace ost{ namespace mol{ namespace mm{

void Modeller::GenerateDisulfidBonds(ost::mol::EntityHandle& handle){

  ost::mol::ResidueHandleList res_list = handle.GetResidueList();
  ost::mol::XCSEditor ed = handle.EditXCS();

  for(ost::mol::ResidueHandleList::iterator i = res_list.begin();
      i != res_list.end(); ++i){

    if(!IsCysteine(*i)) continue;
    ost::mol::AtomHandle s = i->FindAtom("SG");
    if(!s.IsValid()) continue;
    if(BoundToSG(s)) continue; //it's already connected to another sulfur

    //we have now a cysteine gamma-sulfur and can look around for potential
    //disulfid bridges...
    ost::mol::AtomHandleList in_reach = handle.FindWithin(s.GetPos(),2.5);
    for(ost::mol::AtomHandleList::iterator j = in_reach.begin();
        j != in_reach.end(); ++j){
      if(j->GetResidue() != *i && j->GetName() == "SG" &&
         IsCysteine(j->GetResidue())){
        //the other atom is a gamma-sulfur as well!
        if(BoundToSG(*j)) continue; //it's already connected to another sulfur
        if(!ost::mol::BondExists(s,*j)){
          ed.Connect(s,*j);
        }
        ed.RenameResidue(*i,"CYS2");
        ed.RenameResidue(j->GetResidue(),"CYS2");
        break; //sg can only be bound to one other sg
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
