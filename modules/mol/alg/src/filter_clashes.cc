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
#include <ost/log.hh>
#include <ost/mol/mol.hh>
#include <sstream>
#include <iomanip>
#include <math.h>
#include "filter_clashes.hh"
#include <ost/units.hh>



namespace {

// helper function
String bond_string(const ost::mol::AtomView& atom1, const ost::mol::AtomHandle& atom2) {
  String atom1_str = atom1.GetName();
  String atom2_str = atom2.GetName();
  String string1,string2;
  if (atom1_str < atom2_str) {
    string1 = atom1_str;
    string2 = atom2_str;
  } else {
    string1 = atom2_str;
    string2 = atom1_str;
  }
  std::stringstream stkey;
  stkey << string1 << "-" << string2;
  return stkey.str();  
}

// helper function
String bond_string_elems(String& ele1, String ele2) {
  String string1,string2;
  if (ele1 < ele2) {
    string1 = ele1;
    string2 = ele2;
  } else {
    string1 = ele2;
    string2 = ele1;
  }
  std::stringstream stkey;
  stkey << string1 << "-" << string2;
  return stkey.str();  
}

// helper function
String angle_string(const ost::mol::AtomHandle& atom1, const ost::mol::AtomView& atom, const ost::mol::AtomHandle& atom2 ) {
  String atom1_str = atom1.GetName();
  String atom2_str = atom2.GetName();
  String string1,string2;
  if (atom1_str < atom2_str) {
    string1 = atom1_str;
    string2 = atom2_str;
  } else {
    string1 = atom2_str;
    string2 = atom1_str;
  }
  std::stringstream stkey;
  stkey << string1 << "-" << atom.GetName() << "-" << string2;
  return stkey.str();  
}


}  

namespace ost { namespace mol { namespace alg {

std::vector<ClashEvent> ClashingInfo::GetClashList() const
{
  std::vector<ClashEvent> return_list;
  std::vector<UAtomIdentifiers> seen_list;
  for (std::vector<ClashEvent>::const_iterator cl_it=clash_list_.begin();cl_it!=clash_list_.end();++cl_it) {
    UniqueAtomIdentifier atom1 = cl_it->GetFirstAtom();
    UniqueAtomIdentifier atom2 = cl_it->GetSecondAtom();
    UAtomIdentifiers check = std::make_pair(atom2,atom1);
    if (std::find(seen_list.begin(),seen_list.end(),check)==seen_list.end()) {
      return_list.push_back(*cl_it);
      seen_list.push_back(std::make_pair(atom1,atom2));
    }
  }
  return return_list;
}

void ClashingDistances::SetClashingDistance(const String& ele1,const String& ele2, Real min_distance, Real tolerance)
{
  std::stringstream stkey;
  stkey << ele1 << "--" << ele2;
  String key=stkey.str();
  min_distance_[key]=std::make_pair(min_distance,tolerance);
}

std::pair<Real,Real> ClashingDistances::GetClashingDistance(const String& ele1,const String& ele2) const
{
  std::stringstream stkey;
  if (ele1<ele2) {
    stkey << ele1 << "--" << ele2;
  } else {
    stkey << ele2 << "--" << ele1;
  }  
  String key=stkey.str();
  std::map <String,std::pair<Real,Real> >::const_iterator find_ci= min_distance_.find(key);
  if (find_ci == min_distance_.end()) {
      std::stringstream serr;
      serr << "Entry for distance " << key <<  " not found in the parameter table";
      throw Error(serr.str());
  }    
  return find_ci->second;
}

Real ClashingDistances::GetAdjustedClashingDistance(const String& ele1,const String& ele2) const
{
  std::pair <Real,Real> clash_dist = GetClashingDistance(ele1,ele2);
  return clash_dist.first-clash_dist.second;
}

void ClashingDistances::PrintAllDistances() const
{
   for (std::map <String,std::pair<Real,Real> >::const_iterator index = min_distance_.begin();index != min_distance_.end();++index) {
     std::cout << index->first << "\t" << index->second.first << "\t" << index->second.second << std::endl;
   }    
}

Real ClashingDistances::GetMaxAdjustedDistance() const
{
  Real max_adjusted_distance=0;
  for (std::map <String,std::pair<Real,Real> >::const_iterator index = min_distance_.begin();index != min_distance_.end();++index) {
    Real distance = index->second.first; 
    Real tolerance = index->second.second; 
    if ((distance-tolerance) > max_adjusted_distance) {
      max_adjusted_distance=distance-tolerance;
    }
  }  
  return max_adjusted_distance;  
}

bool ClashingDistances::IsEmpty() const
{
  if (min_distance_.size()==0) {
    return  true;
  }
  return false;
}  

void StereoChemicalParams::SetParam(const String& param, const String& residue, Real value, Real st_dev)
{
  std::pair<String,String> key = std::make_pair(param,residue);
  params_[key]=std::make_pair(value,st_dev);
}

std::pair<Real,Real> StereoChemicalParams::GetParam(const String& param,const String& residue) const
{
  std::pair<String,String> key = std::make_pair(param,residue);
  std::map<std::pair<String,String>,std::pair<Real,Real> >::const_iterator find_ci = params_.find(key);
  if (find_ci == params_.end()) {
      std::stringstream serr;
      serr << "Entry " << param <<  " for residue " << residue << " not found in the parameter table";   
      throw Error(serr.str());
  }    
  return find_ci->second;
}

bool StereoChemicalParams::ContainsParam(const String& param,const String& residue) const
{
  std::pair<String,String> key = std::make_pair(param,residue);
  std::map<std::pair<String,String>,std::pair<Real,Real> >::const_iterator find_ci = params_.find(key);
  if (find_ci == params_.end()) {
    return false;
  }    
  return true;
}

void StereoChemicalParams::PrintAllParameters() const 
{
   for (std::map <std::pair<String,String>,std::pair<Real,Real> >::const_iterator index = params_.begin();index != params_.end();++index) {
     std::cout << index->first.first << "\t" << index->first.second << "\t" << index->second.first << "\t" << index->second.second << std::endl;
   }    
};

bool StereoChemicalParams::IsEmpty() const
{
  if (params_.size()==0) {
    return  true;
  }
  return false;
}  

StereoChemicalParams FillStereoChemicalParams(const String& header, std::vector<String>& stereo_chemical_props_file, bool check)
{
  StereoChemicalParams table;
  bool found=false;
  std::vector<String>::const_iterator line_iter=stereo_chemical_props_file.begin();
  while (line_iter!=stereo_chemical_props_file.end()) {
    if ((*line_iter).length() > 1) {
      StringRef line_string_ref(line_iter->data(),(*line_iter).length());
      std::vector<StringRef> line_str_vec = line_string_ref.split();
      if (line_str_vec[0].str()==header) {
        found=true;
        line_iter++;
        while ((*line_iter)[0]!='-') {
          if ((*line_iter)[0]!='#') {
            StringRef second_line_string_ref(line_iter->data(),(*line_iter).length());
            std::vector<StringRef> second_line_str_vec = second_line_string_ref.split();
            if (second_line_str_vec.size()!=4) {
              std::cout << "The number of elements in one of the lines is wrong" << std::endl;
              return StereoChemicalParams();
            }
            StringRef item = second_line_str_vec[0];
            String res = second_line_str_vec[1].str();          
            std::pair<bool,float> parse_value = second_line_str_vec[2].to_float();
            std::pair<bool,float> parse_stddev = second_line_str_vec[3].to_float();
            Real value,stddev;
            if (parse_value.first==true) {
              value=static_cast<Real>(parse_value.second);
            } else {
              std::cout << "One of the values in the third column is not a number" << std::endl;
              return StereoChemicalParams();
            }
            if (parse_stddev.first==true) {
              stddev=static_cast<Real>(parse_stddev.second);
            } else {
              std::cout << "One of the values in the fourth column is not a number" << std::endl;
              return StereoChemicalParams();
            }
            std::vector<StringRef> split_item = item.split('-');
            String rearranged_item;
            if (split_item.size() == 2) {
              String atom1 = split_item[0].str();
              String atom2 = split_item[1].str();
              if (atom2 < atom1) {
                 std::stringstream srearr;
                 srearr << atom2 << "-" << atom1;
                 rearranged_item=srearr.str();                     
              } else {
                 rearranged_item = item.str();
              }          
            } else if (split_item.size() == 3) {
              String atom1 = split_item[0].str();
              String atom = split_item[1].str();
              String atom2 = split_item[2].str();
              if (atom2 < atom1) {
                 std::stringstream srearr;
                 srearr << atom2 << "-" << atom << "-" << atom1;
                 rearranged_item=srearr.str();                
              } else {
                 rearranged_item = item.str();
              }                
            } else {
              std::cout << "One of the strings describing the parameter has the wrong format" << std::endl;
              return StereoChemicalParams();
            }
            table.SetParam(rearranged_item,res,value,stddev);
          }
          line_iter++;
        }
      }  
    }
    line_iter++;    
  }
  if (found==false) {
    std::cout << "Could not find the relevant section in the stereo-chemical parameter file" << std::endl;
    return StereoChemicalParams();
  };
  if (check) {
    if (table.IsEmpty()) {
      std::stringstream serr;
      serr << "Error reading the " << header << " section of the stereo-chemical parameter file.";   
      throw ost::Error(serr.str());
    }
  }  
  return table;
};  

ClashingDistances FillClashingDistances(std::vector<String>& stereo_chemical_props_file, bool check)
{
  ClashingDistances table;
  bool found=false;
  std::vector<String>::const_iterator line_iter=stereo_chemical_props_file.begin();
  while (line_iter!=stereo_chemical_props_file.end()) {
    if ((*line_iter).length() > 1) {
      StringRef line_string_ref(line_iter->data(),(*line_iter).length());
      std::vector<StringRef> line_str_vec = line_string_ref.split();
      if (line_str_vec[0].str()=="Non-bonded") {
        found=true;
        line_iter++;
        while ((*line_iter)[0]!='-') {
          if ((*line_iter)[0]!='#') {
            StringRef second_line_string_ref(line_iter->data(),(*line_iter).length());
            std::vector<StringRef> second_line_str_vec = second_line_string_ref.split();
            if (second_line_str_vec.size()!=3) {
              std::cout << "The number of elements in one of the lines is wrong" << std::endl;
              return ClashingDistances();
            }
            String item = second_line_str_vec[0].str();

            std::pair<bool,float> parse_value = second_line_str_vec[1].to_float();
            std::pair<bool,float> parse_stddev = second_line_str_vec[2].to_float();
            Real value,stddev;
            if (parse_value.first==true) {
              value=static_cast<Real>(parse_value.second);
            } else {
              std::cout << "One of the distance values is not a number" << std::endl;
              return ClashingDistances();
            }
            if (parse_stddev.first==true) {
              stddev=static_cast<Real>(parse_stddev.second);
            } else {
              std::cout << "One of the tolerance values is not a number" << std::endl;
              return ClashingDistances();
            }
            StringRef itemsr(item.data(),item.length());
            std::vector<StringRef> eles = itemsr.split('-');
            if (itemsr.size() != 3) {
              std::cout << "One of the strings describing the interacting atoms has the wrong format" << std::endl;
              return ClashingDistances();
            }
            String ele1=eles[0].str();
            String ele2=eles[1].str();
            if (ele2 < ele1) {
              table.SetClashingDistance(ele2,ele1,value,stddev);
            } else {
              table.SetClashingDistance(ele1,ele2,value,stddev);
            }
          }
          line_iter++;
        }
      }  
    }
    line_iter++;    
  }
  if (found==false) {
    std::cout << "Could not find the relevant section in the stereo-chemical parameter file" << std::endl;
    return ClashingDistances();
  }
  if (check) {
    if (table.IsEmpty()) {
      std::stringstream serr;
      serr << "Error reading the Clashing section of the stereo-chemical parameter file.";   
      throw ost::Error(serr.str());
    }
  }
  return table;
}  


std::pair<EntityView,StereoChemistryInfo> CheckStereoChemistry(const EntityView& ent, const StereoChemicalParams& bond_table, const StereoChemicalParams& angle_table, Real bond_tolerance, Real angle_tolerance, bool always_remove_bb)
{
  std::vector<StereoChemicalBondViolation> bond_violation_list;
  std::vector<StereoChemicalAngleViolation> angle_violation_list;
  Real running_sum_zscore_bonds=0.0;
  Real running_sum_zscore_angles=0.0;
  int bond_count = 0;
  int bad_bond_count = 0;
  int angle_count = 1;
  int bad_angle_count = 0;
  std::map<String,Real> bond_length_sum;
  std::map<String,Real> bond_zscore_sum;
  std::map<String,int> bond_counter_sum;
  LOG_INFO("Checking stereo-chemistry")
  EntityView filtered=ent.CreateEmptyView();
  ResidueViewList residues=ent.GetResidueList();
  for (ResidueViewList::iterator i=residues.begin(), e=residues.end(); i!=e; ++i) {
    bool remove_sc=false, remove_bb=false;
    ResidueView res=*i;
    if (res.GetOneLetterCode()=='?') {
      filtered.AddResidue(res, ViewAddFlag::INCLUDE_ATOMS);
      continue;
    }  
    String residue_str = res.GetName();
    const AtomViewList& atoms=res.GetAtomList();
    for (AtomViewList::const_iterator j=atoms.begin(), e2=atoms.end(); j!=e2; ++j) {
      AtomView atom=*j;
      String ele1=atom.GetElement();
      if (ele1=="H" || ele1=="D") {
        continue;
      }
      BondHandleList bonds = atom.GetBondList();
      
      for (BondHandeList::const_iterator bi = bonds.begin();bi!=bonds.end();++bi) {
          BondHandle bond = *bi;          
          AtomHandle other_atom = bond.GetOther(atom.GetHandle());    
          if (other_atom.GetResidue()!=res.GetHandle()) {
            continue;     
          }         
          String ele2 = other_atom.GetElement();
          if (ele2=="H" || ele2=="D") {
            continue;
          }      
          if (other_atom.GetHashCode() > atom.GetHandle().GetHashCode()) {
            Real blength = bond.GetLength();
            String bond_str = bond_string(atom, other_atom);
            std::pair<Real,Real> length_stddev = bond_table.GetParam(bond_str, residue_str);
            Real ref_length = length_stddev.first;
            Real ref_stddev = length_stddev.second;           
            Real min_length = ref_length - bond_tolerance*ref_stddev;
            Real max_length = ref_length + bond_tolerance*ref_stddev;
            Real zscore = (blength - ref_length)/ref_stddev;
            if (blength < min_length || blength > max_length) {
              LOG_INFO("BOND:" << " " << res.GetChain() << " "
                       << res.GetName() << " " << res.GetNumber() << " "
                       << bond_str << " " << min_length << " " << max_length
                       << " " << blength << " " << zscore << " " << "FAIL")
              bad_bond_count++;
              UniqueAtomIdentifier atom_ui(atom.GetResidue().GetChain().GetName(),
                                           atom.GetResidue().GetNumber(),
                                           atom.GetResidue().GetName(), atom.GetName());
              UniqueAtomIdentifier other_atom_ui(other_atom.GetResidue().GetChain().GetName(),
                                                 other_atom.GetResidue().GetNumber(),
                                                 other_atom.GetResidue().GetName(),
                                                 other_atom.GetName());
              StereoChemicalBondViolation bond_v(atom_ui, other_atom_ui, blength,
                                                 std::make_pair(min_length, max_length));
              bond_violation_list.push_back(bond_v);
              remove_sc=true;
              if (always_remove_bb==true) {
                remove_bb=true;
              } else {
                // we need to check both atom names since the order is random!
                // -> for angles and clashes this is not needed
                String name1 = atom.GetName();
                String name2 = other_atom.GetName();
                if (name1=="CA" || name1=="N" || name1=="O" || name1=="C" ||
                    name2=="CA" || name2=="N" || name2=="O" || name2=="C") {
                  remove_bb=true;
                }
              }
            } else {
              LOG_VERBOSE("BOND:" << " " << res.GetChain() << " "
                          << res.GetName() << " " << res.GetNumber() << " "
                          << bond_str << " " << min_length << " " << max_length
                          << " " << blength << " " << zscore << " " << "PASS")
            }
            bond_count++;
            running_sum_zscore_bonds+=zscore;
            String bond_elems=bond_string_elems(ele1,ele2);
            std::map<String,Real>::const_iterator find_be = bond_length_sum.find(bond_elems);  
            if (find_be==bond_length_sum.end()) {
                bond_length_sum[bond_elems]=blength;
                bond_zscore_sum[bond_elems]=zscore;
                bond_counter_sum[bond_elems]=1;
            } else {
                bond_length_sum[bond_elems]+=blength;
                bond_zscore_sum[bond_elems]+=zscore;
                bond_counter_sum[bond_elems]+=1;
            }
          }  
      }
      
      for (BondHandeList::const_iterator bond_iter1=bonds.begin(); bond_iter1!=bonds.end(); ++bond_iter1) {
        BondHandle bond1=*bond_iter1;
        AtomHandle atom1= bond1.GetOther(atom.GetHandle());
        String ele_atom1=atom1.GetElement();
        if (ele_atom1=="H" || ele_atom1=="D") {
          continue;
        }
        if (atom1.GetResidue()!=res.GetHandle()) {
          continue;        
        }        
        for (BondHandeList::const_iterator bond_iter2=bonds.begin(); bond_iter2!=bonds.end(); ++bond_iter2) {
          BondHandle bond2=*bond_iter2;
          AtomHandle atom2 = bond2.GetOther(atom.GetHandle());
          String ele_atom2=atom2.GetElement();
          if (ele_atom2=="H" || ele_atom2=="D") {
            continue;
          }
          if (atom2.GetResidue()!=res.GetHandle()) {
            continue;        
          }
          if (atom1.GetHashCode() > atom2.GetHashCode()) {
            Real awidth;
            if (atom1.GetName()<atom2.GetName()) {
              awidth = ent.GetAngle(atom1,atom.GetHandle(),atom2);
            } else {
              awidth = ent.GetAngle(atom2,atom.GetHandle(),atom1);
            }    
            awidth/=(ost::Units::deg);
            String angle_str = angle_string(atom1,atom,atom2);
            std::pair<Real,Real> width_stddev = angle_table.GetParam(angle_str,residue_str);
            Real ref_width = width_stddev.first;  
            Real ref_stddev = width_stddev.second;           
            Real min_width = ref_width - angle_tolerance*ref_stddev;
            Real max_width = ref_width + angle_tolerance*ref_stddev;
            Real zscore = (awidth - ref_width)/ref_stddev;
            if (awidth < min_width || awidth > max_width) {
              LOG_INFO("ANGLE:" << " " << res.GetChain() << " "
                       << res.GetName() << " " << res.GetNumber() << " "
                       << angle_str << " " << min_width << " " << max_width
                       << " " << awidth << " " << zscore << " " << "FAIL")
              bad_angle_count++;
              UniqueAtomIdentifier atom1_ui(atom1.GetResidue().GetChain().GetName(),
                                            atom1.GetResidue().GetNumber(),
                                            atom1.GetResidue().GetName(),
                                            atom1.GetName());
              UniqueAtomIdentifier atom_ui(atom.GetResidue().GetChain().GetName(),
                                           atom.GetResidue().GetNumber(),
                                           atom.GetResidue().GetName(),
                                           atom.GetName());
              UniqueAtomIdentifier atom2_ui(atom2.GetResidue().GetChain().GetName(),
                                            atom2.GetResidue().GetNumber(),
                                            atom2.GetResidue().GetName(),
                                            atom2.GetName());
              StereoChemicalAngleViolation angle_v(atom1_ui,atom_ui,atom2_ui,awidth,
                                                   std::make_pair(min_width,max_width));
              angle_violation_list.push_back(angle_v);
              remove_sc=true;
              if (always_remove_bb==true) {
                remove_bb=true;
              }
              String name=atom.GetName();
              if (name=="CA" || name=="N" || name=="O" || name=="C") {
                remove_bb=true;
              }
            } else {
              LOG_VERBOSE("ANGLE:" << " " << res.GetChain() << " "
                          << res.GetName() << " " << res.GetNumber() << " "
                          << angle_str << " " << min_width << " " << max_width
                          << " " << awidth << " " << zscore << " " << "PASS")
            }
            angle_count++;
            running_sum_zscore_angles+=zscore;  
          }  
        }
      }         
    }
    
    if (remove_bb) {
      LOG_INFO("ACTION: removing whole residue " << res);
      res.SetBoolProp("stereo_chemical_violation_backbone",true);
      continue;
    }
    if (remove_sc) {
      LOG_INFO("ACTION: removing sidechain of residue " << res);
      for (AtomViewList::const_iterator 
           j=atoms.begin(), e2=atoms.end(); j!=e2; ++j) {
       AtomView atom=*j;
       String name=atom.GetName();
       if (name=="CA" || name=="N" || name=="O" || name=="C") {
         filtered.AddAtom(atom);
       }
      }
      res.SetBoolProp("stereo_chemical_violation_sidechain",true);
      continue;
    }
    filtered.AddResidue(res, ViewAddFlag::INCLUDE_ATOMS);
  }
  Real avg_zscore_bonds = running_sum_zscore_bonds/static_cast<float>(bond_count);
  Real avg_zscore_angles = running_sum_zscore_angles/static_cast<float>(angle_count);
  std::map<String,BondLengthInfo> avg_bond_length_info;
  for (std::map<String,Real>::const_iterator bls_it=bond_length_sum.begin();bls_it!=bond_length_sum.end();++bls_it) {
    String key = (*bls_it).first;
    int counter=bond_counter_sum[key];
    Real sum_bond_length=(*bls_it).second;
    Real sum_bond_zscore=bond_zscore_sum[key];
    Real avg_length=sum_bond_length/static_cast<Real>(counter);
    Real avg_zscore=sum_bond_zscore/static_cast<Real>(counter);
    BondLengthInfo bond_info(avg_length,avg_zscore,counter);
    avg_bond_length_info[key]=bond_info;
  }
  StereoChemistryInfo info(avg_zscore_bonds, bad_bond_count, bond_count,avg_zscore_angles,
                           bad_angle_count, angle_count,avg_bond_length_info,
                           bond_violation_list,angle_violation_list);
  filtered.AddAllInclusiveBonds();
  return std::make_pair(filtered,info);
}


std::pair<EntityView,StereoChemistryInfo> CheckStereoChemistry(const EntityHandle& ent, const StereoChemicalParams& bond_table, const StereoChemicalParams& angle_table, Real bond_tolerance, Real angle_tolerance, bool always_remove_bb)
{
  return CheckStereoChemistry(ent.CreateFullView(), bond_table, angle_table, bond_tolerance, angle_tolerance, always_remove_bb);
}


std::pair<EntityView,ClashingInfo> FilterClashes(const EntityView& ent, const ClashingDistances& min_distances, bool always_remove_bb)
{
  std::vector<ClashEvent> clash_list;
  int distance_count = 0;
  int bad_distance_count = 0;
  Real average_offset_sum = 0.0;
  LOG_INFO("Filtering non-bonded clashes")
  EntityView filtered=ent.CreateEmptyView();
  ResidueViewList residues=ent.GetResidueList();
  for (ResidueViewList::iterator 
       i=residues.begin(), e=residues.end(); i!=e; ++i) {
    bool remove_sc=false, remove_bb=false;
    ResidueView res=*i;
    if (res.GetOneLetterCode()=='?') {
      filtered.AddResidue(res, ViewAddFlag::INCLUDE_ATOMS);
      continue;
    }  
    const AtomViewList& atoms=res.GetAtomList();
    for (AtomViewList::const_iterator 
         j=atoms.begin(), e2=atoms.end(); j!=e2; ++j) {
      AtomView atom=*j;

      String ele1=atom.GetElement();
      if (ele1=="H" || ele1=="D") {
        continue;
      }
      AtomViewList within=ent.FindWithin(atom.GetPos(),min_distances.GetMaxAdjustedDistance());
      for (AtomViewList::iterator 
           k=within.begin(), e3=within.end(); k!=e3; ++k) {
        AtomView atom2=*k;
        if (atom2==atom) {
          continue;
        }
        String ele2=atom2.GetElement();
        if (ele2=="H" || ele2=="D") {
          continue;
        }

        // In theory, this should also trigger for disulfide bonds, but 
        // since we don't detect disulfides correctly, we can't count on that 
        // and we instead allow S-S distances down to 1.8.       
        if (atom.GetHandle().FindBondToAtom(atom2.GetHandle()).IsValid()) {
          continue;
        }

        Real d=geom::Length2(atom.GetPos()-atom2.GetPos());
        std::pair <Real,Real> distance_tolerance=min_distances.GetClashingDistance(ele1, ele2);
        Real distance=distance_tolerance.first;
        Real tolerance=distance_tolerance.second;
        Real threshold=distance-tolerance;

        if (d<threshold*threshold) {
          LOG_INFO("CLASH:" << " " << atom.GetResidue().GetChain() << " "
                   << atom.GetResidue().GetName() << " "
                   << atom.GetResidue().GetNumber() << " " << atom.GetName()
                   << " " << atom2.GetResidue().GetChain() << " "
                   << atom2.GetResidue().GetName() << " "
                   << atom2.GetResidue().GetNumber() << " " << atom2.GetName()
                   << " " << threshold << " " << sqrt(d) << " "
                   << sqrt(d)-threshold << " " << "FAIL")
          bad_distance_count++; 
          UniqueAtomIdentifier atom_ui(atom.GetResidue().GetChain().GetName(),atom.GetResidue().GetNumber(),atom.GetResidue().GetName(),atom.GetName());
          UniqueAtomIdentifier atom2_ui(atom2.GetResidue().GetChain().GetName(),atom2.GetResidue().GetNumber(),atom2.GetResidue().GetName(),atom2.GetName());
          ClashEvent clash_event(atom_ui,atom2_ui,sqrt(d),threshold);
          clash_list.push_back(clash_event);
          average_offset_sum+=sqrt(d)-threshold;
          remove_sc=true;
          if (always_remove_bb==true) {
            remove_bb=true;
          }
          String name=atom.GetName();
          if (name=="CA" || name=="N" || name=="O" || name=="C") {
            remove_bb=true;
          }
        } else {
          LOG_VERBOSE("CLASH:" << " " << atom.GetResidue().GetChain() << " "
                      << atom.GetResidue().GetName() << " "
                      << atom.GetResidue().GetNumber() << " " << atom.GetName()
                      << " " << atom2.GetResidue().GetChain() << " "
                      << atom2.GetResidue().GetName() << " "
                      << atom2.GetResidue().GetNumber() << " " << atom2.GetName()
                      << " " << threshold << " " << sqrt(d) << " "
                      << sqrt(d)-threshold << " " << "PASS")
        }
        distance_count++;
      }
    }
    
    if (remove_bb) {
      LOG_VERBOSE("ACTION: removing whole residue " << res);
      res.SetBoolProp("steric_clash_backbone", true);
      continue;
    }
    if (remove_sc) {
      LOG_VERBOSE("ACTION: removing sidechain of residue " << res);
      for (AtomViewList::const_iterator 
           j=atoms.begin(), e2=atoms.end(); j!=e2; ++j) {
       AtomView atom=*j;
       String name=atom.GetName();
       if (name=="CA" || name=="N" || name=="O" || name=="C") {
         filtered.AddAtom(atom);
       }
      }
      res.SetBoolProp("steric_clash_sidechain", true);
      continue;
    }
    filtered.AddResidue(res, ViewAddFlag::INCLUDE_ATOMS);
  }
  Real average_offset = 0;
  if (bad_distance_count!=0) {
    average_offset = average_offset_sum / static_cast<Real>(bad_distance_count);
  }
  ClashingInfo info(bad_distance_count,average_offset,clash_list);
  filtered.AddAllInclusiveBonds();
  return std::make_pair(filtered,info);
}


std::pair<EntityView,ClashingInfo> FilterClashes(const EntityHandle& ent,
                                                 const ClashingDistances& min_distances, bool always_remove_bb)
{
  return FilterClashes(ent.CreateFullView(), min_distances, always_remove_bb);
}


}}}
