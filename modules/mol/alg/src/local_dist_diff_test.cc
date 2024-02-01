#include <iomanip>
#include <sstream>
#include <ost/log.hh>
#include <ost/mol/mol.hh>
#include <ost/platform.hh>
#include "local_dist_diff_test.hh"
#include <boost/concept_check.hpp>
#include <boost/filesystem/convenience.hpp>
#include <ost/mol/alg/consistency_checks.hh>

namespace ost { namespace mol { namespace alg {

namespace {

// helper_function
String vector_to_string(std::vector<Real> vec) {
  std::ostringstream str;
  for (unsigned int i = 0; i < vec.size(); ++i) {
      str << vec[i];
      if (i+1 != vec.size()) {
          str << ", ";
      }
  }
  return str.str();
}

// helper function
bool within_tolerance(Real mdl_dist, const std::pair<Real,Real>& values, Real tol)
{
  // this considers min & max distance observed in all reference structures
  // -> for single reference case values.first and values.second are the same
  return (values.first-tol)<=mdl_dist && (values.second+tol)>=mdl_dist;
}      

// helper function
std::pair<long int, long int> calc_overlap1(const ResidueRDMap& res_distance_list, const ResNum& rnum,
                                    ChainView mdl_chain, int sequence_separation,
                                    std::vector<Real>& tol_list, bool only_fixed, 
                                    bool swap,std::vector<std::pair<long int, long int> >& overlap_list, bool log )
{
  std::pair<long int, long int> overlap(0, 0);
  ResidueView mdl_res=mdl_chain.FindResidue(rnum); 
  for (ResidueRDMap::const_iterator ai=res_distance_list.begin(); ai!=res_distance_list.end(); ++ai) {
    const UAtomIdentifiers& uais = ai->first; 
    const std::pair <Real,Real>& values = ai->second;
    const UniqueAtomIdentifier& first_atom=uais.first;
    const UniqueAtomIdentifier& second_atom=uais.second;
    String name=swap ? SwappedName(first_atom.GetAtomName()) : first_atom.GetAtomName();
    AtomView av1=mdl_res ? mdl_res.FindAtom(name) : AtomView();
 
    if (only_fixed) {
       if (std::abs(first_atom.GetResNum().GetNum()-second_atom.GetResNum().GetNum())<=sequence_separation) {
          continue;
        }
      if (Swappable(second_atom.GetResidueName(), second_atom.GetAtomName())) {
          continue;
      }
    }
    if (!only_fixed) {   
      if (first_atom.GetResNum().GetNum()<=(second_atom.GetResNum().GetNum()+sequence_separation)) {
        continue;
      }    
    }
    ResidueView rv2=mdl_chain.FindResidue(second_atom.GetResNum());
    overlap.second+=tol_list.size();
    int rindex2=0;
    int rindex1=mdl_res ? mdl_res.GetIndex() : -1;
    if (!only_fixed && rindex1!=-1)
      overlap_list[rindex1].second+=tol_list.size();
    if (!rv2) {
      continue;
    }
    rindex2=rv2.GetIndex();
    if (!only_fixed)
      overlap_list[rindex2].second+=tol_list.size();
    AtomView av2=rv2.FindAtom(second_atom.GetAtomName());
    if (!(av1 && av2)) {
     continue;
    }
    Real mdl_dist=geom::Length(av1.GetPos()-av2.GetPos());
    std::vector<Real>::const_reverse_iterator rend_it=tol_list.rend();
    for (std::vector<Real>::const_reverse_iterator tol_list_it=tol_list.rbegin();tol_list_it!=rend_it;++tol_list_it) {
      Real tol = * tol_list_it; 
      if (within_tolerance(mdl_dist,values,tol)) {
        if (log) {
          LOG_VERBOSE("LDDT:" << " " << av1.GetResidue().GetChain() << " "
                      << av1.GetResidue().GetName() << " "
                      << av1.GetResidue().GetNumber() << " " << av1.GetName()
                      << " " << av2.GetResidue().GetChain() << " "
                      << av2.GetResidue().GetName() << " "
                      << av2.GetResidue().GetNumber() << " " << av2.GetName()
                      << " " << mdl_dist << " " << values.first << " "
                      << values.second << " " << tol << " " << "PASS")
        }
        overlap.first+=1;
        if (!only_fixed) {
          overlap_list[rindex1].first+=1;
          overlap_list[rindex2].first+=1;
        }
      } else {
        if (log) {
          LOG_VERBOSE("LDDT:" << " " << av1.GetResidue().GetChain() << " "
                      << av1.GetResidue().GetName() << " "
                      << av1.GetResidue().GetNumber() << " " << av1.GetName()
                      << " " << av2.GetResidue().GetChain() << " "
                      << av2.GetResidue().GetName() << " "
                      << av2.GetResidue().GetNumber() << " " << av2.GetName()
                      << " " << mdl_dist << " " << values.first << " "
                      << values.second << " " << tol << " " << "FAIL");
        }
        break;  
      }
    }    
  }
  return overlap;
}  

// helper function used by the alignment-based Local Distance Difference Test
std::pair<Real, Real> calc_overlap2(const seq::ConstSequenceHandle& ref_seq,
                                    const seq::ConstSequenceHandle& mdl_seq,
                                    int pos, Real tol, Real max_dist, 
                                    bool only_fixed, bool swap)
{
  std::pair<Real, Real> overlap(0.0, 0.0);
  EntityView ref=ref_seq.GetAttachedView();
  ResidueView ref_res=ref_seq.GetResidue(pos);
  if (!ref_res.IsValid()) {
    return std::pair<Real,Real>(0.0, 0.0);
  }
  AtomViewList ref_atoms=ref_res.GetAtomList();
  ResidueView mdl_res=mdl_seq.GetResidue(pos);
  AtomViewList within;
  if (max_dist<0) {
    within=ref.GetAtomList();
  }  
  for (AtomViewList::iterator ai=ref_atoms.begin(),
       ae=ref_atoms.end(); ai!=ae; ++ai) {
    if (ai->GetElement()=="H") { continue; }
    String name=swap ? SwappedName(ai->GetName()) : ai->GetName();
    AtomView av1=mdl_res ? mdl_res.FindAtom(name) : AtomView();
    if (max_dist>=0){ 
      within=ref.FindWithin(ai->GetPos(), max_dist);
    }
    for (AtomViewList::iterator aj=within.begin(),
         ae2=within.end(); aj!=ae2; ++aj) {
      if (aj->GetElement()=="H" || 
          aj->GetResidue().GetChain()!=ai->GetResidue().GetChain()) { 
          continue; 
      }
      if (only_fixed) {
        if (aj->GetResidue().GetNumber()==ref_res.GetNumber()) {
          continue;
        }
        if (Swappable(aj->GetResidue().GetName(), aj->GetName())) {
          continue;
        }
        overlap.second+=1.0;
        // map from residue index to position in alignment
        try {
          int aln_pos=ref_seq.GetPos(aj->GetResidue().GetIndex());
          ResidueView r2=mdl_seq.GetResidue(aln_pos);
          if (!r2.IsValid()) {
            continue;
          }
          AtomView av2=r2.FindAtom(aj->GetName());
          if (!(av1 && av2)) {
            continue;
          }
          Real mdl_dist=geom::Length(av1.GetPos()-av2.GetPos());
          Real ref_dist=geom::Length(ai->GetPos()-aj->GetPos());
          if (std::abs(mdl_dist-ref_dist)<tol) {
            overlap.first+=1;
          }
        } catch(...) { }
        continue;
      } else {
        if (aj->GetResidue().GetNumber()>ref_res.GetNumber()) {
          overlap.second+=1.0;    

          try {
           int aln_pos=ref_seq.GetPos(aj->GetResidue().GetIndex());
            ResidueView r2=mdl_seq.GetResidue(aln_pos);
            if (!r2.IsValid()) {
              continue;
            }
            AtomView av2=r2.FindAtom(aj->GetName());
            if (!(av1 && av2)) {
              continue;
            }
            Real mdl_dist=geom::Length(av1.GetPos()-av2.GetPos());
            Real ref_dist=geom::Length(ai->GetPos()-aj->GetPos());
            if (std::abs(mdl_dist-ref_dist)<tol) {
              overlap.first+=1;
            }
          } catch (...) { }
        }  
      }
    }      
  }
  return overlap;
}

// for each residue with multiple possible nomenclature conventions, checks which choice (switched or not) 
// of atom nomenclature gives the highest lddt score then changes the naming convention of the input 
// entity view accordingly
void check_and_swap(const GlobalRDMap& glob_dist_list, const EntityView& mdl, std::vector<Real> cutoff_list, int sequence_separation,  std::vector<std::pair<long int, long int> > overlap_list)
{
  ChainView mdl_chain=mdl.GetChainList()[0]; 
  XCSEditor edi=mdl.GetHandle().EditXCS(BUFFERED_EDIT);
  for (GlobalRDMap::const_iterator i=glob_dist_list.begin(); i!=glob_dist_list.end(); ++i) {  
    ResNum rnum =  i->first;
    if (i->second.size()==0) {
      continue;
    }
    ResidueView mdl_res=mdl_chain.FindResidue(rnum); 
    if (!mdl_res) {   
      continue;
    }
    String rname = mdl_res.GetName();
    if (!(rname=="GLU" || rname=="ASP" || rname=="VAL" || rname=="TYR" ||
          rname=="PHE" || rname=="LEU" || rname=="ARG")) {
      continue;
    }
    std::pair<long int, long int> ov1=calc_overlap1(i->second, rnum,mdl_chain, sequence_separation,
                                          cutoff_list, true, 
                                          false, overlap_list,false);
 
    std::pair<long int, long int> ov2=calc_overlap1(i->second, rnum, mdl_chain, sequence_separation,
                                          cutoff_list, true, 
                                          true, overlap_list,false); 

    if (static_cast<Real>(ov1.first)/ov1.second<
        static_cast<Real>(ov2.first)/ov2.second) {
      AtomViewList atoms=mdl_res.GetAtomList();
      for (AtomViewList::iterator j=atoms.begin(), e2=atoms.end(); j!=e2; ++j) {
        if (Swappable(rname, j->GetName())) {
          edi.RenameAtom(j->GetHandle(), SwappedName(j->GetName()));
        }
      }
    } 
  }
}

// helper function to update existence map entries for multiple reference input structures
void update_existence_map (ExistenceMap& ex_map, const EntityView& ev, int ref_counter)
{
  AtomViewList ev_atom=ev.GetAtomList();
  for (AtomViewList::iterator ev_atom_it=ev_atom.begin(); ev_atom_it!=ev_atom.end();++ev_atom_it) {
    UniqueAtomIdentifier uai (ev_atom_it->GetResidue().GetChain().GetName(),ev_atom_it->GetResidue().GetNumber(),ev_atom_it->GetResidue().GetName(),ev_atom_it->GetName());
    ex_map[uai] |=1 << ref_counter;
  }
}

// helper function for super-fast lookup of atom existence in multiple reference input structures
int in_existence_map(const ExistenceMap& ex_map, const UniqueAtomIdentifier& uai) 
{
  ExistenceMap::const_iterator find_uai_ci = ex_map.find(uai);
  return find_uai_ci!=ex_map.end() ? find_uai_ci->second : 0;
}    

// merges distance lists from multiple reference structures. The login is described in the code
void merge_distance_lists(GlobalRDMap& ref_dist_map, const GlobalRDMap& new_dist_map, ExistenceMap& ex_map, const EntityView& ref,int ref_counter)
{
  // iterate over the residues in the ref_dist_map   
  for (GlobalRDMap::iterator ref_dist_map_it=ref_dist_map.begin();ref_dist_map_it!=ref_dist_map.end();++ref_dist_map_it) {
    ResNum ref_resnum = ref_dist_map_it->first;
    GlobalRDMap::const_iterator find_new_res_ci = new_dist_map.find(ref_resnum);
    //if the residue is found in new_dist_map,  
    if (find_new_res_ci != new_dist_map.end()) {
      //iterate over the the reference distances in the ResidueDistanceMap  
      //It's on purpose that we don't increase the loop variable inside
      //the for statement! This is required to make iteration work when erasing
      //an element from the map while iterating over it.
      for (ResidueRDMap::iterator
           ref_res_map_it = ref_dist_map_it->second.begin();
           ref_res_map_it!=ref_dist_map_it->second.end();) {
        const UAtomIdentifiers& ref_rd = ref_res_map_it->first;
        bool erased=false;
        std::pair<Real,Real>& ref_minmax = ref_res_map_it->second;
        ResidueRDMap::const_iterator find_new_rd_ci = find_new_res_ci->second.find(ref_rd);
        // if you find the distance in the residue new, udate min and max  
        if (find_new_rd_ci != find_new_res_ci->second.end()) {
          if (find_new_rd_ci->second.first < ref_minmax.first) {
            ref_minmax.first=find_new_rd_ci->second.first;
          } else if (find_new_rd_ci->second.second > ref_minmax.second) {
            ref_minmax.second=find_new_rd_ci->second.second;
          }
        } else {              
          // if you don't find it in the residue new, check that it is not missing because it is too long  
          UniqueAtomIdentifier first_atom_to_find = ref_rd.first;
          UniqueAtomIdentifier second_atom_to_find = ref_rd.second;  
          // if both atoms are there, remove the distance from the ref_dist_map,
          if ((ref.FindAtom(first_atom_to_find.GetChainName(),first_atom_to_find.GetResNum(),first_atom_to_find.GetAtomName()).IsValid() &&
               ref.FindAtom(second_atom_to_find.GetChainName(),second_atom_to_find.GetResNum(),second_atom_to_find.GetAtomName()).IsValid()) ) {
            erased=true;
            ref_dist_map_it->second.erase(ref_res_map_it++);
          }
        }
        if (!erased) {
          ++ref_res_map_it;
        }
      }
      // now iterate over the the new  reference distances in residue new
      for (ResidueRDMap::const_iterator new_res_map_it = find_new_res_ci->second.begin(); new_res_map_it!=find_new_res_ci->second.end();++new_res_map_it) {
        UAtomIdentifiers new_rd = new_res_map_it->first;
        std::pair<Real,Real> new_minmax = new_res_map_it->second;
        ResidueRDMap::const_iterator find_ref_rd_ci = ref_dist_map_it->second.find(new_rd);
        // if the distance is found in the residue ref, 
        // it has been taken care of before. If not
        if (find_ref_rd_ci==ref_dist_map_it->second.end()) {
          UniqueAtomIdentifier first_atom_to_find = new_rd.first;
          UniqueAtomIdentifier second_atom_to_find = new_rd.second;
          // check that there isn't a structure already processed where both atoms are in
          // if there is none, add the distance to the residue ref map 
          if (!(in_existence_map(ex_map,first_atom_to_find) & in_existence_map(ex_map,second_atom_to_find))) {
            ref_dist_map_it->second[new_rd]= new_minmax; 
          }
        }      
      } 
    }
  // if the residue was not found in the new list, it means that it is 
  // absent in the new structure, no new information
  }
  // now iterate over the residues in the new_list
  for (GlobalRDMap::const_iterator new_dist_map_it=new_dist_map.begin();new_dist_map_it!=new_dist_map.end();++new_dist_map_it) {
    ResNum new_resnum = new_dist_map_it->first;
    GlobalRDMap::const_iterator find_ref_res_ci = ref_dist_map.find(new_resnum);
    // if the residue is found in new_dist_map, it has been taken care before,
    // if not, add it to the res_dist_map:
      if (find_ref_res_ci == ref_dist_map.end()) {
      ref_dist_map[new_resnum] = new_dist_map_it->second;           
    }
  }
}

}

// Computes coverage
std::pair<int,int> ComputeCoverage(const EntityView& v,const GlobalRDMap& glob_dist_list)
{
  int second=0;
  int first=0;
  if (v.GetResidueList().size()==0) {
    if (glob_dist_list.size()==0) {
      return std::make_pair(0,-1);
    } else {    
      return std::make_pair(0,glob_dist_list.size());
    }  
  }
  ChainView vchain=v.GetChainList()[0];
  for (GlobalRDMap::const_iterator i=glob_dist_list.begin();i!=glob_dist_list.end();++i)
  {
    ResNum rnum = (*i).first;
    second++;
    if (vchain.FindResidue(rnum)) {
      first++;
    }
  }
  return std::make_pair(first,second);
}

bool IsResnumInGlobalRDMap(const ResNum& resnum, const GlobalRDMap& glob_dist_list)
{
  for (GlobalRDMap::const_iterator i=glob_dist_list.begin(), e=glob_dist_list.end(); i!=e; ++i) {
    ResNum rn = i->first;
    if (rn==resnum) {
      return true;
    }
  }
  return false;
}


// helper function
bool IsStandardResidue(String rn)
{
  String upper_rn=rn;
  std::transform(rn.begin(),rn.end(),rn.begin(),toupper);
  if (upper_rn == "ALA" ||    
  upper_rn == "ARG" ||
  upper_rn == "ASN" ||
  upper_rn == "ASP" ||
  upper_rn == "GLN" ||
  upper_rn == "GLU" ||
  upper_rn == "LYS" ||
  upper_rn == "SER" ||
  upper_rn == "CYS" ||      
  upper_rn == "TYR" ||  
  upper_rn == "TRP" || 
  upper_rn == "THR" ||
  upper_rn == "VAL" ||
  upper_rn == "ILE" ||
  upper_rn == "MET" ||
  upper_rn == "LEU" ||  
  upper_rn == "GLY" ||  
  upper_rn == "PRO" ||  
  upper_rn == "HIS" ||  
  upper_rn == "PHE") {
    return true;
  }  
  return false;
}

StereoChemicalProps::StereoChemicalProps():
    is_valid(false) {}

StereoChemicalProps::StereoChemicalProps(
  ost::mol::alg::StereoChemicalParams& init_bond_table,
  ost::mol::alg::StereoChemicalParams& init_angle_table,
  ost::mol::alg::ClashingDistances& init_nonbonded_table):
    is_valid(true),
    bond_table(init_bond_table),
    angle_table(init_angle_table),
    nonbonded_table(init_nonbonded_table) {}

lDDTSettings::lDDTSettings(): radius(15.0), 
                              sequence_separation(0),
                              label("locallddt") {
    cutoffs.push_back(0.5);
    cutoffs.push_back(1.0);
    cutoffs.push_back(2.0);
    cutoffs.push_back(4.0);
  }

lDDTSettings::lDDTSettings(Real init_radius, 
                           int init_sequence_separation,
                           std::vector<Real>& init_cutoffs,
                           String init_label):
                    radius(init_radius), 
                    sequence_separation(init_sequence_separation),
                    cutoffs(init_cutoffs),
                    label(init_label) {}

std::string lDDTSettings::ToString() {
  std::ostringstream rep;
  rep << "Inclusion Radius: " << radius << "\n";
  rep << "Sequence separation: " << sequence_separation << "\n";
  rep << "Cutoffs: " << vector_to_string(cutoffs) << "\n";
  rep << "Residue properties label: " << label << "\n";

  return rep.str();
}

void lDDTSettings::PrintParameters() {
  std::cout << ToString();
}

// Default constructor is neccessary for boost exposure
lDDTLocalScore::lDDTLocalScore(): cname(""),
                                  rname(""),
                                  rnum(-1),
                                  is_assessed(""),
                                  quality_problems(""),
                                  local_lddt(-1.0),
                                  conserved_dist(-1),
                                  total_dist(-1) {}

lDDTLocalScore::lDDTLocalScore(String init_cname,
                               String init_rname,
                               int init_rnum,
                               String init_is_assessed,
                               String init_quality_problems,
                               Real init_local_lddt,
                               int init_conserved_dist,
                               int init_total_dist): 
                                       cname(init_cname),
                                       rname(init_rname),
                                       rnum(init_rnum),
                                       is_assessed(init_is_assessed),
                                       quality_problems(init_quality_problems),
                                       local_lddt(init_local_lddt),
                                       conserved_dist(init_conserved_dist),
                                       total_dist(init_total_dist) {}

String lDDTLocalScore::ToString(bool structural_checks) const {
  std::stringstream stkeylddt;
  std::stringstream outstr;
  stkeylddt << "(" << conserved_dist << "/" << total_dist << ")";
  String dist_string = stkeylddt.str();
  if (structural_checks) {
    outstr << cname << "\t" << rname << "\t" << rnum << '\t' << is_assessed  << '\t' << quality_problems << '\t' << local_lddt << "\t" << dist_string;
  } else {
    outstr << cname << "\t" << rname << "\t" << rnum << '\t' << is_assessed  << '\t' << local_lddt << "\t" << dist_string;
  }
  return outstr.str();
}

String lDDTLocalScore::Repr() const {
  std::stringstream outstr;
  outstr << "<lDDTLocalScore " << cname << "." << rname << "." << rnum << ">";
  return outstr.str();
}

String lDDTLocalScore::GetHeader(bool structural_checks, int cutoffs_length) {
  std::stringstream outstr;
  if (structural_checks) {
    outstr << "Chain\tResName\tResNum\tAsses.\tQ.Prob.\tScore\t(Conserved/Total, over " << cutoffs_length << " thresholds)";
  } else {
    outstr << "Chain\tResName\tResNum\tAsses.\tScore\t(Conserved/Total, over " << cutoffs_length << " thresholds)";
  }
  return outstr.str();
}

lDDTScorer::lDDTScorer(std::vector<EntityView>& init_references,
                       ost::mol::EntityView& init_model,
                       lDDTSettings& init_settings):
  settings(init_settings),
  model_view(init_model),
  references_view(init_references)
  {
    _score_calculated = false;
    _score_valid = false;
    _has_local_scores = false;
    _num_cons_con = -1;
    _num_tot_con = -1;
    _global_score = -1.0;
    CleanlDDTReferences(references_view);
    _PrepareGlobalRDMap();  
  }

Real lDDTScorer::GetGlobalScore(){
  if (!_score_calculated) {
    _ComputelDDT();
  }
  return _global_score;
}

int lDDTScorer::GetNumConservedContacts(){
  if (!_score_calculated) {
    _ComputelDDT();
  }
  return _num_cons_con;
}

int lDDTScorer::GetNumTotalContacts(){
  if (!_score_calculated) {
    _ComputelDDT();
  }
  return _num_tot_con;
}

std::vector<lDDTLocalScore> lDDTScorer::GetLocalScores(){
  if (!_has_local_scores) {
    _GetLocallDDT();
  }
  return _local_scores;
}

void lDDTScorer::PrintPerResidueStats(){
  if (!_has_local_scores) {
    _GetLocallDDT();
  }
  PrintlDDTPerResidueStats(_local_scores,
                           false,
                           settings.cutoffs.size());
}

std::vector<EntityView> lDDTScorer::GetReferences(){
  return references_view;
}

void lDDTScorer::_PrepareGlobalRDMap(){
  glob_dist_list = PreparelDDTGlobalRDMap(references_view,
                                          settings.cutoffs,
                                          settings.sequence_separation,
                                          settings.radius);
}

bool lDDTScorer::IsValid(){
  return _score_valid;
}

void lDDTScorer::_ComputelDDT(){
  std::pair<int,int> cov = ComputeCoverage(model_view,glob_dist_list);
  if (cov.second == -1) {
    LOG_INFO("Coverage: 0 (0 out of 0 residues)");
  } else {
    std::stringstream sout;
    sout << "Coverage: " << (float(cov.first)/float(cov.second)) << " (" << cov.first << " out of " << cov.second << " residues)";
    LOG_INFO(sout.str());
  }

  if (cov.first == 0) {
    _num_tot_con = 0;
    _num_cons_con = 0;
    _global_score = 0.0;
    _score_calculated = true;
    _score_valid = false;
  }

  std::pair<int,int> total_ov=alg::LocalDistDiffTest(model_view, glob_dist_list, settings.cutoffs, settings.sequence_separation, settings.label);
  Real lddt = static_cast<Real>(total_ov.first)/(static_cast<Real>(total_ov.second) ? static_cast<Real>(total_ov.second) : 1);
  _num_tot_con = total_ov.second ? total_ov.second : 1;
  _num_cons_con = total_ov.first;
  _global_score = lddt;
  _score_calculated = true;
  _score_valid = true;
}

void lDDTScorer::_GetLocallDDT(){
  if (!_score_calculated){
    _ComputelDDT();
  }
  _local_scores = GetlDDTPerResidueStats(model_view,
                                         glob_dist_list,
                                         false, // do not print structural checks
                                         settings.label);
  _has_local_scores = true;
}


GlobalRDMap CreateDistanceList(const EntityView& ref,Real max_dist)
{
 GlobalRDMap dist_list; 
 if (!ref.GetChainCount()) {
   return dist_list;
 }
 ResidueViewList ref_residues=ref.GetChainList()[0].GetResidueList(); 
 for (ResidueViewList::iterator i=ref_residues.begin(), e=ref_residues.end(); i!=e; ++i) {
   if (IsStandardResidue(i->GetName())) {
     ResidueRDMap res_dist_list;
     ResNum rnum = i->GetNumber();  
     AtomViewList ref_atoms=i->GetAtomList();
     AtomViewList within;
     if (max_dist<0){
       within=ref.GetAtomList();
     }  
     for (AtomViewList::iterator ai=ref_atoms.begin(), ae=ref_atoms.end(); ai!=ae; ++ai) {
       UniqueAtomIdentifier first_atom(ai->GetResidue().GetChain().GetName(),ai->GetResidue().GetNumber(),ai->GetResidue().GetName(),ai->GetName());
       if (ai->GetElement()=="H") { continue; }
       if (max_dist>=0){ 
         within=ref.FindWithin(ai->GetPos(), max_dist);
       }
       for (AtomViewList::iterator aj=within.begin(), ae2=within.end(); aj!=ae2; ++aj) {      
         UniqueAtomIdentifier second_atom(aj->GetResidue().GetChain().GetName(),aj->GetResidue().GetNumber(),aj->GetResidue().GetName(),aj->GetName());
         if (aj->GetElement()=="H" ||
             aj->GetResidue().GetChain()!=ai->GetResidue().GetChain()) {
             continue;
         }
         Real dist=geom::Length(ai->GetPos()-aj->GetPos());
         UAtomIdentifiers atoms = std::make_pair(first_atom,second_atom); 
         std::pair<Real,Real> values = std::make_pair(dist,dist);  
         res_dist_list[atoms]=values;
       }
     }   
     dist_list[rnum]=res_dist_list;
   }
 } 
 return dist_list;
} 

GlobalRDMap CreateDistanceListFromMultipleReferences(const std::vector<EntityView>& ref_list, std::vector<Real>& cutoff_list, int sequence_separation, Real max_dist)
{
  int ref_counter=0;  
  ExistenceMap ex_map;  
  GlobalRDMap glob_dist_list = CreateDistanceList(ref_list[0],max_dist);
  update_existence_map (ex_map,ref_list[0],ref_counter);
  ref_counter++;  
  for (std::vector<EntityView>::const_iterator ref_list_it=ref_list.begin()+1;ref_list_it!=ref_list.end();++ref_list_it) {
       EntityView ref = *ref_list_it;
       std::vector<std::pair<long int, long int> > overlap_list(ref.GetResidueCount(), std::pair<long int, long int>(0, 0));
       check_and_swap(glob_dist_list,ref,cutoff_list,sequence_separation,overlap_list);
       GlobalRDMap new_dist_list=CreateDistanceList(ref,max_dist);
       merge_distance_lists(glob_dist_list,new_dist_list,ex_map,ref,ref_counter);
       update_existence_map (ex_map,ref,ref_counter); 
       ref_counter++;
  }    
  return glob_dist_list;  
}

void PrintResidueRDMap(const ResidueRDMap& res_dist_list)
{
  for (ResidueRDMap::const_iterator res_dist_list_it = res_dist_list.begin();res_dist_list_it!=res_dist_list.end();++res_dist_list_it) {
    UAtomIdentifiers uais = res_dist_list_it->first;
    std::pair<Real,Real> minmax = res_dist_list_it->second;
    std::cout << uais.first.GetChainName() << " " << uais.first.GetResNum() << " " << uais.first.GetResidueName() << " " << uais.first.GetAtomName() << " " <<
                 uais.second.GetChainName() << " " << uais.second.GetResNum() << " " << uais.second.GetResidueName() << " " << uais.second.GetAtomName() << " " <<
                 minmax.first << " " << minmax.second << std::endl;       
  }
}


void PrintGlobalRDMap(const GlobalRDMap& glob_dist_list){
  for (GlobalRDMap::const_iterator glob_dist_list_it = glob_dist_list.begin();glob_dist_list_it!=glob_dist_list.end();++glob_dist_list_it) {
    if (glob_dist_list_it->second.size()!=0) {  
      PrintResidueRDMap(glob_dist_list_it->second);
    }    
  }    
}



std::pair<long int,long int> LocalDistDiffTest(const EntityView& mdl, const GlobalRDMap& glob_dist_list,
                   std::vector<Real> cutoff_list, int sequence_separation, const String& local_lddt_property_string)
{
  if (!mdl.GetResidueCount()) {
    LOG_WARNING("model structures doesn't contain any residues");
    return std::make_pair(0,0);
  }
  if (glob_dist_list.size()==0) {
    LOG_WARNING("global reference list is empty");
    return std::make_pair(0,0);
  }
  std::vector<std::pair<long int, long int> > overlap_list(mdl.GetResidueCount(), std::pair<long int, long int>(0, 0));
  check_and_swap(glob_dist_list,mdl,cutoff_list,sequence_separation,overlap_list);
  ChainView mdl_chain=mdl.GetChainList()[0];  
  overlap_list.clear();
  std::pair<long int, long int> total_ov(0, 0);
  for (GlobalRDMap::const_iterator i=glob_dist_list.begin(), e=glob_dist_list.end(); i!=e; ++i) {
    ResNum rn = i->first;
    if (i->second.size()!=0) {
      std::pair<long int, long int> ov1=calc_overlap1(i->second, rn, mdl_chain, sequence_separation, cutoff_list, 
                                            false, false, overlap_list,true);
      total_ov.first+=ov1.first;
      total_ov.second+=ov1.second;
    }

  }
  for (GlobalRDMap::const_iterator i=glob_dist_list.begin(),
       e=glob_dist_list.end();i!=e; ++i) {
    ResNum rn = i->first;
    if(local_lddt_property_string!="") {
      ResidueView mdlr=mdl_chain.FindResidue(rn);
      if (mdlr.IsValid()) {
        int mdl_res_index =mdlr.GetIndex();
        Real local_lddt=static_cast<Real>(overlap_list[mdl_res_index].first)/(static_cast<Real>(overlap_list[mdl_res_index].second) ? static_cast<Real>(overlap_list[mdl_res_index].second) : 1);
        mdlr.SetFloatProp(local_lddt_property_string, local_lddt);
        mdlr.SetIntProp(local_lddt_property_string+"_conserved", overlap_list[mdl_res_index].first);
        mdlr.SetIntProp(local_lddt_property_string+"_total", overlap_list[mdl_res_index].second);
      }
    }
  }
  overlap_list.clear();
  return std::make_pair(total_ov.first,total_ov.second);
}

Real LocalDistDiffTest(const EntityView& mdl, const EntityView& target, Real cutoff, Real max_dist, const String& local_lddt_property_string)
{
   std::vector<Real> cutoffs;
   cutoffs.push_back(cutoff);
   GlobalRDMap glob_dist_list = CreateDistanceList(target,max_dist);
   int sequence_separation = 0; 
   std::pair<long int,long int>  total_ov = LocalDistDiffTest(mdl, glob_dist_list, cutoffs, sequence_separation, local_lddt_property_string);
   return static_cast<Real>(total_ov.first)/(static_cast<Real>(total_ov.second) ? static_cast<Real>(total_ov.second) : 1);
}

Real LocalDistDiffTest(const EntityView& v,
                       std::vector<EntityView>& ref_list,
                       const GlobalRDMap& glob_dist_list,
                       lDDTSettings& settings) {

  std::pair<int,int> cov = ComputeCoverage(v,glob_dist_list);
  if (cov.second == -1) {
    std::cout << "Coverage: 0 (0 out of 0 residues)" << std::endl;
  } else {
    std::cout << "Coverage: " << (float(cov.first)/float(cov.second)) << " (" << cov.first << " out of " << cov.second << " residues)" << std::endl;
  }

  if (cov.first==0) {
    std::cout << "Global LDDT score: 0.0" << std::endl;
    return 0.0;
  }

  std::pair<int,int> total_ov=alg::LocalDistDiffTest(v, glob_dist_list, settings.cutoffs, settings.sequence_separation, settings.label);
  Real lddt = static_cast<Real>(total_ov.first)/(static_cast<Real>(total_ov.second) ? static_cast<Real>(total_ov.second) : 1);
  std::cout << "Global LDDT score: " << std::setprecision(4) << lddt << std::endl;
  std::cout << "(" << std::fixed << total_ov.first << " conserved distances out of " << total_ov.second
            << " checked, over " << settings.cutoffs.size() << " thresholds)" << std::endl;
  return lddt;
}

Real LocalDistDiffTest(const ost::seq::AlignmentHandle& aln,
                   Real cutoff, Real max_dist, int ref_index, int mdl_index)
{
  seq::ConstSequenceHandle ref_seq=aln.GetSequence(ref_index);
  seq::ConstSequenceHandle mdl_seq=aln.GetSequence(mdl_index);  
  if (!ref_seq.HasAttachedView()) {
    LOG_ERROR("reference sequence doesn't have a view attached.");
    return 0.0;
  }
  if (!mdl_seq.HasAttachedView()) {
    LOG_ERROR("model sequence doesn't have a view attached");
    return 0.0;
  }
  XCSEditor edi=ref_seq.GetAttachedView().GetHandle().EditXCS(BUFFERED_EDIT);
  for (int i=0; i<aln.GetLength(); ++i) {
    ResidueView mdl_res=mdl_seq.GetResidue(i);
    if (!mdl_res) {
      continue;
    }
    String rname=mdl_res.GetName();
    if (!(rname=="GLU" || rname=="ASP" || rname=="VAL" || rname=="TYR" || 
         rname=="PHE" || rname=="LYS" || rname=="ARG")) {
      continue;
    }
    std::pair<Real, Real> ov1=calc_overlap2(ref_seq, mdl_seq, i,
                                            cutoff, max_dist, true, 
                                            false);
    std::pair<Real, Real> ov2=calc_overlap2(ref_seq, mdl_seq, i,
                                            cutoff, max_dist, true, 
                                            true);
    if (ov1.first/ov1.second<ov2.first/ov2.second) {
     AtomViewList atoms=mdl_res.GetAtomList();
     for (AtomViewList::iterator j=atoms.begin(), 
          e2=atoms.end(); j!=e2; ++j) {
       if (Swappable(rname, j->GetName())) {
         edi.RenameAtom(j->GetHandle(), SwappedName(j->GetName()));
       }
     }
    }
  }
  std::pair<Real, Real> total_ov(0.0, 0.0);
  for (int i=0; i<aln.GetLength(); ++i) {
     std::pair<Real, Real> ov1=calc_overlap2(ref_seq, mdl_seq, i, cutoff, 
                                             max_dist, false, false);
     total_ov.first+=ov1.first;
     total_ov.second+=ov1.second;
  }
   return total_ov.first/(total_ov.second ? total_ov.second : 1);
}

Real LDDTHA(EntityView& v, const GlobalRDMap& global_dist_list, int sequence_separation)
{
    std::vector<Real> cutoffs;
    cutoffs.push_back(0.5);
    cutoffs.push_back(1.0);
    cutoffs.push_back(2.0);
    cutoffs.push_back(4.0);
    String label="locallddt";
    std::pair<long int,long int> total_ov=alg::LocalDistDiffTest(v, global_dist_list, cutoffs,  sequence_separation, label);
    return static_cast<Real>(total_ov.first)/(static_cast<Real>(total_ov.second) ? static_cast<Real>(total_ov.second) : 1);
}


void CleanlDDTReferences(std::vector<EntityView>& ref_list){
  for (unsigned int i=0;i<ref_list.size();i++) {
      if (ref_list[0].GetChainList()[0].GetName()!=ref_list[i].GetChainList()[0].GetName()) {
        std::cout << "ERROR: First chains in the reference structures have different names" << std::endl;
        exit(-1);
    }
    ref_list[i] = ref_list[i].GetChainList()[0].Select("peptide=true");
  }
}

GlobalRDMap PreparelDDTGlobalRDMap(const std::vector<EntityView>& ref_list,
                                   std::vector<Real>& cutoff_list,
                                   int sequence_separation,
                                   Real max_dist){
  GlobalRDMap glob_dist_list;
  if (ref_list.size()==1) {
    glob_dist_list = CreateDistanceList(ref_list[0], max_dist);
  } else {
    glob_dist_list = CreateDistanceListFromMultipleReferences(ref_list,
                                                              cutoff_list,
                                                              sequence_separation,
                                                              max_dist);
  }

  return glob_dist_list;
}

void CheckStructure(EntityView& ent,
                    StereoChemicalParams& bond_table,
                    StereoChemicalParams& angle_table,
                    ClashingDistances& nonbonded_table,
                    Real bond_tolerance,
                    Real angle_tolerance){
  // performs structural checks and filters the structure
  StereoChemistryInfo stereo_chemistry_info;
  try { 
    std::pair<EntityView,StereoChemistryInfo> csc_result = alg::CheckStereoChemistry(ent,bond_table,angle_table,bond_tolerance,angle_tolerance);
    ent = csc_result.first;
    stereo_chemistry_info = csc_result.second;
  } catch (std::exception& e) {       
    std::cout << "An error occurred during the structure quality checks, stage 1:" << std::endl;    
    std::cout << e.what() << std::endl;
    exit(-1);
  }
  std::cout << "Average Z-Score for bond lengths: " << std::fixed << std::setprecision(5) << stereo_chemistry_info.GetAvgZscoreBonds() << std::endl;
  std::cout << "Bonds outside of tolerance range: " << stereo_chemistry_info.GetBadBondCount() << " out of " << stereo_chemistry_info.GetBondCount() << std::endl;
  std::cout << "Bond\tAvg Length\tAvg zscore\tNum Bonds" << std::endl;
  std::map<String,BondLengthInfo> avg_bond_length_info = stereo_chemistry_info.GetAvgBondLengthInfo();
  for (std::map<String,BondLengthInfo>::const_iterator abli_it=avg_bond_length_info.begin();abli_it!=avg_bond_length_info.end();++abli_it) {
    String key = (*abli_it).first;
    BondLengthInfo bond_length_info = (*abli_it).second;
    std::cout << key << "\t" << std::fixed << std::setprecision(5) << std::left << std::setw(10) <<
                 bond_length_info.GetAvgLength() << "\t" << std::left << std::setw(10) << bond_length_info.GetAvgZscore() << "\t" << bond_length_info.GetCount()  << std::endl;
  }
  std::cout << "Average Z-Score angle widths: " << std::fixed << std::setprecision(5) << stereo_chemistry_info.GetAvgZscoreAngles() << std::endl;
  std::cout << "Angles outside of tolerance range: " << stereo_chemistry_info.GetBadAngleCount() << " out of " << stereo_chemistry_info.GetAngleCount() << std::endl;
  ClashingInfo clash_info;
  try {
    std::pair<EntityView,ClashingInfo> fc_result = alg::FilterClashes(ent,nonbonded_table);
    ent = fc_result.first;
    clash_info = fc_result.second;
  } catch (std::exception& e) {
    std::stringstream serr;
    serr << "An error occurred during the structure quality checks, stage 2: " << e.what();    
    throw ost::Error(serr.str());
  }
  std::cout << clash_info.GetClashCount() << " non-bonded short-range distances shorter than tolerance distance" << std::endl;
  std::cout << "Distances shorter than tolerance are on average shorter by: " << std::fixed << std::setprecision(5) << clash_info.GetAverageOffset() << std::endl;
}

std::vector<lDDTLocalScore> GetlDDTPerResidueStats(EntityView& model,
                                                   GlobalRDMap& glob_dist_list,
                                                   bool structural_checks,
                                                   String label){
  std::vector<lDDTLocalScore> scores;
  EntityView outv = model;
  for (ChainViewList::const_iterator ci = outv.GetChainList().begin(),
       ce = outv.GetChainList().end(); ci != ce; ++ci) {
    for (ResidueViewList::const_iterator rit = ci->GetResidueList().begin(),
         re = ci->GetResidueList().end(); rit != re; ++rit) {

      ResidueView ritv=*rit;
      ResNum rnum = ritv.GetNumber();
      bool assessed = false;
      String assessed_string="No";
      String quality_problems_string;
      if (structural_checks) {
        quality_problems_string="No";
      } else {
        quality_problems_string="NA";
      }
      Real lddt_local = -1;
      String lddt_local_string="-";
      int conserved_dist = -1;
      int total_dist = -1;
      if (IsResnumInGlobalRDMap(rnum,glob_dist_list)) {
        assessed = true;
        assessed_string="Yes";
      }
      if (ritv.HasProp("stereo_chemical_violation_sidechain") || 
          ritv.HasProp("steric_clash_sidechain")) {
        quality_problems_string="Yes";
      }
      if (ritv.HasProp("stereo_chemical_violation_backbone") || 
          ritv.HasProp("steric_clash_backbone")) {
        quality_problems_string="Yes+";
      }

      if (assessed==true) {
        if (ritv.HasProp(label)) {
          lddt_local=ritv.GetFloatProp(label);
          std::stringstream stkeylddt;
          stkeylddt <<  std::fixed << std::setprecision(4) << lddt_local;
          lddt_local_string=stkeylddt.str();
          conserved_dist=ritv.GetIntProp(label+"_conserved");
          total_dist=ritv.GetIntProp(label+"_total");
        } else {
          //std::cout << label << std::endl;
          lddt_local = 0;
          lddt_local_string="0.0000";
          conserved_dist = 0;
          total_dist = 0;
        }
      }
      // std::tuple<String, String, int, String, String, Real, int, int>
      lDDTLocalScore row(ritv.GetChain().GetName(),
                         ritv.GetName(),
                         ritv.GetNumber().GetNum(),
                         assessed_string,
                         quality_problems_string,
                         lddt_local,
                         conserved_dist,
                         total_dist);
      scores.push_back(row);
    }
  }

  return scores;
}

void PrintlDDTPerResidueStats(std::vector<lDDTLocalScore>& scores, bool structural_checks, int cutoffs_length){
  if (structural_checks) {
    std::cout << "Local LDDT Scores:" << std::endl;
    std::cout << "(A 'Yes' in the 'Quality Problems' column stands for problems" << std::endl;
    std::cout << "in the side-chain of a residue, while a 'Yes+' for problems" << std::endl;
    std::cout << "in the backbone)" << std::endl;
  } else {
    std::cout << "Local LDDT Scores:" << std::endl;
  }
  std::cout << lDDTLocalScore::GetHeader(structural_checks, cutoffs_length) << std::endl;
  for (std::vector<lDDTLocalScore>::const_iterator sit = scores.begin(); sit != scores.end(); ++sit) {
    std::cout << sit->ToString(structural_checks) << std::endl;
  }
  std::cout << std::endl;
}

// debugging code
/*
Real OldStyleLDDTHA(EntityView& v, const GlobalRDMap& global_dist_list)
{
    Real lddt =0;
    std::vector<Real> cutoffs05;
    cutoffs05.push_back(0.5);
    std::pair<Real,Real> lddt05o=alg::LocalDistDiffTest(v, global_dist_list, cutoffs05, "locallddt0.5");
    std::vector<Real> cutoffs1;
    cutoffs1.push_back(1.0);
    std::pair<Real,Real> lddt1o=alg::LocalDistDiffTest(v, global_dist_list, cutoffs1, "locallddt1");
    std::vector<Real> cutoffs2;
    cutoffs2.push_back(2.0);
    std::pair<Real,Real> lddt2o=alg::LocalDistDiffTest(v, global_dist_list, cutoffs2, "locallddt2");
    std::vector<Real> cutoffs4;
    cutoffs4.push_back(4.0);
    std::pair<Real,Real> lddt4o=alg::LocalDistDiffTest(v, global_dist_list, cutoffs4, "locallddt4");
    Real lddt05 = lddt05o.first/(lddt05o.second ? lddt05o.second : 1);
    Real lddt1 = lddt1o.first/(lddt1o.second ? lddt1o.second : 1);
    Real lddt2 = lddt2o.first/(lddt2o.second ? lddt2o.second : 1);
    Real lddt4 = lddt4o.first/(lddt4o.second ? lddt4o.second : 1);            
    lddt = (lddt05+lddt1+lddt2+lddt4)/4.0;
    for (ResidueViewIter rit=v.ResiduesBegin();rit!=v.ResiduesEnd();++rit){
      ResidueView ritv = *rit;
      if (ritv.HasProp("locallddt0.5")) {
        Real lddt_local = (ritv.GetFloatProp("locallddt0.5")+ritv.GetFloatProp("locallddt1")+ritv.GetFloatProp("locallddt2")+ritv.GetFloatProp("locallddt4"))/4.0;
        ritv.SetFloatProp("locallddt",lddt_local);
      }    
    }    
    return lddt;
}
*/

}}}
