#include <ost/log.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/alg/distance_rmsd_test.hh>
#include <ost/mol/alg/distance_test_common.hh>
#include <boost/concept_check.hpp>

namespace ost { namespace mol { namespace alg {

namespace {


Real calc_squared_diff(const std::vector<geom::Vec3>& pos_one,
                       const std::vector<geom::Vec3>& pos_two,
                       const std::vector<Real>& exp_dist,
                       Real cap){
  //expect all all vectors to have the same length...
  Real squared_diff = 0.0;
  Real diff;
  for(uint i = 0; i < pos_one.size(); ++i){
    diff = std::min(std::abs(geom::Distance(pos_one[i],pos_two[i]) - exp_dist[i]),cap);
    squared_diff += diff*diff;
  }
  return squared_diff;
}

void fill_residue_positions(const ResidueRDMap& res_dist_list, const ResidueView& res, 
                            const ChainView& mdl_chain, Real cap_distance, int sequence_separation, 
                            std::map<String,geom::Vec3>& positions){

  AtomViewList atom_list = res.GetAtomList();
  //let's already fill in the positions and check later on whether a swap is necessary
  for(AtomViewList::const_iterator i = atom_list.begin(); 
      i != atom_list.end(); ++i){
    positions[i->GetName()] = i->GetPos();
  }

  String rname = res.GetName();
  if(!(rname=="GLU" || rname=="ASP" || rname=="VAL" || rname=="TYR" ||
      rname=="PHE" || rname=="LEU" || rname=="ARG")){
    return; //no ambiguity check necessary for all other amino acids 
  }

  std::map<String,geom::Vec3> original_positions;
  std::map<String,geom::Vec3> swapped_positions;
  AtomView atom;
  for(AtomViewList::const_iterator i = atom_list.begin(); i !=atom_list.end(); ++i){
    if(Swappable(rname,i->GetName())){
      atom = res.FindAtom(SwappedName(i->GetName()));
      if(atom){ //both atoms are valid => note, that only in this case a swap is possible
        original_positions[i->GetName()] = i->GetPos();
        swapped_positions[i->GetName()] = atom.GetPos();
      }
    }
  }

  std::vector<geom::Vec3> or_positions;
  std::vector<geom::Vec3> sw_positions;
  std::vector<geom::Vec3> other_positions;
  std::vector<Real> expected_dist;
  ResidueView other_residue;

  for(ResidueRDMap::const_iterator i = res_dist_list.begin(); 
      i != res_dist_list.end(); ++i){
    //check whether first atom is one of the swappables
    const UniqueAtomIdentifier& first_atom = i->first.first;
    if(original_positions.find(first_atom.GetAtomName()) == original_positions.end()) continue;
    //we only check the second atom of the distance if it is NOT swappable
    const UniqueAtomIdentifier& second_atom = i->first.second;
    if(Swappable(second_atom.GetResidueName(),second_atom.GetAtomName())) continue;
    //check whether other atom is valid and check the sequence separation
    other_residue = mdl_chain.FindResidue(second_atom.GetResNum());
    if(!other_residue) continue;
    if(std::abs(res.GetNumber().GetNum() - other_residue.GetNumber().GetNum()) <= sequence_separation) continue;
    atom = other_residue.FindAtom(second_atom.GetAtomName());
    if(!atom) continue;   

    or_positions.push_back(original_positions[first_atom.GetAtomName()]);
    sw_positions.push_back(swapped_positions[first_atom.GetAtomName()]);
    expected_dist.push_back(i->second.first);
    other_positions.push_back(atom.GetPos());
  }

  Real diff_orig = calc_squared_diff(or_positions,other_positions,expected_dist, cap_distance);
  Real diff_swapped = calc_squared_diff(sw_positions,other_positions,expected_dist, cap_distance);
  //finally swap the stuff if necessary
  if(diff_swapped < diff_orig){
    for(std::map<String,geom::Vec3>::iterator i = swapped_positions.begin();
        i != swapped_positions.end(); ++i){
      positions[i->first] = i->second;
    }
  }
}

void fill_positions(const GlobalRDMap& glob_dist_list, const ChainView& mdl_chain,
                    Real cap_distance, int sequence_separation,
                    std::map<ResNum, std::map<String, geom::Vec3> >& positions){

  ResidueView actual_res;
  for(GlobalRDMap::const_iterator i = glob_dist_list.begin(); 
      i != glob_dist_list.end(); ++i){
    actual_res = mdl_chain.FindResidue(i->first);
    if(actual_res){
      std::map<String,geom::Vec3> residue_positions;
      fill_residue_positions(i->second, actual_res, mdl_chain,
                             cap_distance, sequence_separation,
                             residue_positions);
      if(!residue_positions.empty()) positions[i->first] = residue_positions;
    }
  }
}

void fill_values(const GlobalRDMap& glob_dist_list, 
                 const std::map<ResNum, std::map<String, geom::Vec3> >& positions,
                 Real cap_distance, int sequence_separation,
                 std::map<ResNum, std::pair<Real,long int> >& drmsd_values){

  ResNum rnum_one, rnum_two;
  std::vector<geom::Vec3> positions_one, positions_two;
  std::vector<Real> expected_distances;
  std::map<ResNum, std::map<String, geom::Vec3> >::const_iterator res_positions_one, res_positions_two;
  std::map<String,geom::Vec3>::const_iterator pos_one, pos_two;
  int penalties;
  Real squared_diff;
  int num_interactions;
  Real squared_cap_distance = cap_distance * cap_distance;

  for(GlobalRDMap::const_iterator i = glob_dist_list.begin(); 
      i != glob_dist_list.end(); ++i){
    rnum_one = i->first;
    res_positions_one = positions.find(rnum_one);
    if(res_positions_one == positions.end()){
      //this residue is missing, so we give it the full penalty
      drmsd_values[rnum_one] = std::make_pair(i->second.size()*squared_cap_distance,i->second.size());
      continue;
    }

    positions_one.clear();
    positions_two.clear();
    expected_distances.clear();
    penalties = 0;
    
    for(ResidueRDMap::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
      //check sequence separation
      rnum_two = j->first.second.GetResNum();
      if(std::abs(rnum_one.GetNum() - rnum_two.GetNum()) <= sequence_separation) continue;
      //if residue is not present in the positions map, the positions are also missing in the
      //model => penalty
      res_positions_two = positions.find(rnum_two);
      if(res_positions_two == positions.end()){
        ++penalties;
        continue;
      }
      //check whether both positions are present in the positions map
      pos_one = res_positions_one->second.find(j->first.first.GetAtomName());
      pos_two = res_positions_two->second.find(j->first.second.GetAtomName());
      if(pos_one == res_positions_one->second.end() || pos_two == res_positions_two->second.end()){
        ++penalties; //one of the atoms is missing
        continue;
      }
      positions_one.push_back(pos_one->second);
      positions_two.push_back(pos_two->second);
      expected_distances.push_back(j->second.first);
    }
    squared_diff = calc_squared_diff(positions_one, positions_two,
                                     expected_distances, cap_distance);
    num_interactions = positions_one.size();
    //add the penalty values
    if(penalties > 0){
      squared_diff += squared_cap_distance*penalties;
      num_interactions += penalties;
    }
    drmsd_values[rnum_one] = std::make_pair(squared_diff,num_interactions);
  }
}

}//ns

std::pair<Real,long int> DistanceRMSDTest(const EntityView& mdl, const GlobalRDMap& glob_dist_list,
                                          Real cap_distance, int sequence_separation, 
                                          const String& local_drmsdt_property_string){
  if (!mdl.GetResidueCount()) {
    LOG_WARNING("model structures doesn't contain any residues");
    return std::make_pair<long int,long int>(0,0);
  }
  if (glob_dist_list.size()==0) {
    LOG_WARNING("global reference list is empty");
    return std::make_pair<long int,long int>(0,0);
  }
  
  //every residue gets a map of positions
  std::map<ResNum, std::map<String, geom::Vec3> > positions;
  //every residue will have a sum of squared distances and number of distances
  std::map<ResNum, std::pair<Real,long int> > drmsd_values;
  //currently we're only looking at the first chain
  ChainView mdl_chain=mdl.GetChainList()[0];

  fill_positions(glob_dist_list, mdl_chain, cap_distance, 
                 sequence_separation,positions);
  fill_values(glob_dist_list, positions, cap_distance, 
              sequence_separation, drmsd_values);

  Real total_squared_diff = 0.0;
  long int total_num_interactions = 0;
  Real drmsd;
  std::pair<Real, long int> residue_data;
  ResidueView actual_residue;
  for(GlobalRDMap::const_iterator i = glob_dist_list.begin();
      i != glob_dist_list.end(); ++i){
    residue_data = drmsd_values[i->first];
    total_squared_diff += residue_data.first;
    total_num_interactions += residue_data.second;
    actual_residue = mdl_chain.FindResidue(i->first);
    if(actual_residue){
      drmsd = residue_data.second > 0 ? std::sqrt(residue_data.first/residue_data.second) : 0;
      actual_residue.SetFloatProp(local_drmsdt_property_string, drmsd);
      actual_residue.SetFloatProp(local_drmsdt_property_string+"_sum",residue_data.first);
      actual_residue.SetFloatProp(local_drmsdt_property_string+"_count",residue_data.second);
    }
  }
  drmsd = total_num_interactions > 0 ? std::sqrt(total_squared_diff/total_num_interactions) : 0;
  return std::make_pair(drmsd,total_num_interactions);
}


Real DRMSD(const EntityView& v, const GlobalRDMap& global_dist_list, Real cap_distance, int sequence_separation)
{
    String label="localdrmsd";
    std::pair<Real,long int> result=alg::DistanceRMSDTest(v, global_dist_list, cap_distance, sequence_separation, label);
    return result.first;
}

}}} //ns
