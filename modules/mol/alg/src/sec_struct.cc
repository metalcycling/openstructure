#include <ost/mol/alg/sec_struct.hh>

#include <ost/mol/sec_structure.hh>

#include <limits>
#include <iostream>

namespace {

struct ExtendedBridge {

  ExtendedBridge() {
    // should be enough for most strands
    strand_one.reserve(16);
    strand_two.reserve(16);
  }

  bool parallel; // antiparallel otherwise
  std::vector<int> strand_one;
  std::vector<int> strand_two;
};

bool Connected(const std::vector<int>& connected_to_next,
               int start_idx, int length) {
  if(start_idx < 0) return false;
  if(start_idx + length > static_cast<int>(connected_to_next.size())) {
    return false;
  }
  for(int i = 0; i < (length - 1); ++i) {
    if(connected_to_next[start_idx + i] == 0) return false;
  }
  return true;
}

bool AddIfValid(int idx_one, int idx_two, int start_idx, int end_idx,
                const std::vector<int>& connected_to_next,
                bool parallel,
                std::vector<ExtendedBridge>& extended_bridges) {

  // this function assumes idx_one to be sorted in an increasing manner 
  // when called multiple times!

  // idx_one is expected to be in range [start_idx, end_idx]
  // if idx_two is in the same range we can skip cases where
  // idx_one > idx_two, since the interactions are symmetric and we would
  // see every bridge twice. If idx_two is outside this range 
  // (e.g. when we only calculate the secondary structure for a stretch 
  // of a larger structure) we have to add it anyway.

  if(idx_two >= start_idx && idx_two <= end_idx && idx_two < idx_one) {
    return false;
  }

  // the abs in the following check is still needed for the case when
  // idx_two is outside [start_idx, end_idx]
  if(std::abs(idx_two - idx_one) >= 3 &&
     Connected(connected_to_next, idx_one - 1, 3) && 
     Connected(connected_to_next, idx_two - 1, 3)) {

    // it is a bridge! Let's see whether we can extend
    // an existing bridge...
    bool added = false;

    for(uint i = 0; i < extended_bridges.size(); ++i) {

      if(extended_bridges[i].parallel == parallel && 
         idx_one == extended_bridges[i].strand_one.back() + 1) {

        if(parallel && idx_two == extended_bridges[i].strand_two.back() + 1) {
          extended_bridges[i].strand_one.push_back(idx_one);
          extended_bridges[i].strand_two.push_back(idx_two);
          added = true;
          break;
        }

        if(!parallel && idx_two == extended_bridges[i].strand_two.back() - 1) {
          extended_bridges[i].strand_one.push_back(idx_one);
          extended_bridges[i].strand_two.push_back(idx_two);
          added = true;
          break;         
        }
      }
    }

    if(!added) {
      // haven't seen it... Lets add a new one!
      extended_bridges.push_back(ExtendedBridge());
      extended_bridges.back().parallel = parallel;
      extended_bridges.back().strand_one.push_back(idx_one);
      extended_bridges.back().strand_two.push_back(idx_two);
    }
    return true;
  }
  return false;
}

void ResolveBulges(const std::vector<int>& connected_to_next,
                   std::vector<ExtendedBridge>& extended_bridges) {
  
  for(uint i = 0; i < extended_bridges.size(); ++i) {
    for(uint j = i + 1; j < extended_bridges.size(); ++j) {

      if(extended_bridges[i].parallel == extended_bridges[j].parallel) {
    
        int i_strand_one_begin = extended_bridges[i].strand_one[0];
        int i_strand_one_end   = extended_bridges[i].strand_one.back();
        int i_strand_two_begin = extended_bridges[i].strand_two[0];
        int i_strand_two_end   = extended_bridges[i].strand_two.back();
        int j_strand_one_begin = extended_bridges[j].strand_one[0];
        int j_strand_one_end   = extended_bridges[j].strand_one.back();
        int j_strand_two_begin = extended_bridges[j].strand_two[0];
        int j_strand_two_end   = extended_bridges[j].strand_two.back();

        // first strand must not overlap
        if(i_strand_one_end >= j_strand_one_begin &&
           i_strand_one_begin <= j_strand_one_end) {
             continue;
        }

        // both strands must be connected, the first strand is the
        // same for parallel and anti-parallel bridges
        int begin = std::min(i_strand_one_begin, j_strand_one_begin);
        int end = std::max(i_strand_one_end, j_strand_one_end);
        if(!Connected(connected_to_next, begin, end - begin)) continue;

        // the distances are dependent on parallel / anti-parallel
        int d_one = -1;
        int d_two = -1;

        if(extended_bridges[i].parallel) {
          // check connectivity for second strand
          begin = std::min(i_strand_two_begin, j_strand_two_begin);
          end = std::max(i_strand_two_end, j_strand_two_end);
          if(!Connected(connected_to_next, begin, end - begin)) continue;

          d_one = j_strand_one_begin - i_strand_one_end;
          d_two = j_strand_two_begin - i_strand_two_end;
        }
        else {
          // check connectivity for second strand
          begin = std::min(i_strand_two_end, j_strand_two_end);
          end = std::max(i_strand_two_begin, j_strand_two_begin);
          if(!Connected(connected_to_next, begin, end - begin)) continue;
          
          d_one = j_strand_one_begin - i_strand_one_end;
          d_two = i_strand_two_end - j_strand_two_begin;
        }

        if(d_one < 0 || d_two < 0) continue;

        if((d_one < 6 && d_two < 3) || (d_one < 3 && d_two < 6)) {

          extended_bridges[i].strand_one.insert(
                                  extended_bridges[i].strand_one.end(),
                                  extended_bridges[j].strand_one.begin(),
                                  extended_bridges[j].strand_one.end());
          extended_bridges[i].strand_two.insert(
                                  extended_bridges[i].strand_two.end(),
                                  extended_bridges[j].strand_two.begin(),
                                  extended_bridges[j].strand_two.end());
          extended_bridges.erase(extended_bridges.begin() + j);
          --j;
        }
      }
    }
  }
}

void AssignBridges(int start_idx, int size,
	                 const std::vector<int>& donor_for_one,
                   const std::vector<int>& donor_for_two, 
                   const std::vector<int>& connected_to_next,
                   String& ss_assignment) {

  int full_size = donor_for_one.size();
  int end_idx = start_idx + size - 1;

  std::vector<ExtendedBridge> extended_bridges;

  // let's search for bridges

  std::vector<int> parallel_partners;
  std::vector<int> antiparallel_partners;

  for(int i = std::max(0, start_idx - 5); 
      i < std::min(start_idx + size + 5, full_size); ++i) {

    for(int j = 0; j < full_size; ++j) {

      // check for parallel bridge
      if(i > 0 && i < full_size - 1) {
        if((donor_for_one[i+1] == j && donor_for_one[j] == i - 1) ||
           (donor_for_one[i+1] == j && donor_for_two[j] == i - 1) ||
           (donor_for_two[i+1] == j && donor_for_one[j] == i - 1) ||
           (donor_for_two[i+1] == j && donor_for_two[j] == i - 1)) {
          if(AddIfValid(i,j,start_idx, end_idx, connected_to_next,
                        true, extended_bridges)) continue;
        } 
      }
      if(j > 0 && j < full_size - 1) {
        if((donor_for_one[j+1] == i && donor_for_one[i] == j - 1) ||
           (donor_for_one[j+1] == i && donor_for_two[i] == j - 1) ||
           (donor_for_two[j+1] == i && donor_for_one[i] == j - 1) ||
           (donor_for_two[j+1] == i && donor_for_two[i] == j - 1)) {
          if(AddIfValid(i,j,start_idx, end_idx, connected_to_next,
                        true, extended_bridges)) continue;
        }
      }

      // check for antiparallel bridge
      if(i > 0 && i < full_size - 1 && j > 0 && j < full_size - 1) {
        if((donor_for_one[j+1] == i - 1 && donor_for_one[i+1] == j - 1) ||
           (donor_for_one[j+1] == i - 1 && donor_for_two[i+1] == j - 1) ||
           (donor_for_two[j+1] == i - 1 && donor_for_one[i+1] == j - 1) ||
           (donor_for_two[j+1] == i - 1 && donor_for_two[i+1] == j - 1)) {
          if(AddIfValid(i,j,start_idx, end_idx, connected_to_next,
                        false, extended_bridges)) continue;
        }
      }
      if((donor_for_one[i] == j && donor_for_one[j] == i) ||
         (donor_for_one[i] == j && donor_for_two[j] == i) ||
         (donor_for_two[i] == j && donor_for_one[j] == i) ||
         (donor_for_two[i] == j && donor_for_two[j] == i)) {
        AddIfValid(i,j,start_idx, end_idx, connected_to_next, 
                   false, extended_bridges);
      }
    }
  }

  ResolveBulges(connected_to_next, extended_bridges);

  for(uint i = 0; i < extended_bridges.size(); ++i) {
    if(extended_bridges[i].strand_one.size() > 1) {
      for(int j = extended_bridges[i].strand_one.front(); 
          j <= extended_bridges[i].strand_one.back(); ++j) {
        int idx = j - start_idx;
        if(idx >= 0 && idx < size) ss_assignment[idx] = 'E';
      }
      if(extended_bridges[i].parallel) {
        for(int j = extended_bridges[i].strand_two.front(); 
            j <= extended_bridges[i].strand_two.back(); ++j) {
          int idx = j - start_idx;
          if(idx >= 0 && idx < size) ss_assignment[idx] = 'E';
        }
      }
      else {
        for(int j = extended_bridges[i].strand_two.back(); 
            j <= extended_bridges[i].strand_two.front(); ++j) {
          int idx = j - start_idx;
          if(idx >= 0 && idx < size) ss_assignment[idx] = 'E';
        }
      }
    }
    else {
      int idx = extended_bridges[i].strand_one[0] - start_idx;
      if(idx >= 0 && idx < size && ss_assignment[idx] != 'E') {
        ss_assignment[idx] = 'B';
      }
      idx = extended_bridges[i].strand_two[0] - start_idx;
      if(idx >= 0 && idx < size && ss_assignment[idx] != 'E') {
        ss_assignment[idx] = 'B';
      }
    }
  }
}

void AssignTurns(int start_idx, int size, 
                 const std::vector<geom::Vec3>& ca_positions,
	               const std::vector<int>& donor_for_one,
                 const std::vector<int>& donor_for_two, 
                 const std::vector<int>& connected_to_next,
                 String& ss_assignment) {

  int full_size = donor_for_one.size();

  // in turns we have following numbering scheme
  // 1: its the beginning of a turn
  // 2: its in the middle of a turn
  // 3: its the end of a turn
  // 4: its both, start and end

  std::vector<std::vector<int> > turn_arrays(3);
  turn_arrays[0].assign(full_size, 0);
  turn_arrays[1].assign(full_size, 0);
  turn_arrays[2].assign(full_size, 0);
  std::vector<int>& turn_3 = turn_arrays[0];
  std::vector<int>& turn_4 = turn_arrays[1];
  std::vector<int>& turn_5 = turn_arrays[2];

  // search for turns
  for(int turn_type = 3; turn_type <= 5; ++turn_type) {
    std::vector<int>& turn_array = turn_arrays[turn_type - 3];
    for(int i = std::max(0, start_idx - turn_type); 
        i < std::min(start_idx + size, full_size - turn_type); ++i) {

      if((donor_for_one[i + turn_type] == i ||
          donor_for_two[i + turn_type] == i) &&
          Connected(connected_to_next, i, turn_type))  {
            
        // last element is end of turn
        turn_array[i + turn_type] = 3;

        // set all the stuff in between to 2 (except it already is a start...)
        for(int j = i + 1; j < i + turn_type; ++j) {
          if(turn_array[j] == 0) turn_array[j] = 2;
        }

        // set first element to start (might be start AND end)
        if(turn_array[i] == 3) turn_array[i] = 4;
        else turn_array[i] = 1;
      }
    }
  }

  // to fill in all the turns in the relevant range
  // its easier to work in a separate array and then fill everything
  // into ss_assignment in the end to avoid indexing nightmare...
  //
  // coding:
  // 0: C
  // 1: T
  // 2: S
  // 3: G
  // 4: H
  // 5: I
  // 6: other stuff (all kind of bridges)
  std::vector<int> turn_array(full_size, 0);

  // first map the "other stuff" from the ss_assignment we already have
  for(int i = 0; i < size; ++i) {
    if(ss_assignment[i] != 'C') {
      turn_array[start_idx + i] = 6;
    }
  }

  // do H
  for(int i = std::max(start_idx - 3, 1); i < std::min(start_idx + size, full_size - 3); ++i) {
    if((turn_4[i] == 1 || turn_4[i] == 4) && 
       (turn_4[i - 1] == 1 || turn_4[i - 1] == 4)) {
      turn_array[i] = 4;
      turn_array[i+1] = 4;
      turn_array[i+2] = 4;
      turn_array[i+3] = 4;
    }
  }

  // do G
  for(int i = std::max(start_idx - 2, 1); i < std::min(start_idx + size, full_size - 2); ++i) {
    if((turn_3[i] == 1 || turn_3[i] == 4) && 
       (turn_3[i - 1] == 1 || turn_3[i - 1] == 4)) {
      // check whether we can still add the turn
      if((turn_array[i] == 0 || turn_array[i] == 3) &&
         (turn_array[i+1] == 0 || turn_array[i+1] == 3) &&
         (turn_array[i+2] == 0 || turn_array[i+2] == 3)) {
        turn_array[i] = 3;
        turn_array[i+1] = 3;
        turn_array[i+2] = 3;
      }
    }
  }

  // do I
  for(int i = std::max(start_idx - 4, 1); i < std::min(start_idx + size, full_size - 4); ++i) {
    if((turn_5[i] == 1 || turn_5[i] == 4) && 
       (turn_5[i - 1] == 1 || turn_5[i - 1] == 4)) {

      // check whether we can still add the turn
      bool add_it = true;
      for(int j = 0; j < 5; ++j) {
        if(!(turn_array[i+j] == 0 ||
             turn_array[i+j] == 5 ||
             turn_array[i+j] == 4)) { // we prefer pi helices!! 
                                        // (default in dssp)
          add_it = false;
          break;
        }
      }

      if(add_it) {
        turn_array[i] = 5;
        turn_array[i+1] = 5;
        turn_array[i+2] = 5;
        turn_array[i+3] = 5;
        turn_array[i+4] = 5;
      }
    }
  }

  // do turns(T) 

  // do turn_3
  for(int i = start_idx; i < start_idx + size; ++i) {
    if(turn_array[i] == 0) {
      if((i > 0 && (turn_3[i - 1] == 1 || turn_3[i - 1] == 4)) ||
         (i > 1 && (turn_3[i - 2] == 1 || turn_3[i - 2] == 4))) {
        turn_array[i] = 1;
      }
    }
  }

  // do turn_4
  for(int i = start_idx; i < start_idx + size; ++i) {
    if(turn_array[i] == 0) {
      if((i > 0 && (turn_4[i - 1] == 1 || turn_4[i - 1] == 4)) ||
         (i > 1 && (turn_4[i - 2] == 1 || turn_4[i - 2] == 4)) ||
         (i > 2 && (turn_4[i - 3] == 1 || turn_4[i - 3] == 4))) {
        turn_array[i] = 1;
      }
    }
  }

  // do turn_5
  for(int i = start_idx; i < start_idx + size; ++i) {
    if(turn_array[i] == 0) {
      if((i > 0 && (turn_5[i - 1] == 1 || turn_5[i - 1] == 4)) ||
         (i > 1 && (turn_5[i - 2] == 1 || turn_5[i - 2] == 4)) ||
         (i > 2 && (turn_5[i - 3] == 1 || turn_5[i - 3] == 4)) ||
         (i > 3 && (turn_5[i - 4] == 1 || turn_5[i - 4] == 4))) {
        turn_array[i] = 1;
      }
    }
  }

  // do bends (S)
  Real angle;
  for(int i = std::max(start_idx, 2); 
      i < std::min(start_idx + size, full_size - 2); ++i){
    if(turn_array[i] == 0) {
      angle = geom::Angle(ca_positions[i]-ca_positions[i-2],
                          ca_positions[i+2]-ca_positions[i]);
      if(angle > 1.2217 &&
         Connected(connected_to_next, i - 2, 5)) {
        turn_array[i] = 2;
      }
    }   
  }

  for(int i = start_idx; i < start_idx + size; ++i) {
    int idx = i - start_idx;
    if(turn_array[i] == 1) ss_assignment[idx] = 'T';
    else if(turn_array[i] == 2) ss_assignment[idx] = 'S';
    else if(turn_array[i] == 3) ss_assignment[idx] = 'G';
    else if(turn_array[i] == 4) ss_assignment[idx] = 'H';
    else if(turn_array[i] == 5) ss_assignment[idx] = 'I';
  }
}

} // anon ns

namespace ost{ namespace mol{ namespace alg{

String RawEstimateSS(const std::vector<geom::Vec3>& ca_positions, 
	                   int start_idx, int size,
	                   const std::vector<int>& donor_for_one,
                     const std::vector<int>& donor_for_two, 
                     const std::vector<int>& connected_to_next) {

  if(size < 1){
    throw ost::Error("Size of input is too small in dssp calculation!");
  }

  int full_size = ca_positions.size();

  if(start_idx + size > full_size) {
    throw ost::Error("Selected stretch is not fully covered by provided data in dssp calculation!");
  }

  if(full_size != static_cast<int>(donor_for_one.size())){
    throw ost::Error("Position and donor data are inconsistent for dssp calculation!");
  }

  if(full_size != static_cast<int>(donor_for_two.size())){
    throw ost::Error("Position and donor data are inconsistent for dssp calculation!");
  }

  String ss_assignment;
  ss_assignment.assign(size,'C');

  AssignBridges(start_idx, size, donor_for_one, donor_for_two,
                connected_to_next, ss_assignment);
                
  AssignTurns(start_idx, size, ca_positions, 
              donor_for_one, donor_for_two, 
              connected_to_next, ss_assignment);

  return ss_assignment;
}

void PrepareSSData(const ost::mol::ResidueViewList& res_list,
                   std::vector<int>& res_indices,
                   std::vector<geom::Vec3>& ca_positions,
                   std::vector<int>& donor_for_one,
                   std::vector<int>& donor_for_two,
                   std::vector<int>& connected_to_next) {

  res_indices.clear();
  ca_positions.clear();
  donor_for_one.clear();
  donor_for_two.clear();
  connected_to_next.clear();

  std::vector<geom::Vec3> n_positions;
  std::vector<geom::Vec3> c_positions;
  std::vector<geom::Vec3> o_positions;

  int total_num_residues = res_list.size();
  n_positions.reserve(total_num_residues);
  ca_positions.reserve(total_num_residues);
  c_positions.reserve(total_num_residues);
  o_positions.reserve(total_num_residues);
  res_indices.reserve(total_num_residues);
  connected_to_next.reserve(total_num_residues);

  for(int res_idx = 0; res_idx < total_num_residues; ++res_idx) {

    const ost::mol::ResidueView& res = res_list[res_idx];
    
    // add information if all atoms are present
    ost::mol::AtomView n  = res.FindAtom("N");
    ost::mol::AtomView ca = res.FindAtom("CA");
    ost::mol::AtomView c  = res.FindAtom("C");
    ost::mol::AtomView o  = res.FindAtom("O");
      
    if(n.IsValid() && ca.IsValid() && c.IsValid() && o.IsValid()) {
      n_positions.push_back(n.GetPos());
      ca_positions.push_back(ca.GetPos());
      c_positions.push_back(c.GetPos());
      o_positions.push_back(o.GetPos());
      res_indices.push_back(res_idx);
    }
  }

  int num_residues = res_indices.size();

  for(int i = 0; i < num_residues - 1; ++i) {

    if(res_list[res_indices[i]].GetChain() == 
       res_list[res_indices[i+1]].GetChain() &&
       geom::Length2(n_positions[i+1] - c_positions[i]) <= 6.25) {
      connected_to_next.push_back(1);
    }
    else {
      connected_to_next.push_back(0);
    }
  }

  // the last residue needs to be handled as well...
  connected_to_next.push_back(0);

  // estimate the hydrogen positions and also keep track which ones are valid
  // (e.g. in case of chain breaks / prolines...)
  // please note, that DSSP simply uses the directional vector given by
  // carbon and oxygen of the previous residue to determine the hydrogen
  // position
  std::vector<geom::Vec3> h_positions(num_residues);
  std::vector<int> valid_hydrogens(num_residues, 0);

  for(int i = 0; i < num_residues - 1; ++i) {
    if(connected_to_next[i] && res_list[res_indices[i+1]].GetName() != "PRO") {
      geom::Vec3 dir_vec = geom::Normalize(c_positions[i] - o_positions[i]);
      h_positions[i+1] = n_positions[i+1] + dir_vec;
      valid_hydrogens[i+1] = 1;
    }
  }

  std::pair<Real, Real> max_e_pair = 
                std::make_pair(std::numeric_limits<Real>::max(),
                               std::numeric_limits<Real>::max());

  std::vector<std::pair<Real,Real> > donor_energies(num_residues, max_e_pair);
  std::vector<std::pair<int,int> > acceptors(num_residues, std::make_pair(-1, -1));
  Real e;

  for(int i = 0; i < num_residues - 1; ++i){
    for(int j = i + 1; j < num_residues; ++j){
      if(geom::Length2(ca_positions[i] - ca_positions[j]) < 81.0){

        //i as acceptor, j as donor
        if(j > i + 1 && valid_hydrogens[j]){
          e = DSSPHBondEnergy(h_positions[j], n_positions[j],
                              c_positions[i], o_positions[i]);
          // stupid rounding to be consistent with dssp
          e = std::floor(e * 1000 + 0.5) * 0.001;

          if(e < donor_energies[j].first) {
            donor_energies[j].second = donor_energies[j].first;
            acceptors[j].second = acceptors[j].first;
            donor_energies[j].first = e;
            acceptors[j].first = i;
          }
          else if(e < donor_energies[j].second) {
            donor_energies[j].second = e;
            acceptors[j].second = i;
          }
        }

        //i as donor, j as acceptor
        if(valid_hydrogens[i]){

          e = DSSPHBondEnergy(h_positions[i], n_positions[i],
                              c_positions[j], o_positions[j]);
          // stupid rounding to be consistent with dssp
          e = std::floor(e * 1000 + 0.5) * 0.001;

          if(e < donor_energies[i].first) {
            donor_energies[i].second = donor_energies[i].first;
            acceptors[i].second = acceptors[i].first;
            donor_energies[i].first = e;
            acceptors[i].first = j;
          }
          else if(e < donor_energies[i].second) {
            donor_energies[i].second = e;
            acceptors[i].second = j; 
          }
        }
      }
    }
  } 

  donor_for_one.assign(num_residues, -1);
  donor_for_two.assign(num_residues, -1);

  for(int i = 0; i < num_residues; ++i) {
    if(donor_energies[i].first < -0.5) {
      donor_for_one[i] = acceptors[i].first;
    }
    if(donor_energies[i].second < -0.5) {
      donor_for_two[i] = acceptors[i].second;
    }
  }
}

void AssignSecStruct(ost::mol::EntityView& ent) {

  ost::mol::ResidueViewList res_list = ent.GetResidueList();

  for(uint res_idx = 0; res_idx < res_list.size(); ++res_idx) {
    // set all residues to coil to start with
    res_list[res_idx].SetSecStructure(ost::mol::SecStructure('C'));
  }

  std::vector<int> res_indices;
  std::vector<geom::Vec3> ca_positions;
  std::vector<int> donor_for_one;
  std::vector<int> donor_for_two;
  std::vector<int> connected_to_next;

  PrepareSSData(res_list, res_indices, ca_positions, 
                donor_for_one, donor_for_two, connected_to_next);

  String ss = RawEstimateSS(ca_positions, 0, ca_positions.size(),
                            donor_for_one, donor_for_two,
                            connected_to_next);

  for(uint i = 0; i < ca_positions.size(); ++i) {
    res_list[res_indices[i]].SetSecStructure(ost::mol::SecStructure(ss[i]));
  }
}

void AssignSecStruct(ost::mol::EntityHandle& ent) {
  ost::mol::EntityView view = ent.CreateFullView();
  AssignSecStruct(view);
}

}}} //ns
