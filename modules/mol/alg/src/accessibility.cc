#include <ost/mol/alg/accessibility.hh>

#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include <ost/mol/chain_view.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/atom_view.hh>

namespace{

struct Cube {
  
  Cube() {
    particle_indices.reserve(20);
  }

  void AddIndex(int idx) { particle_indices.push_back(idx); }
  const std::vector<int>& GetIndices() const { return particle_indices; };  

  std::vector<int> particle_indices;
};

struct CubeGrid {
  CubeGrid(Real cel, int x_c, int y_c, int z_c,
           Real x_min, Real y_min, Real z_min): cube_edge_length(cel),
                                                x_min(x_min),
                                                y_min(y_min),
                                                z_min(z_min) {

    x_cubes = std::max(1, x_c);
    y_cubes = std::max(1, y_c);
    z_cubes = std::max(1, z_c);
    int num_cubes = x_cubes * y_cubes * z_cubes;
    cubes = new Cube*[num_cubes];
    // assign NULL to each cube
    for(int i = 0; i < num_cubes; ++i) {
      cubes[i] = NULL;
    }
  }

  ~CubeGrid() {
    // cleanup cubes
    int num_cubes = x_cubes * y_cubes * z_cubes;
    for(int i = 0; i < num_cubes; ++i) {
      if(cubes[i] != NULL) delete cubes[i];
    }
    delete [] cubes;
  }

  static long int NumCubes(int x_c, int y_c, int z_c) {
    int x_cubes = std::max(1, x_c);
    int y_cubes = std::max(1, y_c);
    int z_cubes = std::max(1, z_c);
    long int n = x_cubes;
    n *= y_cubes;
    n *= z_cubes;
    return n;
  }

  int GetCubeIdx(Real x, Real y, Real z) {
    int x_cube = std::min(static_cast<int>((x - x_min) / cube_edge_length), 
                          x_cubes - 1);
    int y_cube = std::min(static_cast<int>((y - y_min) / cube_edge_length), 
                          y_cubes - 1);
    int z_cube = std::min(static_cast<int>((z - z_min) / cube_edge_length), 
                          z_cubes - 1);    
    return x_cube * y_cubes * z_cubes + y_cube * z_cubes + z_cube;
  }

  void AddIndex(Real x, Real y, Real z, int idx) {
    int cube_idx = this->GetCubeIdx(x, y, z);
    if(cubes[cube_idx] == NULL) cubes[cube_idx] = new Cube; 
    cubes[cube_idx]->AddIndex(idx);
  }

  void GetNeighbouringCubes(int cube_idx,
                            std::vector<Cube*>& neighbouring_cubes) const {

    int temp = cube_idx;
    int x_cube = temp / (y_cubes * z_cubes);
    temp -= x_cube * y_cubes * z_cubes;
    int y_cube = temp / z_cubes;
    temp -= y_cube * z_cubes;
    int z_cube =  temp;

    // let's try all neighbours... only add if they actually exist!
    int x_min = std::max(x_cube - 1, 0);
    int x_max = std::min(x_cube + 1, x_cubes - 1);
    int y_min = std::max(y_cube - 1, 0);
    int y_max = std::min(y_cube + 1, y_cubes - 1);
    int z_min = std::max(z_cube - 1, 0);
    int z_max = std::min(z_cube + 1, z_cubes - 1);
   
    neighbouring_cubes.clear();
    for(int i = x_min; i <= x_max; ++i) {
      for(int j = y_min; j <= y_max; ++j) {
        for(int k = z_min; k <= z_max; ++k){
          if(!(i == x_cube && j == y_cube && k == z_cube)){
            // its not the central cube, so lets add it if initialized
            cube_idx = i * y_cubes * z_cubes + j * z_cubes + k;
            if(cubes[cube_idx] != NULL) {
              neighbouring_cubes.push_back(cubes[cube_idx]);
            }
          }
        }
      }
    }
  }

  int GetNumCubes() const { return x_cubes*y_cubes*z_cubes; }

  bool IsInitialized(int cube_idx) const { return cubes[cube_idx] != NULL; }

  Cube* GetCube(int cube_idx) const { return cubes[cube_idx]; }

  Real cube_edge_length;
  int x_cubes;
  int y_cubes;
  int z_cubes;
  Real x_min;
  Real y_min;
  Real z_min;
  Cube** cubes;
};

// this struct might occur a bit stupid. The only thing it does is to keep track
// of occupied cubes with the exactly Lookup scheme than the CubeGrid. 
// It does not store anything else.
struct CubeOccupationGrid{

  CubeOccupationGrid(const CubeGrid& cube_grid): 
                     cube_edge_length(cube_grid.cube_edge_length),
                     x_cubes(cube_grid.x_cubes),
                     y_cubes(cube_grid.y_cubes),
                     z_cubes(cube_grid.z_cubes),
                     x_min(cube_grid.x_min),
                     y_min(cube_grid.y_min),
                     z_min(cube_grid.z_min) {
    int num_cubes = x_cubes * y_cubes * z_cubes;
    occupied_cubes = new bool[num_cubes];
    memset(occupied_cubes, 0, num_cubes * sizeof(bool));
  }

  ~CubeOccupationGrid() { delete [] occupied_cubes; }

  void Occupy(Real x_pos, Real y_pos, Real z_pos) {
    int x_cube = (x_pos - x_min) / cube_edge_length;
    int y_cube = (y_pos - y_min) / cube_edge_length;
    int z_cube = (z_pos - z_min) / cube_edge_length;
    int cube_idx = x_cube * y_cubes * z_cubes + y_cube * z_cubes + z_cube;
    occupied_cubes[cube_idx] = true;
  }

  bool IsOccupied(int cube_idx) {
    return occupied_cubes[cube_idx];
  }

  bool HasOccupiedNeighbour(int cube_idx) {
    int temp = cube_idx;
    int x_cube = temp / (y_cubes * z_cubes);
    temp -= x_cube * y_cubes * z_cubes;
    int y_cube = temp / z_cubes;
    temp -= y_cube * z_cubes;
    int z_cube =  temp;

    // let's try all neighbours... only add if they actually exist!
    int x_min = std::max(x_cube - 1, 0);
    int x_max = std::min(x_cube + 1, x_cubes - 1);
    int y_min = std::max(y_cube - 1, 0);
    int y_max = std::min(y_cube + 1, y_cubes - 1);
    int z_min = std::max(z_cube - 1, 0);
    int z_max = std::min(z_cube + 1, z_cubes - 1);
   
    for(int i = x_min; i <= x_max; ++i) {
      for(int j = y_min; j <= y_max; ++j) {
        for(int k = z_min; k <= z_max; ++k){
          if(!(i == x_cube && j == y_cube && k == z_cube)){
            // its not the central cube, so lets add it if initialized
            cube_idx = i * y_cubes * z_cubes + j * z_cubes + k;
            if(occupied_cubes[cube_idx]) return true; 
          }
        }
      }
    }
    return false;
  }

  Real cube_edge_length;
  int x_cubes;
  int y_cubes;
  int z_cubes;
  Real x_min;
  Real y_min;
  Real z_min;
  bool* occupied_cubes;
};

Real GetAtomAccessibilityNACCESS(Real x_pos, Real y_pos, Real z_pos,
                                 Real radius,
                                 const std::vector<Real>& x,
                                 const std::vector<Real>& y,
                                 const std::vector<Real>& z, 
                                 const std::vector<Real>& radii) {
                                  
  int num_close_atoms = x.size();

  if(num_close_atoms == 0) {
    // return area of full sphere
    return Real(4.0) * M_PI * radius * radius;
  } 

  std::vector<Real> dx(num_close_atoms);
  std::vector<Real> dy(num_close_atoms);
  std::vector<Real> dsqr(num_close_atoms);
  std::vector<Real> d(num_close_atoms);

  for(int i = 0; i < num_close_atoms; ++i) {
    Real a = x_pos - x[i];
    Real b = y_pos - y[i];  
    Real c = a*a + b*b;
    dx[i] = a;
    dy[i] = b;
    dsqr[i] = c;
    d[i] = std::sqrt(c);
  }

  Real area = 0.0;
  int num_z_slices = 20;
  Real z_res = Real(2.0) * radius / num_z_slices;
  Real z_grid = z_pos - radius + Real(0.5) * z_res;
  Real pi_two = Real(2.0) * M_PI;
  std::vector<std::pair<Real, Real> > arcs;

  for(int z_slice_idx = 0; z_slice_idx < num_z_slices; ++z_slice_idx) {
    
    arcs.clear();
    Real dz = z_grid - z_pos;

    // the radius of the circle given by the intersection of the sphere
    // with the current z plane (remember: r*r = x*x + y*y)
    Real z_slice_squared_r = radius*radius - dz*dz;

    Real z_slice_r = std::sqrt(z_slice_squared_r);
    bool fully_enclosed = false;

    for(int idx = 0; idx < num_close_atoms; ++idx) {

      Real close_atom_dz = z_grid - z[idx];
      Real close_atom_radius = radii[idx];
      

      // the radius of the circle given by the intersection of the sphere
      // around the close atom and the current z plane
      Real close_atom_z_slice_squared_r = close_atom_radius*close_atom_radius - 
                                           close_atom_dz*close_atom_dz;

      if(close_atom_z_slice_squared_r < 0.0) {
        continue;
      }

      Real close_atom_z_slice_r = std::sqrt(close_atom_z_slice_squared_r);

      Real summed_r = close_atom_z_slice_r + z_slice_r;

      if(d[idx] < summed_r) {
        // the atoms are close!

        Real r_diff = z_slice_r - close_atom_z_slice_r;

        if(d[idx] < std::abs(r_diff)) {
          // one circle is inside the other!
          if(r_diff <= Real(0.0)) {
            // the circle of the close atom fully encloses circle of currently
            // processed atom!
            fully_enclosed = true;
            break;
          } 
          else {
            // the circle of the currently processed atom encloses circle of the
            // close atom => simply skip close atom...
            continue;
          }
        }

        // estimate alpha: the angle between l1 and l2.
        // l1: line between the two circle centers
        // l2: line between circle center of reference atom and intersection point
        Real trig_test = (dsqr[idx] + z_slice_squared_r - 
                           close_atom_z_slice_squared_r) / 
                          (Real(2.0) * d[idx] * z_slice_r);
    
        if(trig_test >= 1.0) trig_test = Real(0.99999);
        if(trig_test <=-1.0) trig_test = Real(-0.99999);
        Real alpha = std::acos(trig_test);

        // estimate beta: the angle between l and the x-axis
        // l: line between the two circle centers
        Real beta = std::atan2(dy[idx], dx[idx]) + M_PI;

        // start and end of arc relative to x-axis
        Real arc_start = beta - alpha;
        Real arc_end = beta + alpha;
    
        // enforce range [0.0, 2*pi]
        if(arc_start < 0.0) arc_start += pi_two;
        if(arc_end > pi_two) arc_end -= pi_two;

        if(arc_end >= arc_start) {
          arcs.push_back(std::make_pair(arc_start, arc_end));
        }
        else {
          arcs.push_back(std::make_pair(arc_start, pi_two));
          arcs.push_back(std::make_pair(0.0, arc_end));
        }
      }
    }

    Real arc_sum = 0.0;

    if(!fully_enclosed) {

      if(arcs.empty()) {
        arc_sum = pi_two;
      }
      else{
        std::sort(arcs.begin(), arcs.end());
        arc_sum = arcs[0].first;
        Real end = arcs[0].second;
        for(uint i = 1; i < arcs.size(); ++i) {
          if(end < arcs[i].first) arc_sum += (arcs[i].first - end);
          if(arcs[i].second > end) end = arcs[i].second;
        }
        arc_sum += (pi_two - end);
      }
    }

    area += (arc_sum * z_res);
    z_grid += z_res;
  }

  // scale to VdW shell
  return area*radius;
}


struct DSSPHelper{
  
  DSSPHelper(Real x, Real y, Real z, Real r, Real d): x_pos(x), y_pos(y),
                                                      z_pos(z), 
                                                      squared_radius(r),
                                                      squared_distance(d) { }
  
  Real x_pos;
  Real y_pos;
  Real z_pos;
  Real squared_radius;
  Real squared_distance;
  
  bool operator<(const DSSPHelper& rhs) const { return squared_distance < 
                                                       rhs.squared_distance; }
};


Real GetAtomAccessibilityDSSP(Real x_pos, Real y_pos, Real z_pos,
                              Real radius,
                              const std::vector<Real>& x,
                              const std::vector<Real>& y,
                              const std::vector<Real>& z, 
                              const std::vector<Real>& radii) {

  int num_close_atoms = x.size();
  
  if(num_close_atoms == 0) {
    // return area of full sphere
    return Real(4.0) * M_PI * radius * radius;
  } 

  std::vector<DSSPHelper> helpers;
  helpers.reserve(num_close_atoms);

  Real dx, dy, dz, dist_squared, summed_radii;
  for(int i = 0; i < num_close_atoms; ++i) {
    dx = x_pos - x[i];
    dy = y_pos - y[i];
    dz = z_pos - z[i];
    dist_squared = dx*dx + dy*dy + dz*dz;
    summed_radii = radius + radii[i];
    summed_radii *= summed_radii;
    if(dist_squared < summed_radii) {
      helpers.push_back(DSSPHelper(x[i] - x_pos, y[i] - y_pos, z[i] - z_pos, 
                                   radii[i] * radii[i], dist_squared));
    }
  }

  std::sort(helpers.begin(), helpers.end());
  
  const ost::mol::alg::DSSPAccessibilityParam& param = 
  ost::mol::alg::DSSPAccessibilityParam::GetInstance();
  const std::vector<Real>& fibonacci_x = param.GetFibonacciX();
  const std::vector<Real>& fibonacci_y = param.GetFibonacciY();
  const std::vector<Real>& fibonacci_z = param.GetFibonacciZ();
  Real fibo_weight = param.GetPointWeight();
  Real summed_fibo_weight = 0.0;
  Real fibo_x, fibo_y, fibo_z;
  uint num_helpers = helpers.size();

  for(uint i = 0; i < fibonacci_x.size(); ++i) {
    fibo_x = fibonacci_x[i] * radius;
    fibo_y = fibonacci_y[i] * radius;
    fibo_z = fibonacci_z[i] * radius;

    bool covered = false;

    for(uint j = 0; j < num_helpers; ++j) {
      dx = fibo_x - helpers[j].x_pos;
      dy = fibo_y - helpers[j].y_pos;
      dz = fibo_z - helpers[j].z_pos;
      dist_squared = dx*dx + dy*dy + dz*dz;
      if(dist_squared < helpers[j].squared_radius) {
        covered = true;
        break;
      }
    }
    if(!covered) {
      summed_fibo_weight += fibo_weight;
    } 
  }

  return summed_fibo_weight * radius * radius;
}

Real GetAtomAccessibility(Real x_pos, Real y_pos, Real z_pos,
                          Real radius,
                          const std::vector<Real>& x,
                          const std::vector<Real>& y,
                          const std::vector<Real>& z, 
                          const std::vector<Real>& radii,
                          ost::mol::alg::AccessibilityAlgorithm algorithm) {
  switch(algorithm){
    case ost::mol::alg::NACCESS: return GetAtomAccessibilityNACCESS(x_pos, y_pos,
                                                                    z_pos, radius,
                                                                    x, y, z, radii);

    case ost::mol::alg::DSSP: return GetAtomAccessibilityDSSP(x_pos, y_pos,
                                                              z_pos, radius,
                                                              x, y, z, radii);

  }

  // silence compiler warning...
  return Real(0);
}


void SolveCube(const CubeGrid& cube_grid, int cube_idx,
               const std::vector<Real>& x, const std::vector<Real>& y, 
               const std::vector<Real>& z, const std::vector<Real>& radii, 
               ost::mol::alg::AccessibilityAlgorithm algorithm,
               std::vector<Real>& asa) {

  //prepare some stuff
  std::vector<Real> close_atom_x, close_atom_y, close_atom_z, close_atom_radii;

  close_atom_x.reserve(200);
  close_atom_y.reserve(200);
  close_atom_z.reserve(200);
  close_atom_radii.reserve(200);

  Cube* central_cube = cube_grid.GetCube(cube_idx);;
  std::vector<Cube*> neighbouring_cubes;
  cube_grid.GetNeighbouringCubes(cube_idx, neighbouring_cubes);
  const std::vector<int>& cube_atoms = central_cube->GetIndices();

  for(uint idx = 0; idx < cube_atoms.size(); ++idx) {

    int atom_idx = cube_atoms[idx];
    close_atom_x.clear();
    close_atom_y.clear();
    close_atom_z.clear();
    close_atom_radii.clear();

    Real current_x = x[atom_idx];
    Real current_y = y[atom_idx];
    Real current_z = z[atom_idx];
    Real current_radius = radii[atom_idx];

    // gather the stuff from all neighbouring cubes
    for(std::vector<Cube*>::iterator cube_it = neighbouring_cubes.begin();
        cube_it != neighbouring_cubes.end(); ++cube_it) {

      const std::vector<int>& cube_indices = (*cube_it)->GetIndices();

      for(std::vector<int>::const_iterator it = cube_indices.begin();
          it != cube_indices.end(); ++it) {
        int close_atom_idx = *it;
        close_atom_x.push_back(x[close_atom_idx]);
        close_atom_y.push_back(y[close_atom_idx]);
        close_atom_z.push_back(z[close_atom_idx]);
        close_atom_radii.push_back(radii[close_atom_idx]);
      }
    }

    // gather stuff from central cube
    const std::vector<int>& cube_indices = central_cube->GetIndices();
    for(std::vector<int>::const_iterator it = cube_indices.begin();
        it != cube_indices.end(); ++it) {

      int close_atom_idx = *it;

      if(atom_idx != close_atom_idx){
        close_atom_x.push_back(x[close_atom_idx]);
        close_atom_y.push_back(y[close_atom_idx]);
        close_atom_z.push_back(z[close_atom_idx]);
        close_atom_radii.push_back(radii[close_atom_idx]);
      }
    }
    
    // DOIT DOIT DOIT!!!
    asa[atom_idx] = GetAtomAccessibility(current_x, current_y, current_z, 
                                         current_radius, 
                                         close_atom_x, close_atom_y, 
                                         close_atom_z, close_atom_radii,
                                         algorithm);
  }
}


void AddProbeToRadii(const std::vector<Real>& radii,
                     Real probe_radius,
                     std::vector<Real>& full_radii) {
  int size = radii.size();
  full_radii.resize(size);
  for(int i = 0; i < size; ++i) {
    full_radii[i] = radii[i] + probe_radius;
  }
}


CubeGrid SetupCubeGrid(const std::vector<Real>& x,
                       const std::vector<Real>& y,
                       const std::vector<Real>& z,
                       const std::vector<Real>& full_radii) {

  Real x_min = x[0];
  Real x_max = x[0];
  Real y_min = y[0];
  Real y_max = y[0];
  Real z_min = z[0];
  Real z_max = z[0];
  Real r_max = full_radii[0];

  for(uint i = 1; i < x.size(); ++i){
    x_min = std::min(x_min, x[i]);
    x_max = std::max(x_max, x[i]);
    y_min = std::min(y_min, y[i]);
    y_max = std::max(y_max, y[i]);
    z_min = std::min(z_min, z[i]);
    z_max = std::max(z_max, z[i]);
    r_max = std::max(r_max, full_radii[i]);
  }
  
  Real x_range = x_max - x_min;
  Real y_range = y_max - y_min;
  Real z_range = z_max - z_min;
  Real cube_edge_length = 2 * r_max;
  int x_cubes = std::ceil(x_range / cube_edge_length);
  int y_cubes = std::ceil(y_range / cube_edge_length);
  int z_cubes = std::ceil(z_range / cube_edge_length);
  
  while(true) {
    long int n = CubeGrid::NumCubes(x_cubes, y_cubes, z_cubes);
    if(n <= 125000000) {
      // n is small enough, 125000000 corresponds to a grid of size 500x500x500
      break;
    }
    // we get too many cubes, let's increase edge length...
    cube_edge_length *= 1.1;
    x_cubes = std::ceil(x_range / cube_edge_length);
    y_cubes = std::ceil(y_range / cube_edge_length);
    z_cubes = std::ceil(z_range / cube_edge_length);
  }

  CubeGrid cube_grid(cube_edge_length, x_cubes, y_cubes, z_cubes,
                     x_min, y_min, z_min);

  return cube_grid;
}


void CalculateASAOligo(std::vector<Real>& x, std::vector<Real>& y, 
                       std::vector<Real>& z, std::vector<Real>& radii,
                       std::vector<int> chain_indices, Real probe_radius, 
                       ost::mol::alg::AccessibilityAlgorithm algorithm,
                       std::vector<Real>& asa,
                       std::vector<Real>& asa_single_chain) {

  int num_atoms = x.size();
  asa.resize(num_atoms);

  std::vector<Real> full_radii;
  AddProbeToRadii(radii, probe_radius, full_radii);

  CubeGrid cube_grid = SetupCubeGrid(x,y,z,full_radii);

  for(int i = 0; i < num_atoms; ++i) {
    cube_grid.AddIndex(x[i], y[i], z[i], i);
  }

  int num_cubes = cube_grid.GetNumCubes();
  for(int cube_idx = 0; cube_idx < num_cubes; ++cube_idx) {

    if(cube_grid.IsInitialized(cube_idx)) {
      // there is at least one atom!
      SolveCube(cube_grid, cube_idx, x, y, z, full_radii, 
                algorithm, asa);
    }
  }

  // let's start with the single chain stuff
  // the asa are assumed to be exactly the same for the majority of the atoms
  asa_single_chain = asa;

  // let's find all unique chains
  std::set<int> unique_indices;
  for(std::vector<int>::const_iterator it = chain_indices.begin(); 
      it != chain_indices.end(); ++it) {
    unique_indices.insert(*it);
  }

  for(std::set<int>::iterator set_it = unique_indices.begin(); 
      set_it != unique_indices.end(); ++set_it) {

    int current_chain_idx = *set_it;

    // let's generate a cube grid with exactly the same size as before 
    // but only fill it with stuff from the current chain
    CubeGrid single_chain_cube_grid(cube_grid.cube_edge_length, 
                                    cube_grid.x_cubes, 
                                    cube_grid.y_cubes, 
                                    cube_grid.z_cubes,
                                    cube_grid.x_min, 
                                    cube_grid.y_min, 
                                    cube_grid.z_min);

    // in here we store the cube grids that are occupied by the atoms NOT
    // belonging to the current chain
    CubeOccupationGrid occupation_grid(single_chain_cube_grid);

    for(int i = 0; i < num_atoms; ++i) {
      if(current_chain_idx == chain_indices[i]) {
        single_chain_cube_grid.AddIndex(x[i], y[i], z[i], i);
      }
      else {
        occupation_grid.Occupy(x[i], y[i], z[i]);
      }
    }

    for(int cube_idx = 0; cube_idx < single_chain_cube_grid.GetNumCubes(); 
        ++cube_idx) {

      if(single_chain_cube_grid.IsInitialized(cube_idx) && 
         (occupation_grid.IsOccupied(cube_idx) ||
          occupation_grid.HasOccupiedNeighbour(cube_idx))) {
        // there is at least one atom in the cube AND the environment has 
        // changes compared to the previous asa calculation
        SolveCube(single_chain_cube_grid, cube_idx, x, y, z, full_radii, 
                  algorithm, asa_single_chain);
      }
    }
  }
}


void CalculateASA(const std::vector<Real>& x, 
                  const std::vector<Real>& y, 
                  const std::vector<Real>& z, 
                  const std::vector<Real>& radii, 
                  Real probe_radius, 
                  ost::mol::alg::AccessibilityAlgorithm algorithm,
                  std::vector<Real>& asa) {

  int num_atoms = x.size();
  asa.resize(num_atoms);

  std::vector<Real> full_radii;
  AddProbeToRadii(radii, probe_radius, full_radii);

  CubeGrid cube_grid = SetupCubeGrid(x, y, z, full_radii);

  for(int i = 0; i < num_atoms; ++i) {
    cube_grid.AddIndex(x[i], y[i], z[i], i);
  }

  int num_cubes = cube_grid.GetNumCubes();
  for(int cube_idx = 0; cube_idx < num_cubes; ++cube_idx) {
    if(cube_grid.IsInitialized(cube_idx)) {
      // there is at least one atom!
      SolveCube(cube_grid, cube_idx, x, y, z, full_radii, 
                algorithm, asa);
    }
  }
}


void ASAParamFromAtomList(const ost::mol::AtomViewList& atom_list,
                          ost::mol::alg::AccessibilityAlgorithm algorithm,
                          std::vector<Real>& x_pos,
                          std::vector<Real>& y_pos,
                          std::vector<Real>& z_pos,
                          std::vector<Real>& radii) {

  if(algorithm == ost::mol::alg::NACCESS) {                          
    const ost::mol::alg::NACCESSAccessibilityParam& param = 
    ost::mol::alg::NACCESSAccessibilityParam::GetInstance();
    String rname, aname, ele;

    for(ost::mol::AtomViewList::const_iterator at_it = atom_list.begin();
        at_it != atom_list.end(); ++at_it) {
      rname = at_it->GetResidue().GetName();
      aname = at_it->GetName();
      ele = at_it->GetElement();
      geom::Vec3 at_pos = at_it->GetPos();
      x_pos.push_back(at_pos[0]);
      y_pos.push_back(at_pos[1]);
      z_pos.push_back(at_pos[2]);
      radii.push_back(param.GetVdWRadius(rname, aname, ele));
    } 
  } else if(algorithm == ost::mol::alg::DSSP) {
    const ost::mol::alg::DSSPAccessibilityParam& param = 
    ost::mol::alg::DSSPAccessibilityParam::GetInstance();
    String aname;

    for(ost::mol::AtomViewList::const_iterator at_it = atom_list.begin();
        at_it != atom_list.end(); ++at_it) {
      aname = at_it->GetName();
      geom::Vec3 at_pos = at_it->GetPos();
      x_pos.push_back(at_pos[0]);
      y_pos.push_back(at_pos[1]);
      z_pos.push_back(at_pos[2]);
      radii.push_back(param.GetVdWRadius(aname));
    } 
  }
}

void ChainIndicesFromAtomList(const ost::mol::AtomViewList& atom_list,
                              std::vector<int>& chain_indices) {

  // since no function to directly access a chain index, we have to do
  // that with an ugly hack... thats stupid and should be replaced...
  ost::mol::ChainViewList chain_list = 
  atom_list[0].GetEntity().GetChainList();
  ost::mol::ChainViewList::iterator chain_list_begin = chain_list.begin();
  ost::mol::ChainViewList::iterator chain_list_end = chain_list.end();

  for(ost::mol::AtomViewList::const_iterator at_it = atom_list.begin();
      at_it != atom_list.end(); ++at_it) {
    ost::mol::ChainView chain = at_it->GetResidue().GetChain();
    ost::mol::ChainViewList::iterator found_it = 
    std::find(chain_list_begin, chain_list_end, chain);
    chain_indices.push_back(static_cast<int>(found_it - chain_list_begin));
  } 
}  

Real SetAccessibilityProps(ost::mol::EntityView& ent,
                           ost::mol::AtomViewList& atom_list,
                           const std::vector<Real>& asa,
                           const String& asa_atom,
                           const String& asa_abs,
                           const String& asa_rel,
                           ost::mol::alg::AccessibilityAlgorithm algorithm) {

  // first set accessibilities of all involved residues to 0.0
  for(uint idx = 0; idx < asa.size(); ++idx) {
    ost::mol::ResidueView res = atom_list[idx].GetResidue();
    res.SetFloatProp(asa_abs, 0.0);
  }

  // assign absolute accessibilities
  Real summed_asa = 0.0;
  for(uint idx = 0; idx < asa.size(); ++idx) {
    Real val = asa[idx];
    atom_list[idx].SetFloatProp(asa_atom, val);
    ost::mol::ResidueView res = atom_list[idx].GetResidue();
    Real current_asa = res.GetFloatProp(asa_abs);
    res.SetFloatProp(asa_abs, current_asa + val);
    summed_asa += val;
  }

  String rname;
  ost::mol::ResidueViewList res_list = ent.GetResidueList();
  if(algorithm == ost::mol::alg::NACCESS) {
    // go over residue and assign relative accessibilities
    for(uint idx = 0; idx < res_list.size(); ++idx) {
      rname = res_list[idx].GetName();
      Real tot_acc = 
      ost::mol::alg::NACCESSAccessibilityParam::GetInstance().GetResidueAccessibility(rname);
      if(tot_acc == Real(-1.0)) {
        // no accessibility found... 
        // let's mimic NACCESS behaviour
        res_list[idx].SetFloatProp(asa_rel, Real(-99.9));
      }
      else {
        // the fraction gets multiplied by 100 (Thats how NACCESS does it...)
        Real val = res_list[idx].GetFloatProp(asa_abs, 0.0) / tot_acc * 100.0;
        res_list[idx].SetFloatProp(asa_rel, val);
      }
    }
  } else if(algorithm == ost::mol::alg::DSSP) {
    // go over residue and assign relative accessibilities
    for(uint idx = 0; idx < res_list.size(); ++idx) {
      rname = res_list[idx].GetName();
      Real tot_acc = 
      ost::mol::alg::DSSPAccessibilityParam::GetInstance().GetResidueAccessibility(rname);
      if(tot_acc == Real(-1.0)) {
        // no accessibility found... 
        // let's mimic DSSP binding behaviour
        // => do nothing
      }
      else {
        Real val = res_list[idx].GetFloatProp(asa_abs, 0.0) / tot_acc;
        res_list[idx].SetFloatProp(asa_rel, val);
      }
    }
  }

  return summed_asa;
}


} //ns



namespace ost { namespace mol { namespace alg {

NACCESSAccessibilityParam::NACCESSAccessibilityParam() {

  std::map<String, Real> ala_map;
  std::map<String, Real> arg_map;
  std::map<String, Real> asp_map;
  std::map<String, Real> asn_map;
  std::map<String, Real> cys_map;
  std::map<String, Real> glu_map;
  std::map<String, Real> gln_map;
  std::map<String, Real> gly_map;
  std::map<String, Real> his_map;
  std::map<String, Real> ile_map;
  std::map<String, Real> leu_map;
  std::map<String, Real> lys_map;
  std::map<String, Real> met_map;
  std::map<String, Real> phe_map;
  std::map<String, Real> pro_map;
  std::map<String, Real> ser_map;
  std::map<String, Real> thr_map;
  std::map<String, Real> trp_map;
  std::map<String, Real> tyr_map;
  std::map<String, Real> val_map; 
  std::map<String, Real> asx_map;
  std::map<String, Real> glx_map; 
  std::map<String, Real> ace_map;
  std::map<String, Real> __a_map;
  std::map<String, Real> __c_map;  
  std::map<String, Real> __g_map;
  std::map<String, Real> __t_map;
  std::map<String, Real> __u_map;
  std::map<String, Real> twomg_map;
  std::map<String, Real> h2u_map;
  std::map<String, Real> hem_map;
  std::map<String, Real> hoh_map;
 
  ala_map["N"]  = 1.65;  
  ala_map["CA"] = 1.87;  
  ala_map["C"]  = 1.76;  
  ala_map["O"]  = 1.40;  
  ala_map["CB"] = 1.87;

  arg_map["N"]  = 1.65; 
  arg_map["CA"] = 1.87; 
  arg_map["C"]  = 1.76; 
  arg_map["O"]  = 1.40; 
  arg_map["CB"] = 1.87; 
  arg_map["CG"] = 1.87; 
  arg_map["CD"] = 1.87; 
  arg_map["NE"] = 1.65; 
  arg_map["CZ"] = 1.76; 
  arg_map["NH1"]= 1.65; 
  arg_map["NH2"]= 1.65; 

  asp_map["N"]  = 1.65;
  asp_map["CA"] = 1.87;
  asp_map["C"]  = 1.76;
  asp_map["O"]  = 1.40;
  asp_map["CB"] = 1.87;
  asp_map["CG"] = 1.76;
  asp_map["OD1"]= 1.40;
  asp_map["OD2"]= 1.40;

  asn_map["N"]  = 1.65;  
  asn_map["CA"] = 1.87;  
  asn_map["C"]  = 1.76;  
  asn_map["O"]  = 1.40;  
  asn_map["CB"] = 1.87;  
  asn_map["CG"] = 1.76;  
  asn_map["OD1"]= 1.40;  
  asn_map["ND2"]= 1.65;  
         
  cys_map["N"]  = 1.65;  
  cys_map["CA"] = 1.87;  
  cys_map["C"]  = 1.76;  
  cys_map["O"]  = 1.40;  
  cys_map["CB"] = 1.87;  
  cys_map["SG"] = 1.85;  

  glu_map["N"]  = 1.65; 
  glu_map["CA"] = 1.87; 
  glu_map["C"]  = 1.76; 
  glu_map["O"]  = 1.40; 
  glu_map["CB"] = 1.87; 
  glu_map["CG"] = 1.87; 
  glu_map["CD"] = 1.76; 
  glu_map["OE1"]= 1.40; 
  glu_map["OE2"]= 1.40;
     
  gln_map["N"]  = 1.65; 
  gln_map["CA"] = 1.87; 
  gln_map["C"]  = 1.76; 
  gln_map["O"]  = 1.40; 
  gln_map["CB"] = 1.87; 
  gln_map["CG"] = 1.87; 
  gln_map["CD"] = 1.76; 
  gln_map["OE1"]= 1.40; 
  gln_map["NE2"]= 1.65; 

  gly_map["N"]  = 1.65;  
  gly_map["CA"] = 1.87;  
  gly_map["C"]  = 1.76;  
  gly_map["O"]  = 1.40;  

  his_map["N"]  = 1.65; 
  his_map["CA"] = 1.87; 
  his_map["C"]  = 1.76; 
  his_map["O"]  = 1.40; 
  his_map["CB"] = 1.87; 
  his_map["CG"] = 1.76; 
  his_map["ND1"]= 1.65; 
  his_map["CD2"]= 1.76; 
  his_map["CE1"]= 1.76; 
  his_map["NE2"]= 1.65; 
          
  ile_map["N"]  = 1.65; 
  ile_map["CA"] = 1.87; 
  ile_map["C"]  = 1.76;
  ile_map["O"]  = 1.40; 
  ile_map["CB"] = 1.87;
  ile_map["CG1"]= 1.87; 
  ile_map["CG2"]= 1.87; 
  ile_map["CD1"]= 1.87; 

  leu_map["N"]  = 1.65;
  leu_map["CA"] = 1.87; 
  leu_map["C"]  = 1.76; 
  leu_map["O"]  = 1.40; 
  leu_map["CB"] = 1.87; 
  leu_map["CG"] = 1.87; 
  leu_map["CD1"]= 1.87; 
  leu_map["CD2"]= 1.87; 

  lys_map["N"]  = 1.65; 
  lys_map["CA"] = 1.87; 
  lys_map["C"]  = 1.76; 
  lys_map["O"]  = 1.40; 
  lys_map["CB"] = 1.87; 
  lys_map["CG"] = 1.87; 
  lys_map["CD"] = 1.87; 
  lys_map["CE"] = 1.87; 
  lys_map["NZ"] = 1.50; 
           
  met_map["N"]  = 1.65; 
  met_map["CA"] = 1.87; 
  met_map["C"]  = 1.76; 
  met_map["O"]  = 1.40; 
  met_map["CB"] = 1.87; 
  met_map["CG"] = 1.87; 
  met_map["SD"] = 1.85; 
  met_map["CE"] = 1.87; 
            
  phe_map["N"]  = 1.65; 
  phe_map["CA"] = 1.87;
  phe_map["C"]  = 1.76; 
  phe_map["O"]  = 1.40; 
  phe_map["CB"] = 1.87; 
  phe_map["CG"] = 1.76; 
  phe_map["CD1"]= 1.76; 
  phe_map["CD2"]= 1.76; 
  phe_map["CE1"]= 1.76; 
  phe_map["CE2"]= 1.76; 
  phe_map["CZ"]= 1.76; 
        
  pro_map["N"]  = 1.65; 
  pro_map["CA"] = 1.87; 
  pro_map["C"]  = 1.76; 
  pro_map["O"]  = 1.40; 
  pro_map["CB"] = 1.87; 
  pro_map["CG"] = 1.87; 
  pro_map["CD"] = 1.87;

  ser_map["N"]  = 1.65; 
  ser_map["CA"] = 1.87; 
  ser_map["C"]  = 1.76; 
  ser_map["O"]  = 1.40; 
  ser_map["CB"] = 1.87; 
  ser_map["OG"] = 1.40; 

  thr_map["N"]  = 1.65; 
  thr_map["CA"] = 1.87; 
  thr_map["C"]  = 1.76; 
  thr_map["O"]  = 1.40; 
  thr_map["CB"] = 1.87; 
  thr_map["OG1"]= 1.40; 
  thr_map["CG2"]= 1.87; 

  trp_map["N"]  = 1.65; 
  trp_map["CA"] = 1.87; 
  trp_map["C"]  = 1.76; 
  trp_map["O"]  = 1.40; 
  trp_map["CB"] = 1.87; 
  trp_map["CG"] = 1.76; 
  trp_map["CD1"]= 1.76; 
  trp_map["CD2"]= 1.76; 
  trp_map["NE1"]= 1.65; 
  trp_map["CE2"]= 1.76; 
  trp_map["CE3"]= 1.76; 
  trp_map["CZ2"]= 1.76; 
  trp_map["CZ3"]= 1.76; 
  trp_map["CH2"]= 1.76; 

  tyr_map["N"]  = 1.65; 
  tyr_map["CA"] = 1.87; 
  tyr_map["C"]  = 1.76; 
  tyr_map["O"]  = 1.40; 
  tyr_map["CB"] = 1.87; 
  tyr_map["CG"] = 1.76; 
  tyr_map["CD1"]= 1.76; 
  tyr_map["CD2"]= 1.76; 
  tyr_map["CE1"]= 1.76; 
  tyr_map["CE2"]= 1.76; 
  tyr_map["CZ"]= 1.76; 
  tyr_map["OH"]= 1.40; 

  val_map["N"] = 1.65; 
  val_map["CA"] = 1.87; 
  val_map["C"] = 1.76; 
  val_map["O"] = 1.40; 
  val_map["CB"] = 1.87; 
  val_map["CG1"] = 1.87; 
  val_map["CG2"] = 1.87;
            
  asx_map["N"] = 1.65; 
  asx_map["CA"] = 1.87; 
  asx_map["C"] = 1.76; 
  asx_map["O"] = 1.40; 
  asx_map["CB"] = 1.87; 
  asx_map["CG"] = 1.76; 
  asx_map["AD1"] = 1.50; 
  asx_map["AD2"] = 1.50;
            
  glx_map["N"] = 1.65; 
  glx_map["CA"] = 1.87; 
  glx_map["C"] = 1.76; 
  glx_map["O"] = 1.40; 
  glx_map["CB"] = 1.87; 
  glx_map["CG"] = 1.76; 
  glx_map["CD"] = 1.87; 
  glx_map["AE1"] = 1.50; 
  glx_map["AE2"] = 1.50; 

  ace_map["C"] = 1.76; 
  ace_map["O"] = 1.40; 
  ace_map["CH3"] = 1.87; 

  __a_map["P"] = 1.90; 
  __a_map["O1P"] = 1.40; 
  __a_map["O2P"] = 1.40; 
  __a_map["O5*"] = 1.40; 
  __a_map["C5*"] = 1.80; 
  __a_map["C4*"] = 1.80; 
  __a_map["O4*"] = 1.40; 
  __a_map["C3*"] = 1.80; 
  __a_map["O3*"] = 1.40; 
  __a_map["C2*"] = 1.80; 
  __a_map["C1*"] = 1.80; 
  __a_map["N9"] = 1.60; 
  __a_map["C8"] = 1.80; 
  __a_map["N7"] = 1.60; 
  __a_map["C5"] = 1.80; 
  __a_map["C6"] = 1.80; 
  __a_map["N6"] = 1.60; 
  __a_map["N1"] = 1.60; 
  __a_map["C2"] = 1.80; 
  __a_map["N3"] = 1.60; 
  __a_map["C4"] = 1.80; 

  __c_map["P"] = 1.90; 
  __c_map["O1P"] = 1.40; 
  __c_map["O2P"] = 1.40; 
  __c_map["O5*"] = 1.40; 
  __c_map["C5*"] = 1.80; 
  __c_map["C4*"] = 1.80; 
  __c_map["O4*"] = 1.40; 
  __c_map["C3*"] = 1.80; 
  __c_map["O3*"] = 1.40; 
  __c_map["C2*"] = 1.80; 
  __c_map["C1*"] = 1.80; 
  __c_map["N1"] = 1.60; 
  __c_map["C2"] = 1.80; 
  __c_map["O2"] = 1.40; 
  __c_map["N3"] = 1.60; 
  __c_map["C4"] = 1.80; 
  __c_map["N4"] = 1.60; 
  __c_map["C5"] = 1.80; 
  __c_map["C6"] = 1.80; 

  __g_map["P"] = 1.90; 
  __g_map["O1P"] = 1.40; 
  __g_map["O2P"] = 1.40; 
  __g_map["O5*"] = 1.40; 
  __g_map["C5*"] = 1.80; 
  __g_map["C4*"] = 1.80; 
  __g_map["O4*"] = 1.40; 
  __g_map["C3*"] = 1.80; 
  __g_map["O3*"] = 1.40; 
  __g_map["C2*"] = 1.80; 
  __g_map["C1*"] = 1.80; 
  __g_map["N9"] = 1.60; 
  __g_map["C8"] = 1.80; 
  __g_map["N7"] = 1.60; 
  __g_map["C5"] = 1.80; 
  __g_map["C6"] = 1.80; 
  __g_map["O6"] = 1.40; 
  __g_map["N1"] = 1.60; 
  __g_map["C2"] = 1.80; 
  __g_map["N2"] = 1.60; 
  __g_map["N3"] = 1.60; 
  __g_map["C4"] = 1.80; 

  __t_map["P"] = 1.90; 
  __t_map["O1P"] = 1.40; 
  __t_map["O2P"] = 1.40; 
  __t_map["O5*"] = 1.40; 
  __t_map["C5*"] = 1.80; 
  __t_map["C4*"] = 1.80; 
  __t_map["O4*"] = 1.40; 
  __t_map["C3*"] = 1.80; 
  __t_map["O3*"] = 1.40; 
  __t_map["C2*"] = 1.80; 
  __t_map["C1*"] = 1.80; 
  __t_map["N1"] = 1.60; 
  __t_map["C2"] = 1.80; 
  __t_map["O2"] = 1.40; 
  __t_map["N3"] = 1.60; 
  __t_map["C4"] = 1.80; 
  __t_map["O4"] = 1.40; 
  __t_map["C5"] = 1.80; 
  __t_map["C5M"] = 1.80; 
  __t_map["C6"] = 1.80;

  __u_map["P"] = 1.90; 
  __u_map["O1P"] = 1.40; 
  __u_map["O2P"] = 1.40; 
  __u_map["O5*"] = 1.40; 
  __u_map["C5*"] = 1.80; 
  __u_map["C4*"] = 1.80; 
  __u_map["O4*"] = 1.40; 
  __u_map["C3*"] = 1.80; 
  __u_map["O3*"] = 1.40; 
  __u_map["C2*"] = 1.80; 
  __u_map["O2*"] = 1.40; 
  __u_map["C1*"] = 1.80; 
  __u_map["N1"] = 1.60; 
  __u_map["C2"] = 1.80; 
  __u_map["O2"] = 1.40; 
  __u_map["N3"] = 1.60; 
  __u_map["C4"] = 1.80; 
  __u_map["O4"] = 1.40; 
  __u_map["C5"] = 1.80; 
  __u_map["C6"] = 1.80;

  twomg_map["P"] = 1.90;
  twomg_map["O1P"] = 1.40; 
  twomg_map["O2P"] = 1.40; 
  twomg_map["O5*"] = 1.40; 
  twomg_map["C5*"] = 1.80; 
  twomg_map["C4*"] = 1.80; 
  twomg_map["O4*"] = 1.40; 
  twomg_map["C3*"] = 1.80; 
  twomg_map["O3*"] = 1.40; 
  twomg_map["C2*"] = 1.80; 
  twomg_map["O2*"] = 1.40; 
  twomg_map["C1*"] = 1.80; 
  twomg_map["N9"] = 1.60;
  twomg_map["C8"] = 1.80;
  twomg_map["N7"] = 1.60;
  twomg_map["C5"] = 1.80;
  twomg_map["C6"] = 1.80;
  twomg_map["O6"] = 1.40;
  twomg_map["N1"] = 1.60;
  twomg_map["C2"] = 1.80;
  twomg_map["N2"] = 1.60;
  twomg_map["C2A"] = 1.80; 
  twomg_map["N3"] = 1.60;
  twomg_map["C4"] = 1.80;

  h2u_map["P"] = 1.90; 
  h2u_map["O1P"] = 1.40; 
  h2u_map["O2P"] = 1.40; 
  h2u_map["O5*"] = 1.40; 
  h2u_map["C5*"] = 1.80; 
  h2u_map["C4*"] = 1.80; 
  h2u_map["O4*"] = 1.40; 
  h2u_map["C3*"] = 1.80; 
  h2u_map["O3*"] = 1.40; 
  h2u_map["C2*"] = 1.80; 
  h2u_map["O2*"] = 1.40; 
  h2u_map["C1*"] = 1.80; 
  h2u_map["N1"] = 1.60; 
  h2u_map["C2"] = 1.80; 
  h2u_map["O2"] = 1.40; 
  h2u_map["N3"] = 1.60; 
  h2u_map["C4"] = 1.80; 
  h2u_map["O4"] = 1.40; 
  h2u_map["C5"] = 1.80; 
  h2u_map["C6"] = 1.80; 

  hem_map["FE"] =  1.47; 
  hem_map["CHA"] = 2.00; 
  hem_map["CHB"] = 2.00; 
  hem_map["CHC"] = 2.00; 
  hem_map["CHD"] = 2.00; 
  hem_map["N A"] = 1.55; 
  hem_map["C1A"] = 1.78; 
  hem_map["C2A"] = 1.78; 
  hem_map["C3A"] = 1.78; 
  hem_map["C4A"] = 1.78; 
  hem_map["CMA"] = 1.90; 
  hem_map["CAA"] = 1.90; 
  hem_map["CBA"] = 1.90; 
  hem_map["CGA"] = 1.90; 
  hem_map["N B"] = 1.55; 
  hem_map["C1B"] = 1.78; 
  hem_map["C2B"] = 1.78; 
  hem_map["C3B"] = 1.78; 
  hem_map["C4B"] = 1.78; 
  hem_map["CMB"] = 1.90; 
  hem_map["CAB"] = 1.90; 
  hem_map["CBB"] = 1.90; 
  hem_map["N C"] = 1.55; 
  hem_map["C1C"] = 1.78; 
  hem_map["C2C"] = 1.78; 
  hem_map["C3C"] = 1.78; 
  hem_map["C4C"] = 1.78; 
  hem_map["CMC"] = 1.90; 
  hem_map["CAC"] = 1.90; 
  hem_map["CBC"] = 1.90; 
  hem_map["N D"] = 1.55; 
  hem_map["C1D"] = 1.78; 
  hem_map["C2D"] = 1.78; 
  hem_map["C3D"] = 1.78; 
  hem_map["C4D"] = 1.78; 
  hem_map["CMD"] = 1.90; 
  hem_map["CAD"] = 1.90; 
  hem_map["CBD"] = 1.90; 
  hem_map["CGD"] = 1.90; 
  hem_map["O1A"] = 1.35; 
  hem_map["O2A"] = 1.35; 
  hem_map["O1D"] = 1.35; 
  hem_map["O2D"] = 1.35; 

  hoh_map["O"] = 1.40; 


  vdw_radii_["ALA"] = ala_map;
  vdw_radii_["ARG"] = arg_map;
  vdw_radii_["ASP"] = asp_map;
  vdw_radii_["ASN"] = asn_map;
  vdw_radii_["CYS"] = cys_map;
  vdw_radii_["GLU"] = glu_map;
  vdw_radii_["GLN"] = gln_map;
  vdw_radii_["GLY"] = gly_map;
  vdw_radii_["HIS"] = his_map;
  vdw_radii_["ILE"] = ile_map;
  vdw_radii_["LEU"] = leu_map;
  vdw_radii_["LYS"] = lys_map;
  vdw_radii_["MET"] = met_map;
  vdw_radii_["PHE"] = phe_map;
  vdw_radii_["PRO"] = pro_map;
  vdw_radii_["SER"] = ser_map;
  vdw_radii_["THR"] = thr_map;
  vdw_radii_["TRP"] = trp_map;
  vdw_radii_["TYR"] = tyr_map;
  vdw_radii_["VAL"] = val_map;
  vdw_radii_["ASX"] = asx_map;
  vdw_radii_["GLX"] = glx_map;
  vdw_radii_["ACE"] = ace_map;
  vdw_radii_["__A"] = __a_map;
  vdw_radii_["__C"] = __c_map;
  vdw_radii_["__G"] = __g_map;
  vdw_radii_["__T"] = __t_map;
  vdw_radii_["__U"] = __u_map;
  vdw_radii_["2MG"] = twomg_map;
  vdw_radii_["H2U"] = h2u_map;
  vdw_radii_["HEM"] = hem_map;
  vdw_radii_["HOH"] = hoh_map;

  accessibilities_["ALA"] = 107.95;
  accessibilities_["CYS"] = 134.28;
  accessibilities_["ASP"] = 140.39;
  accessibilities_["GLU"] = 172.25;
  accessibilities_["PHE"] = 199.48;
  accessibilities_["GLY"] =  80.10;
  accessibilities_["HIS"] = 182.88;
  accessibilities_["ILE"] = 175.12;
  accessibilities_["LYS"] = 200.81;
  accessibilities_["LEU"] = 178.63;
  accessibilities_["MET"] = 194.15;
  accessibilities_["ASN"] = 143.94;
  accessibilities_["PRO"] = 136.13;
  accessibilities_["GLN"] = 178.50;
  accessibilities_["ARG"] = 238.76;
  accessibilities_["SER"] = 116.50;
  accessibilities_["THR"] = 139.27;
  accessibilities_["VAL"] = 151.44;
  accessibilities_["TRP"] = 249.36;
  accessibilities_["TYR"] = 212.76;
}


Real NACCESSAccessibilityParam::GuessRadius(const String& ele) const{
  if(ele == "C") return 1.80;
  if(ele == "N") return 1.60;
  if(ele == "S") return 1.85;
  if(ele == "O") return 1.40;
  if(ele == "P") return 1.90;
  if(ele == "CA") return 2.07;
  if(ele == "FE") return 1.47;
  if(ele == "CU") return 1.78;
  if(ele == "ZN") return 1.39;
  if(ele == "MG") return 1.73;
  if(ele == "H") return 1.0;
  if(ele == "D") return 1.0;
  return 1.80;
}


Real NACCESSAccessibilityParam::GetVdWRadius(const String& rname, 
                                             const String& aname, 
                                             const String& ele) const{
  std::map<String, std::map<String, Real> >::const_iterator res_map_it =
  vdw_radii_.find(rname);
  if(res_map_it != vdw_radii_.end()) {
    std::map<String, Real>::const_iterator at_it = 
    res_map_it->second.find(aname);
    if(at_it != res_map_it->second.end()) return at_it->second;
  }
  return GuessRadius(ele);
}


Real NACCESSAccessibilityParam::GetResidueAccessibility(const String& rname) const{
  std::map<String, Real>::const_iterator it = accessibilities_.find(rname);
  if(it != accessibilities_.end()) return it->second;
  else return Real(-1.0);
}

DSSPAccessibilityParam::DSSPAccessibilityParam() {

  // these values are from the dssp binding. nobody knows
  // the source...
  accessibilities_["ALA"] = Real(118.0);
  accessibilities_["ARG"] = Real(317.0);
  accessibilities_["ASN"] = Real(238.0);
  accessibilities_["ASP"] = Real(243.0);
  accessibilities_["CYS"] = Real(183.0);
  accessibilities_["GLN"] = Real(262.0);
  accessibilities_["GLU"] = Real(286.0);
  accessibilities_["GLY"] = Real(154.0);
  accessibilities_["HIS"] = Real(258.0);
  accessibilities_["ILE"] = Real(228.0);
  accessibilities_["LEU"] = Real(243.0);
  accessibilities_["LYS"] = Real(278.0);
  accessibilities_["MET"] = Real(260.0);
  accessibilities_["PHE"] = Real(271.0);
  accessibilities_["PRO"] = Real(204.0);
  accessibilities_["SER"] = Real(234.0);
  accessibilities_["THR"] = Real(206.0);
  accessibilities_["TRP"] = Real(300.0);
  accessibilities_["TYR"] = Real(303.0);
  accessibilities_["VAL"] = Real(216.0);


  // construct evenly distributed points on a sphere.
  // same as in DSSP...
  fibonacci_x_.reserve(401);
  fibonacci_y_.reserve(401);
  fibonacci_z_.reserve(401);

  Real golden_ratio = (Real(1.0) + std::sqrt(5.0)) * Real(0.5); 
	point_weight_ = (Real(4.0) * Real(M_PI)) / Real(401.0);

	for (int i = -200; i <= 200; ++i) {
		Real lat = std::asin((Real(2.0) * Real(i)) / Real(401.0));
    Real lon = std::fmod(Real(i), golden_ratio) * 
               Real(2.0) * Real(M_PI) / golden_ratio;
    fibonacci_x_.push_back(std::sin(lon) * std::cos(lat));
    fibonacci_y_.push_back(std::cos(lon) * std::cos(lat));
    fibonacci_z_.push_back(std::sin(lat));
	}
}

Real DSSPAccessibilityParam::GetVdWRadius(const String& aname) const{
  if(aname == "N") return Real(1.65);
  if(aname == "CA") return Real(1.87);
  if(aname == "C") return Real(1.76);
  if(aname == "O") return Real(1.4);
  return Real(1.8); 
}

Real DSSPAccessibilityParam::GetResidueAccessibility(const String& rname) const{
  std::map<String, Real>::const_iterator it = accessibilities_.find(rname);
  if(it != accessibilities_.end()) return it->second;
  else return Real(-1.0);
}

Real Accessibility(ost::mol::EntityView& ent, 
                   Real probe_radius, bool include_hydrogens,
                   bool include_hetatm, bool include_water,
                   bool oligo_mode,
                   const String& selection, 
                   const String& asa_abs, 
                   const String& asa_rel,
                   const String& asa_atom,
                   AccessibilityAlgorithm algorithm) {

  String internal_selection = selection;

  if(!include_hydrogens) {
    if(internal_selection == "") {
      internal_selection = "ele!=H,D";
    }
    else {
      internal_selection += " and ele!=H,D";
    }
  }

  if(!include_hetatm) {
    if(internal_selection == "") {
      internal_selection = "ishetatm=false";
    }
    else {
      internal_selection += " and ishetatm=false";
    }
  }

  if(!include_water) {
    if(internal_selection == "") {
      internal_selection = "water=false";
    }
    else {
      internal_selection += " and water=false";
    }
  }

  ost::mol::EntityView selected_ent;
  if(internal_selection != "") { 
    selected_ent = ent.Select(internal_selection);
  }    
  else{
    selected_ent = ent;
  }

  // input for algorithm
  std::vector<Real> x_pos;
  std::vector<Real> y_pos;
  std::vector<Real> z_pos;
  std::vector<Real> radii;

  x_pos.reserve(selected_ent.GetAtomCount());
  y_pos.reserve(selected_ent.GetAtomCount());
  z_pos.reserve(selected_ent.GetAtomCount());
  radii.reserve(selected_ent.GetAtomCount());

  // extract data from ent
  ost::mol::AtomViewList atom_list = selected_ent.GetAtomList();
  ASAParamFromAtomList(atom_list, algorithm, x_pos, y_pos, z_pos, radii);

  if(atom_list.size() == 0) return 0.0;

  if(oligo_mode) {

    // we additionally need the chain index for every atom
    std::vector<int> chain_indices;
    chain_indices.reserve(selected_ent.GetAtomCount());
    ChainIndicesFromAtomList(atom_list, chain_indices);

    // do it! do it! do it!              
    std::vector<Real> asa;
    std::vector<Real> asa_single_chain;
    CalculateASAOligo(x_pos, y_pos, z_pos, radii, chain_indices, 
                      probe_radius, algorithm, asa, asa_single_chain);

    Real summed_asa = SetAccessibilityProps(selected_ent, atom_list, asa, 
                                            asa_atom, asa_abs, asa_rel, 
                                            algorithm);

    // in case of the single chain asa, the property names are
    // moodified!!
    String single_chain_asa_atom = asa_atom + "_single_chain";
    String single_chain_asa_abs = asa_abs + "_single_chain";
    String single_chain_asa_rel = asa_rel + "_single_chain";

    SetAccessibilityProps(selected_ent, atom_list, asa_single_chain, 
                          single_chain_asa_atom, 
                          single_chain_asa_abs, 
                          single_chain_asa_rel, algorithm);

    return summed_asa;
  }
  else {
    
    // do it! do it! do it!
    std::vector<Real> asa;
    CalculateASA(x_pos, y_pos, z_pos, radii, probe_radius, algorithm, asa);

    Real summed_asa = SetAccessibilityProps(selected_ent, atom_list, asa, 
                                            asa_atom, asa_abs, asa_rel,
                                            algorithm);

    return summed_asa;
  }
}


Real Accessibility(ost::mol::EntityHandle& ent, 
                   Real probe_radius, bool include_hydrogens,
                   bool include_hetatm, bool include_water,
                   bool oligo_mode,
                   const String& selection,
                   const String& asa_abs, 
                   const String& asa_rel,
                   const String& asa_atom,
                   AccessibilityAlgorithm algorithm) {

  ost::mol::EntityView entity_view = ent.CreateFullView();
  return Accessibility(entity_view, probe_radius, include_hydrogens,
                       include_hetatm, include_water, oligo_mode,
                       selection, asa_abs, asa_rel, asa_atom, algorithm);
}

}}} //ns
