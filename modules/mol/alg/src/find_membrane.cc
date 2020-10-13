#include <ost/mol/alg/find_membrane.hh>
#include <ost/mol/alg/accessibility.hh>
#include <ost/geom/vecmat3_op.hh>
#include <ost/message.hh>

#include <limits>
#include <exception>
#include <list>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>

namespace{

// Copyright notice of the Levenberg Marquardt minimizer we use...

// Copyright (c) 2007, 2008, 2009 libmv authors.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
// IN THE SOFTWARE.
//
// A simple implementation of levenberg marquardt.
//
// [1] K. Madsen, H. Nielsen, O. Tingleoff. Methods for Non-linear Least
// Squares Problems.
// http://www2.imm.dtu.dk/pubdb/views/edoc_download.php/3215/pdf/imm3215.pdf
//
// TODO(keir): Cite the Lourakis' dogleg paper.

template<typename Function,
         typename Jacobian,
         typename Solver = Eigen::JacobiSVD<
                  Eigen::Matrix<typename Function::FMatrixType::RealScalar,
                  Function::XMatrixType::RowsAtCompileTime,
                  Function::XMatrixType::RowsAtCompileTime> > >
class LevenbergMarquardt {
 public:
  typedef typename Function::XMatrixType::RealScalar Scalar;
  typedef typename Function::FMatrixType FVec;
  typedef typename Function::XMatrixType Parameters;
  typedef Eigen::Matrix<typename Function::FMatrixType::RealScalar,
                        Function::FMatrixType::RowsAtCompileTime,
                        Function::XMatrixType::RowsAtCompileTime> JMatrixType;
  typedef Eigen::Matrix<typename JMatrixType::RealScalar,
                        JMatrixType::ColsAtCompileTime,
                        JMatrixType::ColsAtCompileTime> AMatrixType;

  // TODO(keir): Some of these knobs can be derived from each other and
  // removed, instead of requiring the user to set them.
  enum Status {
    RUNNING,
    GRADIENT_TOO_SMALL,            // eps > max(J'*f(x))
    RELATIVE_STEP_SIZE_TOO_SMALL,  // eps > ||dx|| / ||x||
    ERROR_TOO_SMALL,               // eps > ||f(x)||
    HIT_MAX_ITERATIONS,
  };

  LevenbergMarquardt(const Function &f)
      : f_(f), df_(f) {}

  struct SolverParameters {
   SolverParameters()
       : gradient_threshold(1e-20),
         relative_step_threshold(1e-20),
         error_threshold(1e-16),
         initial_scale_factor(1e-3),
         max_iterations(200) {}
    Scalar gradient_threshold;       // eps > max(J'*f(x))
    Scalar relative_step_threshold;  // eps > ||dx|| / ||x||
    Scalar error_threshold;          // eps > ||f(x)||
    Scalar initial_scale_factor;     // Initial u for solving normal equations.
    int    max_iterations;           // Maximum number of solver iterations.
  };

  struct Results {
    Scalar error_magnitude;     // ||f(x)||
    Scalar gradient_magnitude;  // ||J'f(x)||
    int    iterations;
    Status status;
  };

  Status Update(const Parameters &x, const SolverParameters &params,
                JMatrixType *J, AMatrixType *A, FVec *error, Parameters *g) {
    *J = df_(x);
    *A = (*J).transpose() * (*J);
    *error = -f_(x);
    *g = (*J).transpose() * *error;
    if (g->array().abs().maxCoeff() < params.gradient_threshold) {
      return GRADIENT_TOO_SMALL;
    } else if (error->norm() < params.error_threshold) {
      return ERROR_TOO_SMALL;
    }
    return RUNNING;
  }

  Results minimize(Parameters *x_and_min) {
    SolverParameters params;
    return minimize(params, x_and_min);
  }

  Results minimize(const SolverParameters &params, Parameters *x_and_min) {
    Parameters &x = *x_and_min;
    JMatrixType J;
    AMatrixType A;
    FVec error;
    Parameters g;

    Results results;
    results.status = Update(x, params, &J, &A, &error, &g);

    Scalar u = Scalar(params.initial_scale_factor*A.diagonal().maxCoeff());
    Scalar v = 2;

    Parameters dx, x_new;
    int i;
    for (i = 0; results.status == RUNNING && i < params.max_iterations; ++i) {
//      LOG(INFO) << "iteration: " << i;
//      LOG(INFO) << "||f(x)||: " << f_(x).norm();
//      LOG(INFO) << "max(g): " << g.array().abs().maxCoeff();
//      LOG(INFO) << "u: " << u;
//      LOG(INFO) << "v: " << v;
      AMatrixType A_augmented = A + u*AMatrixType::Identity(J.cols(), J.cols());
      Solver solver(A_augmented, Eigen::ComputeThinU | Eigen::ComputeThinV);
      dx = solver.solve(g);
      if (dx.norm() <= params.relative_step_threshold * x.norm()) {
          results.status = RELATIVE_STEP_SIZE_TOO_SMALL;
          break;
      }

      x_new = x + dx;
      // Rho is the ratio of the actual reduction in error to the reduction
      // in error that would be obtained if the problem was linear.
      // See [1] for details.
      Scalar rho((error.squaredNorm() - f_(x_new).squaredNorm())
                 / dx.dot(u*dx + g));
      if (rho > 0) {
        // Accept the Gauss-Newton step because the linear model fits well.
        x = x_new;
        results.status = Update(x, params, &J, &A, &error, &g);
        Scalar tmp = Scalar(2*rho-1);
        u = u*std::max(Scalar(1/3.), 1 - (tmp*tmp*tmp));
        v = 2;
        continue;
      }
      
      // Reject the update because either the normal equations failed to solve
      // or the local linear model was not good (rho < 0). Instead, increase 
      // to move closer to gradient descent.
      u *= v;
      v *= 2;
    }
    if (results.status == RUNNING) {
      results.status = HIT_MAX_ITERATIONS;
    }
    results.error_magnitude = error.norm();
    results.gradient_magnitude = g.norm();
    results.iterations = i;
    return results;
  }

 private:
  const Function &f_;
  Jacobian df_;
};

geom::Mat3 RotationAroundAxis(geom::Vec3 axis, Real angle) {

  Real aa, ab, ac, ba, bb, bc, ca, cb, cc, one_m_cos, cos_ang, sin_ang;

  cos_ang = std::cos(angle);
  sin_ang = std::sin(angle);
  one_m_cos = 1-cos_ang;

  aa = cos_ang+axis[0]*axis[0]*one_m_cos;
  ab = axis[0]*axis[1]*one_m_cos-axis[2]*sin_ang;
  ac = axis[0]*axis[2]*one_m_cos+axis[1]*sin_ang;

  ba = axis[1]*axis[0]*one_m_cos+axis[2]*sin_ang;
  bb = cos_ang+axis[1]*axis[1]*one_m_cos;
  bc = axis[1]*axis[2]*one_m_cos-axis[0]*sin_ang;

  ca = axis[2]*axis[0]*one_m_cos-axis[1]*sin_ang;
  cb = axis[2]*axis[1]*one_m_cos+axis[0]*sin_ang;
  cc = cos_ang+axis[2]*axis[2]*one_m_cos;

  geom::Mat3 result(aa, ab, ac, ba, bb, bc, ca, cb, cc);
  return result;
}

geom::Vec3 RotateAroundAxis(geom::Vec3 point, geom::Vec3 axis, Real angle) {

  geom::Mat3 rot = RotationAroundAxis(axis, angle);
  geom::Vec3 result = rot*point;
  return result;
}


// Levenberg Marquardt specific objects for the membrane finding algorithm

struct EnergyF {
  EnergyF(const std::vector<geom::Vec3>& p, const std::vector<Real>& t_e, 
          Real l, Real o, const geom::Vec3& ax, const geom::Vec3& tilt_ax):
          positions(p), transfer_energies(t_e), axis(ax),
          tilt_axis(tilt_ax), lambda(l), offset(o) { }

  typedef Eigen::Matrix<Real, 4, 1> XMatrixType;
  typedef Eigen::Matrix<Real, 1, 1> FMatrixType;

  Eigen::Matrix<Real,1,1> operator()(const Eigen::Matrix<Real, 4, 1>& x) const {

    FMatrixType result;
    result(0,0) = 0.0;
    geom::Vec3 tilted_axis = axis;

    tilted_axis = RotateAroundAxis(tilted_axis, tilt_axis, x(0, 0));
    tilted_axis = RotateAroundAxis(tilted_axis, axis, x(1, 0));

    Real pos_on_axis;
    Real distance_to_center;
    Real half_width = Real(0.5) * x(2, 0);
    Real exponent;
    Real one_over_lambda = Real(1.0) / lambda;
  
    int n = transfer_energies.size();

    for(int i = 0; i < n; ++i) {
      pos_on_axis = geom::Dot(tilted_axis, positions[i]);
      distance_to_center = std::abs(x(3, 0) - pos_on_axis);
      exponent = (distance_to_center - half_width) * one_over_lambda;
      result(0, 0) += (1.0 / (1.0 + std::exp(exponent))) * transfer_energies[i];    
    }

    result(0, 0) += offset;

    return result;
  }

  std::vector<geom::Vec3> positions;
  std::vector<Real> transfer_energies;
  geom::Vec3 axis;
  geom::Vec3 tilt_axis;
  Real lambda;

  // define an offset parameter...
  // The levenberg-marquardt algorithm is designed to minimize an error 
  // function, aims to converge towards zero. The offset parameter gets added 
  // at the end of the energy calculation to get a positive result => 
  // forces algorithm to minimize
  Real offset;
};


struct EnergyDF {

  EnergyDF(const EnergyF& f): function(f),d_tilt(0.02),d_angle(0.02),
                         d_width(0.4),d_pos(0.4) { }

  Eigen::Matrix<Real,1,4> operator()(const Eigen::Matrix<Real, 4, 1>& x) const {

    Eigen::Matrix<Real,1,4> result;
    Eigen::Matrix<Real, 4, 1> parameter1 = x;
    Eigen::Matrix<Real, 4, 1> parameter2 = x;

    parameter1(0,0)+=d_tilt;
    parameter2(0,0)-=d_tilt;
    result(0,0) = (function(parameter1)(0, 0) -
                   function(parameter2)(0, 0)) / (2*d_tilt);

    parameter1=x;
    parameter2=x;
    parameter1(1,0)+=d_angle;
    parameter2(1,0)-=d_angle;
    result(0,1) = (function(parameter1)(0,0) -
                   function(parameter2)(0,0)) / (2*d_angle);

    parameter1=x;
    parameter2=x;
    parameter1(2,0)+=d_width;
    parameter2(2,0)-=d_width;
    result(0,2) = (function(parameter1)(0,0) - 
                   function(parameter2)(0,0)) / (2*d_width);

    parameter1=x;
    parameter2=x;
    parameter1(3,0)+=d_pos;
    parameter2(3,0)-=d_pos;
    result(0,3) = (function(parameter1)(0,0) - 
                   function(parameter2)(0,0)) / (2*d_pos);
    return result;
  }

  EnergyF function;
  Real d_tilt;
  Real d_angle;
  Real d_width;
  Real d_pos;
};


void GetRanges(const std::vector<geom::Vec3>& atom_positions,
               Real& min_x, Real& max_x, Real& min_y, Real& max_y,
               Real& min_z, Real& max_z) {

  min_x = std::numeric_limits<Real>::max();
  max_x = -min_x;

  min_y = std::numeric_limits<Real>::max();
  max_y = -min_y;

  min_z = std::numeric_limits<Real>::max();
  max_z = -min_z;

  for(uint i = 0; i < atom_positions.size(); ++i) {
    const geom::Vec3& pos = atom_positions[i];
    min_x = std::min(min_x, pos[0]);
    max_x = std::max(max_x, pos[0]);
    min_y = std::min(min_y, pos[1]);
    max_y = std::max(max_y, pos[1]);
    min_z = std::min(min_z, pos[2]);
    max_z = std::max(max_z, pos[2]);
  }
}


void FloodLevel(char* data, int x_start, int y_start, 
                int x_extent, int y_extent, 
                int orig_value, int dest_value) {

  //http://lodev.org/cgtutor/floodfill.html
  if(orig_value != data[x_start*y_extent + y_start]) {
    return;
  }

  std::vector<std::pair<int,int> > queue;
  queue.push_back(std::make_pair(x_start, y_start));

  int y1,y,x; 
  bool spanLeft, spanRight;
  std::pair<int,int> actual_position;
    
  while(!queue.empty()){    

    actual_position = queue.back();
    queue.pop_back();
    x = actual_position.first;
    y = actual_position.second;

    y1 = y;

    while(y1 >= 0 && data[x*y_extent + y1] == orig_value) {
      --y1;
    }

    y1++;
    spanLeft = spanRight = 0;

    while(y1 < y_extent && 
          data[x*y_extent + y1] == orig_value ) {

      data[x*y_extent + y1] = dest_value;

      if(!spanLeft && x > 0 && 
         data[(x-1)*y_extent + y1] == orig_value) {
        queue.push_back(std::make_pair(x-1,y1));
        spanLeft = 1;
      }
      else if(spanLeft && x > 0 && 
              data[(x-1)*y_extent + y1] != orig_value) {
        spanLeft = 0;
      }

      if(!spanRight && x < x_extent - 1 && 
         data[(x+1)*y_extent + y1] == orig_value) {
        queue.push_back(std::make_pair(x+1,y1));
        spanRight = 1;
      }
      else if(spanRight && x < x_extent - 1 && 
              data[(x+1)*y_extent + y1] != orig_value) {
        spanRight = 0;
      } 
      ++y1;
    }
  }
}


void GetExposedAtoms(const std::vector<geom::Vec3>& atom_positions,
                     const std::vector<Real>& transfer_energies,
                     const std::vector<Real>& asas,
                     std::vector<geom::Vec3>& exposed_atom_positions,
                     std::vector<Real>& exposed_transfer_energies,
                     std::vector<Real>& exposed_asas) {

  // sum of approx. vdw radius of the present heavy atoms (1.8) 
  // plus 1.4 for water.
  Real radius = 3.2;
  Real one_over_radius = Real(1.0) / radius;

  // lets setup a grid in which we place the atoms
  Real min_x, max_x, min_y, max_y, min_z, max_z;
  GetRanges(atom_positions, min_x, max_x, min_y, max_y, min_z, max_z);

  // we guarantee that the thing is properly solvated in the x-y plane and add 
  // some space around this is also necessary to avoid overflow checks in 
  // different places
  min_x -= Real(2.1) * radius;
  min_y -= Real(2.1) * radius;
  min_z -= Real(2.1) * radius;
  max_x += Real(2.1) * radius;
  max_y += Real(2.1) * radius;
  max_z += Real(2.1) * radius;

  int num_xbins = std::ceil((max_x - min_x) * one_over_radius);
  int num_ybins = std::ceil((max_y - min_y) * one_over_radius);
  int num_zbins = std::ceil((max_z - min_z) * one_over_radius);
  int num_bins = num_xbins * num_ybins * num_zbins;
  char* grid = new char[num_bins];
  memset(grid, 0, num_bins);

  for(uint i = 0; i < atom_positions.size(); ++i) {

    const geom::Vec3& pos = atom_positions[i];
    int x_bin = (pos[0] - min_x) * one_over_radius;
    int y_bin = (pos[1] - min_y) * one_over_radius;
    int z_bin = (pos[2] - min_z) * one_over_radius;

    // we're really crude here and simply set all 27 cubes with central
    // cube defined by x_bin, y_bin and z_bin to one
    for(int z = z_bin - 1; z <= z_bin + 1; ++z) {
      for(int x = x_bin - 1; x <= x_bin + 1; ++x) {
        for(int y = y_bin - 1; y <= y_bin + 1; ++y) {
          grid[z*num_xbins*num_ybins + x*num_ybins + y] = 1;
        }
      }
    }
  }

  // lets call flood fill for every layer along the z-axis from every 
  // corner in the x-y plane. 
  for(int z = 0; z < num_zbins; ++z) {
    char* level = &grid[z*num_xbins*num_ybins];
    FloodLevel(level, 0,             0,             num_xbins, num_ybins, 0, 2);
    FloodLevel(level, 0,             num_ybins - 1, num_xbins, num_ybins, 0, 2);
    FloodLevel(level, num_xbins - 1, 0,             num_xbins, num_ybins, 0, 2);
    FloodLevel(level, num_xbins - 1, num_ybins - 1, num_xbins, num_ybins, 0, 2);
  }

  // every cube in every layer that has currently value 1 that has a city-block
  // distance below 3 to any cube with value 2 is considered to be in contact
  // with the outer surface... lets set them to a value of three
  for(int z = 0; z < num_zbins; ++z) {
    char* level = &grid[z*num_xbins*num_ybins];
    for(int x = 0; x < num_xbins; ++x) {
      for(int y = 0; y < num_ybins; ++y) {
        if(level[x*num_ybins + y] == 1) {
          int x_from = std::max(0, x - 3);
          int x_to = std::min(num_xbins-1, x + 3);
          int y_from = std::max(0, y - 3);
          int y_to = std::min(num_ybins-1, y + 3);
          bool exposed = false;
          for(int i = x_from; i <= x_to && !exposed; ++i) {
            for(int j = y_from; j <= y_to; ++j) {
              if(level[i*num_ybins + j] == 2) {
                level[x*num_ybins + y] = 3;
                exposed = true;
                break;
              }
            }
          }
        }
      }
    }
  }

  // all positions that lie in a cube with value 3 are considered to be exposed...
  exposed_atom_positions.clear();
  exposed_transfer_energies.clear();
  exposed_asas.clear();
  for(uint i = 0; i < atom_positions.size(); ++i) {
    const geom::Vec3& pos = atom_positions[i];
    int x_bin = (pos[0] - min_x) * one_over_radius;
    int y_bin = (pos[1] - min_y) * one_over_radius;
    int z_bin = (pos[2] - min_z) * one_over_radius;
    if(grid[z_bin*num_xbins*num_ybins + x_bin*num_ybins + y_bin] == 3) {
      exposed_atom_positions.push_back(pos);
      exposed_transfer_energies.push_back(transfer_energies[i]);
      exposed_asas.push_back(asas[i]);
    }
  }

  // cleanup
  delete [] grid;
}


void ScanAxis(const std::vector<geom::Vec3>& atom_positions, 
              const std::vector<Real>& transfer_energies, 
              const geom::Vec3& axis,
              Real& best_width, Real& best_center, Real& best_energy) {

  int n_pos = atom_positions.size();

  geom::Vec3 normalized_axis = geom::Normalize(axis);
  std::vector<Real> pos_on_axis(n_pos);

  for(int i = 0; i < n_pos; ++i) {
    pos_on_axis[i] = geom::Dot(atom_positions[i], normalized_axis);
  }

  Real min_pos = pos_on_axis[0];
  Real max_pos = min_pos;

  for(int i = 1; i < n_pos; ++i) {
    min_pos = std::min(min_pos, pos_on_axis[i]);
    max_pos = std::max(max_pos, pos_on_axis[i]);
  }

  min_pos = std::floor(min_pos);
  max_pos = std::ceil(max_pos);

  int full_width = int(max_pos - min_pos) + 1;

  //energies representing the energy profile along the axis
  std::vector<Real> mapped_energies(full_width, 0.0);

  for(int i = 0; i < n_pos; ++i) {
    mapped_energies[int(pos_on_axis[i] - min_pos)] += transfer_energies[i];
  }

  best_width = 0;
  best_center = 0;
  best_energy = std::numeric_limits<Real>::max();

  for(int window_width = 10; window_width <= 40; ++window_width) {
    
    if(window_width > full_width) {
      break;
    }
    
    Real energy=0.0;
    for(int i = 0; i < window_width; ++i) {
      energy += mapped_energies[i];
    }

    if(energy < best_energy) {
      best_width = window_width;
      best_center = min_pos + Real(0.5) * window_width + Real(0.5);
      best_energy = energy;
    }
    
    for(int pos = 1; pos < full_width - window_width; ++pos) {
      energy -= mapped_energies[pos-1];
      energy += mapped_energies[pos+window_width-1];
      if(energy < best_energy){
        best_width = window_width;
        best_center = min_pos + pos + Real(0.5) * window_width + Real(0.5);
        best_energy = energy;
      }
    }
  }
}


struct LMInput {
  ost::mol::alg::FindMemParam mem_param;
  geom::Transform initial_transform;
  std::vector<geom::Vec3> exposed_atom_positions;
  std::vector<Real> exposed_transfer_energies;
  std::vector<Real> exposed_asas;
};


void SampleZ(const std::vector<geom::Vec3>& atom_pos, 
             const std::vector<Real>& transfer_energies,
             const std::vector<Real>& asas,
             const geom::Transform& initial_transform,          
             int n_solutions, std::list<LMInput>& top_solutions) {


  std::vector<geom::Vec3> transformed_atom_pos(atom_pos.size());
  for(uint at_idx = 0; at_idx < atom_pos.size(); ++at_idx) {
    transformed_atom_pos[at_idx] = initial_transform.Apply(atom_pos[at_idx]);
  }

  std::vector<geom::Vec3> exposed_atom_positions;
  std::vector<Real> exposed_transfer_energies;
  std::vector<Real> exposed_asas;
  GetExposedAtoms(transformed_atom_pos, transfer_energies, asas,
                  exposed_atom_positions, exposed_transfer_energies,
                  exposed_asas);

  std::vector<Real> tilt_angles;
  std::vector<Real> rotation_angles;
  for(int tilt_deg = 0; tilt_deg <= 45; tilt_deg += 5) {
    if(tilt_deg == 0) {
      tilt_angles.push_back(0.0);
      rotation_angles.push_back(0.0);
    }
    else {
      Real tilt_angle = Real(tilt_deg) / Real(180.) * Real(M_PI);
      for(int angle_deg = 0; angle_deg < 360; angle_deg += 5) {
        tilt_angles.push_back(tilt_angle);
        rotation_angles.push_back(Real(angle_deg) / Real(180.) * Real(M_PI));
      }
    }
  }

  geom::Vec3 normalized_axis(0.0,0.0,1.0);
  geom::Vec3 tilt_axis(1.0,0.0,0.0);

  for(uint i = 0; i < tilt_angles.size(); ++i) {

    Real tilt_angle = tilt_angles[i];
    Real rotation_angle = rotation_angles[i];

    geom::Vec3 tilted_axis = RotateAroundAxis(normalized_axis, tilt_axis, 
                                              tilt_angle);
    geom::Vec3 scan_axis = RotateAroundAxis(tilted_axis, normalized_axis, 
                                            rotation_angle);

    Real actual_width, actual_center, actual_energy;
    ScanAxis(exposed_atom_positions, exposed_transfer_energies, scan_axis,
             actual_width, actual_center, actual_energy);

    if(static_cast<int>(top_solutions.size()) >= n_solutions && 
       actual_energy > top_solutions.back().mem_param.energy) {
      continue;
    }

    LMInput lm_input;
    lm_input.mem_param.axis = normalized_axis;
    lm_input.mem_param.tilt_axis = tilt_axis;
    lm_input.mem_param.tilt = tilt_angle;
    lm_input.mem_param.angle = rotation_angle;
    lm_input.mem_param.width = actual_width;
    lm_input.mem_param.pos = actual_center;
    lm_input.mem_param.energy = actual_energy;
    lm_input.initial_transform = initial_transform;
    lm_input.exposed_atom_positions = exposed_atom_positions;
    lm_input.exposed_transfer_energies = exposed_transfer_energies;
    lm_input.exposed_asas = exposed_asas;

    if(top_solutions.empty()) {
      top_solutions.push_back(lm_input);
    } 
    else {
      bool added = false;
      for(std::list<LMInput>::iterator sol_it = top_solutions.begin(); 
          sol_it != top_solutions.end(); ++sol_it) {
        if(sol_it->mem_param.energy > lm_input.mem_param.energy) {
          top_solutions.insert(sol_it, lm_input);
          added = true;
          break;
        }
      }

      if(!added) {
        top_solutions.push_back(lm_input);
      }

      while(static_cast<int>(top_solutions.size()) > n_solutions) {
        top_solutions.pop_back();
      }
    }
  }
}


ost::mol::alg::FindMemParam GetFinalSolution(const std::list<LMInput>& top_solutions, 
                                             Real lambda) {

  Real best_energy = std::numeric_limits<Real>::max();
  std::list<LMInput>::const_iterator best_sol_it = top_solutions.begin();
  Eigen::Matrix<Real, 4, 1> lm_parameters;
  Eigen::Matrix<Real, 4, 1> best_lm_parameters;
  LevenbergMarquardt<EnergyF, EnergyDF>::Results lm_result; 

  for(std::list<LMInput>::const_iterator sol_it = top_solutions.begin();
      sol_it != top_solutions.end(); ++sol_it) {

    Real offset = std::max(Real(20000.), std::abs(sol_it->mem_param.energy * 2));

    EnergyF en_f(sol_it->exposed_atom_positions, 
                 sol_it->exposed_transfer_energies, 
                 lambda, offset,
                 sol_it->mem_param.axis, 
                 sol_it->mem_param.tilt_axis);

    lm_parameters(0,0) = sol_it->mem_param.tilt;
    lm_parameters(1,0) = sol_it->mem_param.angle;
    lm_parameters(2,0) = sol_it->mem_param.width;
    lm_parameters(3,0) = sol_it->mem_param.pos;

    LevenbergMarquardt<EnergyF,EnergyDF> lm(en_f);
    lm_result = lm.minimize(&lm_parameters);
 
    Real minimized_energy = en_f(lm_parameters)(0, 0) - en_f.offset;
    
    if(minimized_energy < best_energy) {
      best_energy = minimized_energy;
      best_sol_it = sol_it;
      best_lm_parameters = lm_parameters;
    }
  }

  ost::mol::alg::FindMemParam mem_param = best_sol_it->mem_param;
  mem_param.energy = best_energy;
  mem_param.tilt = best_lm_parameters(0,0);
  mem_param.angle = best_lm_parameters(1,0);
  mem_param.width = best_lm_parameters(2,0);
  mem_param.pos = best_lm_parameters(3,0);

  // assign the membrane accessible surface
  // misuses the optimizer energy function and feeds in the 
  // exposed asa instead of the transfer energies
  EnergyF en_f(best_sol_it->exposed_atom_positions, 
               best_sol_it->exposed_asas, 
               lambda, 0.0,
               best_sol_it->mem_param.axis, 
               best_sol_it->mem_param.tilt_axis);
  mem_param.membrane_asa = en_f(best_lm_parameters)(0, 0);

  // the solution is still relative to the initial transform that has
  // been applied when calling the SampleZ funtion!
  geom::Transform t = best_sol_it->initial_transform;
  mem_param.tilt_axis = t.ApplyInverse(mem_param.tilt_axis);
  mem_param.axis = t.ApplyInverse(mem_param.axis);

  return mem_param;
}


ost::mol::EntityHandle CreateMembraneRepresentation(
                                  const std::vector<geom::Vec3>& atom_positions,
                                  const ost::mol::alg::FindMemParam& param,
                                  Real membrane_margin = 15,
                                  Real delta = 2.0) {

  // let's first construct two planes defining the membrane
  geom::Vec3 membrane_axis = param.GetMembraneAxis();
  geom::Vec3 one = param.pos * membrane_axis + 
                   membrane_axis * param.width / 2;
  geom::Vec3 two = param.pos * membrane_axis - 
                   membrane_axis * param.width / 2;
  geom::Plane plane_one = geom::Plane(one, membrane_axis);
  geom::Plane plane_two = geom::Plane(two, membrane_axis);

  // let's find all positions that are somehow close to those planes
  geom::Vec3List close_pos;
  geom::Vec3List close_pos_one;
  geom::Vec3List close_pos_two;

  for(uint i = 0; i < atom_positions.size(); ++i) {

    Real d1 = geom::Distance(plane_one, atom_positions[i]);
    Real d2 = geom::Distance(plane_two, atom_positions[i]);

    if(d1 < Real(3.)) {
      close_pos_one.push_back(atom_positions[i]);
    }

    if(d2 < Real(3.)) {
      close_pos_two.push_back(atom_positions[i]);
    }

    if(d1 < Real(3.) || d2 < Real(3.)) {
      close_pos.push_back(atom_positions[i]);
    }
  }

  // the geometric center of the close pos vector in combination with the 
  // membrane axis define the central "line" of the disks that will represent 
  // the membrane
  geom::Vec3 center_pos = close_pos.GetCenter();
  geom::Line3 center_line = geom::Line3(center_pos, center_pos + membrane_axis);

  // the final radius of the "disks" is based on the maximal distance of any 
  // position in close_pos to the center_line plus the membrane_margin

  Real max_d_to_center_line = 0;
  for(uint i = 0; i < close_pos.size(); ++i) {
    Real d = geom::Distance(center_line, close_pos[i]);
    max_d_to_center_line = std::max(max_d_to_center_line, d);
  }

  Real disk_radius = max_d_to_center_line + membrane_margin;
  int num_sampling_points = (Real(2.) * disk_radius) / delta;

  // reassign the top and bottom positions, that have been only arbitrary 
  // points on the membrane planes
  one = geom::IntersectionPoint(center_line, plane_one);
  two = geom::IntersectionPoint(center_line, plane_two);

  // find a pair of perpendicular vectors, that are on the plane
  geom::Vec3 arbitrary_vec(1.0, 0.0, 0.0);
  if(geom::Angle(membrane_axis, arbitrary_vec) < 0.1) {
    // parallel is not cool in this case
    arbitrary_vec = geom::Vec3(0.0, 1.0, 0.0);
  }
  geom::Vec3 plane_x = geom::Normalize(geom::Cross(membrane_axis, arbitrary_vec));
  geom::Vec3 plane_y = geom::Normalize(geom::Cross(membrane_axis, plane_x));

  // final representing positions come in here
  std::vector<geom::Vec3> final_pos;

  // do plane one
  geom::Vec3 origin = one - delta * num_sampling_points * 0.5 * plane_x - 
                      delta * num_sampling_points * 0.5 * plane_y;

  for(int i = 0; i < num_sampling_points; ++i) {
    for(int j = 0; j < num_sampling_points; ++j) {
      geom::Vec3 pos = origin + i*delta*plane_x + j*delta*plane_y;
      if(geom::Distance(pos, one) < disk_radius) {
        bool far_far_away = true;
        // this is slow...
        for(uint k = 0; k < close_pos_one.size(); ++k) {
          if(geom::Length2(pos - close_pos_one[k]) < Real(16.)) {
            far_far_away = false;
            break;
          }
        }
        if(far_far_away) {
          final_pos.push_back(pos);
        }
      }
    }
  }

  // do plane two
  origin = two - delta * num_sampling_points * 0.5 * plane_x - 
           delta * num_sampling_points * 0.5 * plane_y;

  for(int i = 0; i < num_sampling_points; ++i) {
    for(int j = 0; j < num_sampling_points; ++j) {
      geom::Vec3 pos = origin + i*delta*plane_x + j*delta*plane_y;
      if(geom::Distance(pos, two) < disk_radius) {
        bool far_far_away = true;
        // this is slow...
        for(uint k = 0; k < close_pos_two.size(); ++k) {
          if(geom::Length2(pos - close_pos_two[k]) < Real(16.)) {
            far_far_away = false;
            break;
          }
        }
        if(far_far_away) {
          final_pos.push_back(pos);
        }
      }
    }
  }

  // create hacky entity that contains membrane representing positions and 
  // return
  ost::mol::EntityHandle membrane_ent = ost::mol::CreateEntity();
  ost::mol::XCSEditor ed = membrane_ent.EditXCS();

  ost::mol::ChainHandle chain = ed.InsertChain("M");
  ost::mol::ResidueHandle res = ed.AppendResidue(chain, "MEM");
  String atom_names = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  uint atom_name_idx = 0;
  uint atom_name_secondary_idx = 0;

  for(uint i = 0; i < final_pos.size(); ++i) {

    if(atom_name_secondary_idx == atom_names.size()) {
      ++atom_name_idx;
      atom_name_secondary_idx = 0;
    }
    if(atom_name_idx == atom_names.size()) {
      res = ed.AppendResidue(chain, "MEM");
      atom_name_idx = 0;
      atom_name_secondary_idx = 0;
    }

    String atom_name = "--";
    atom_name[0] = atom_names[atom_name_idx];
    atom_name[1] = atom_names[atom_name_secondary_idx];

    ed.InsertAtom(res, atom_name, final_pos[i]);
    ++atom_name_secondary_idx;
  }

  return membrane_ent;
}


} // anon namespace


namespace ost{ namespace mol{ namespace alg{

geom::Vec3 FindMemParam::GetMembraneAxis() const {  

  geom::Vec3 result = RotateAroundAxis(axis,tilt_axis,tilt);
  result = RotateAroundAxis(result,axis,angle);
  return result;
}


FindMemParam FindMembrane(ost::mol::EntityView& ent, 
                          bool assign_membrane_representation, bool fast) {

  ost::mol::EntityView peptide_view = ent.Select("peptide=true and ele!=H");
  Accessibility(peptide_view);

  std::vector<geom::Vec3> atom_pos;
  std::vector<Real> transfer_energies;
  std::vector<Real> asas;

  atom_pos.reserve(peptide_view.GetAtomCount());
  transfer_energies.reserve(peptide_view.GetAtomCount());

  ost::mol::AtomViewList atoms = peptide_view.GetAtomList();
  String stupid_string("S_N_O_C");

  for(ost::mol::AtomViewList::iterator it = atoms.begin(); 
      it != atoms.end(); ++it) {

    if(!it->HasProp("asaAtom")) {
      continue;
    }

    String element = it->GetElement();
    if(stupid_string.find(element) == std::string::npos) {
      continue;
    }

    Real asa = it->GetFloatProp("asaAtom");
    asas.push_back(asa);
    atom_pos.push_back(it->GetPos());

    if(element == "S") {
      transfer_energies.push_back(asa * Real(10.0));
    }
    else if(element == "N") {
      transfer_energies.push_back(asa * Real(53.0));
    }
    else if(element == "O") {
      transfer_energies.push_back(asa * Real(57.0));
    }
    else if(element == "C") {
      // check whether we find a double bond to distinguish between 
      // hibridization states
      bool assigned_energy = false;
      ost::mol::BondHandleList bond_list = it->GetBondList();
      for(ost::mol::BondHandleList::iterator bond_it = bond_list.begin();
          bond_it != bond_list.end(); ++bond_it){
        unsigned char bond_order = bond_it->GetBondOrder();
        if(bond_order > '1'){
          transfer_energies.push_back(asa * Real(-19.0));
          assigned_energy = true;
          break;;
        }
      }
      if(!assigned_energy) {
        transfer_energies.push_back(asa * Real(-22.6));
      }
    }
  }

  if(atom_pos.size() < 10) {
    throw ost::Error("Cannot detect membrane with such a low number "
                     "of heavy atoms!");
  }

  // we always optimizer along the z-axis. 
  // We therefore have to transform the positions. We use a rotation 
  // around the z-axis with subsequent rotation around the x-axis for this task
  std::vector<geom::Transform> transformations;
  int n_euler_angles = 3;
  int n_transformations = n_euler_angles * n_euler_angles * n_euler_angles;
  std::vector<Real> euler_angles(n_euler_angles);
  euler_angles[0] = 0.0;
  euler_angles[1] = M_PI/3;
  euler_angles[2] = 2*M_PI/3;

  for(int i = 0; i < n_euler_angles; ++i) {
    for(int j = 0; j < n_euler_angles; ++j) {
      for(int k = 0; k < n_euler_angles; ++k) {
        geom::Mat3 rot_matrix = geom::EulerTransformation(euler_angles[i],
                                                          euler_angles[j],
                                                          euler_angles[k]);
        geom::Transform transform;
        transform.SetRot(rot_matrix);
        transformations.push_back(transform);
      }
    }
  }

  // lets use the generated transforms to search for initial solutions that can
  // then be fed into a final minimization...
  std::list<LMInput> top_solutions;
  int n_initial_solutions = fast ? 1 : 20;

  for(int transformation_idx = 0; transformation_idx < n_transformations;
      ++transformation_idx) {
    SampleZ(atom_pos, transfer_energies, asas, transformations[transformation_idx],
            n_initial_solutions, top_solutions);
  }

  // Perform the final minimization and return the best solution.
  // please note, that the returned solution is transformed back in order
  // to match the initial atom positions
  FindMemParam final_solution = GetFinalSolution(top_solutions, 0.9);

  if(assign_membrane_representation) {
    final_solution.membrane_representation = CreateMembraneRepresentation(
                                                                atom_pos,
                                                                final_solution);
  }

  return final_solution;
}


FindMemParam FindMembrane(ost::mol::EntityHandle& ent, 
                          bool assign_membrane_representation,
                          bool fast) {

  ost::mol::EntityView ent_view = ent.CreateFullView();
  return FindMembrane(ent_view, assign_membrane_representation, fast);
}

}}} // ns
