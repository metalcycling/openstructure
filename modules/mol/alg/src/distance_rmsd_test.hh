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
#ifndef OST_MOL_ALG_DISTANCE_RMSD_TEST_HH
#define OST_MOL_ALG_DISTANCE_RMSD_TEST_HH

#include <ost/mol/entity_view.hh>
#include <ost/mol/alg/module_config.hh>
#include <ost/mol/alg/local_dist_diff_test.hh>

namespace ost { namespace mol { namespace alg {
  
/// \brief Calculates the Distance RMSD in a model, given a list of distances with their refence length
///
/// Calculates the two values needed to determine the Distance RMSD Test score for a given model, i.e.
/// the sum of the squared distance deviations and the total number of distances on which the sum was
/// calculated. 
/// The function requires a list of distances for which the distance deviation has to be computed, together 
/// with a reference length for each of them. The function also requires a model on which the returned
/// values are computed. Furthermore, the function requires a default deviation value to be used in the 
/// calculations when a distance is not present in the model.
///
/// The distance information needs to be stored in ian instance of the GlobalRDMap object. This object
/// stores two distance lengths for each distance. This function uses the first of them as reference length 
/// and ignores the second. 
/// 
/// The function only processes standard residues in the first chains of the model and of the reference
/// For residues with symmetric sidechains (GLU, ASP, ARG, VAL, PHE, TYR), the 
/// naming of the atoms is ambigous. For these residues, the function computes the Distance RMSD 
/// score that each naming convention would generate when considering all non-ambigous surrounding atoms.
/// The solution that gives lower score is then picked to compute the final Distance RMSD score for the whole 
/// model.
///
/// A sequence separation parameter can be passed to the function. If this happens, only distances between residues
/// whose separation is higher than the provided parameter are considered when computing the score. 
///
/// If a string is provided as an argument to the function, residue-per-residue statistics are stored as 
/// residue properties. Specifically, the local residue-based Distance RMSD score is stored in a float property 
/// as the provided string, while the residue-based sum of squared distances and the number of distances checked 
/// are saved in two properties named [string]_sum (a float property) and [string]_count (an int property).
std::pair<Real,long int> DLLEXPORT_OST_MOL_ALG 
DistanceRMSDTest(const EntityView& mdl,
                 const GlobalRDMap& glob_dist_list,
                 Real cap_distance, int sequence_separation = 0,
                 const String& local_drmsdt_property_string="");


/// \brief Computes the Distance RMSD Test given a list of distances to check and a model 
///
/// Computes the Local Distance RMSD Test on the provided model. Requires a list of distances to check and a 
/// model for which the score is computed, plus a default deviation value to be used when a distance is not
/// present in the model.
///
/// A sequence separation parameter can be passed to the function. If this happens, only distances between residues
/// whose separation is higher than the provided parameter are considered when computing the score.
Real DLLEXPORT_OST_MOL_ALG DRMSD(const EntityView& v, 
                                 const GlobalRDMap& global_dist_list, 
                                 Real cap_distance, int sequence_separation=0);

}}}

#endif


