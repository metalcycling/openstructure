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
#ifndef OST_MOL_ALG_LOCAL_DIST_TEST_HH
#define OST_MOL_ALG_LOCAL_DIST_TEST_HH

#include <ost/mol/alg/module_config.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/mol/alg/distance_test_common.hh>
#include <ost/mol/alg/filter_clashes.hh>

namespace ost { namespace mol { namespace alg {

struct StereoChemicalProps
{
  bool is_valid;
  ost::mol::alg::StereoChemicalParams bond_table;
  ost::mol::alg::StereoChemicalParams angle_table;
  ost::mol::alg::ClashingDistances nonbonded_table;

  StereoChemicalProps();
  StereoChemicalProps(ost::mol::alg::StereoChemicalParams& init_bond_table,
                      ost::mol::alg::StereoChemicalParams& init_angle_table,
                      ost::mol::alg::ClashingDistances& init_nonbonded_table);
};

struct lDDTSettings {
  Real radius; 
  int sequence_separation;
  std::vector<Real> cutoffs;
  String label;

  lDDTSettings();
  lDDTSettings(Real init_radius, 
               int init_sequence_separation,
               std::vector<Real>& init_cutoffs,
               String init_label);
  void PrintParameters();
  std::string ToString();
};

struct lDDTLocalScore {
  String cname;
  String rname;
  int rnum;
  String is_assessed;
  String quality_problems;
  Real local_lddt;
  int conserved_dist;
  int total_dist;

  lDDTLocalScore();

  lDDTLocalScore(String init_cname,
                 String init_rname,
                 int init_rnum,
                 String init_is_assessed,
                 String init_quality_problems,
                 Real init_local_lddt,
                 int init_conserved_dist,
                 int init_total_dist);

  String ToString(bool structural_checks) const;
  String Repr() const;

  static String GetHeader(bool structural_checks, int cutoffs_length);
};

class lDDTScorer
{
  public:
    lDDTSettings settings;
    EntityView model_view;
    std::vector<EntityView> references_view;
    GlobalRDMap glob_dist_list;

    lDDTScorer(std::vector<EntityView>& init_references,
               ost::mol::EntityView& init_model,
               lDDTSettings& init_settings);
    Real GetGlobalScore();
    std::vector<lDDTLocalScore> GetLocalScores();
    int GetNumConservedContacts(); // number of conserved distances in the model
    int GetNumTotalContacts(); // the number of total distances in the reference structure
    std::vector<EntityView> GetReferences();
    void PrintPerResidueStats();
    bool IsValid();

  private:
    bool _score_calculated;
    bool _score_valid;
    bool _has_local_scores;
    // number of conserved distances in the model and
    // the number of total distances in the reference structure
    int _num_cons_con;
    int _num_tot_con;
    Real _global_score;
    std::vector<lDDTLocalScore> _local_scores;
    void _ComputelDDT();
    void _GetLocallDDT();
    void _PrepareGlobalRDMap();
};

std::pair<int,int> DLLEXPORT_OST_MOL_ALG ComputeCoverage(const EntityView& v,const GlobalRDMap& glob_dist_list);

bool DLLEXPORT_OST_MOL_ALG IsResnumInGlobalRDMap(const ResNum& resnum, const GlobalRDMap& glob_dist_list);
  
/// \brief Calculates number of distances conserved in a model, given a list of distances to check and a model
///
/// Calculates the two values needed to determine the Local Distance Difference Test for a given model, i.e.
/// the number of conserved distances in the model and the number of total distances in the reference structure. 
/// The function requires a list of distances to check, a model on which the distances are checked, and a 
/// list of tolerance thresholds that are used to determine if the distances are conserved. 
/// 
/// The function only processes standard residues in the first chains of the model and of the reference
/// For residues with symmetric sidechains (GLU, ASP, ARG, VAL, PHE, TYR), the 
/// naming of the atoms is ambigous. For these residues, the function computes the Local Distance Difference 
/// Test score that each naming convention would generate when considering all non-ambigous surrounding atoms.
/// The solution that gives higher score is then picked to compute the final Local Difference
/// Distance Score for the whole model.
///
/// A sequence separation parameter can be passed to the function. If this happens, only distances between residues
/// whose separation is higher than the provided parameter are considered when computing the score. 
///
/// If a string is provided as an argument to the function, residue-per-residue statistics are stored as 
/// residue properties. Specifically, the local residue-based lddt score is stored in a float property named
/// as the provided string, while the residue-based number of conserved and total distances are saved in two 
/// int properties named [string]_conserved and [string]_total.
std::pair<long int,long int> DLLEXPORT_OST_MOL_ALG 
LocalDistDiffTest(const EntityView& mdl, const GlobalRDMap& dist_list,
                  std::vector<Real> cutoff_list, int sequence_separation = 0, 
                  const String& local_ldt_property_string="");

/// \brief Calculates the Local Distance Difference Score for a given model with respect to a given target
///
/// Calculates the Local Distance Difference Test score for a given model with respect to a given reference structure. Requires
/// a model, a reference structure, a list of thresholds that are used to determine if distances are conserved, and an inclusion
/// radius value used to determine which distances are checked.
/// 
/// The function only processes standard residues in the first chains of the model and of the reference
/// For residues with symmetric sidechains (GLU, ASP, ARG, VAL, PHE, TYR), the 
/// naming of the atoms is ambigous. For these residues, the function computes the Local Distance Difference 
/// Test score that each naming convention would generate when considering all non-ambigous surrounding atoms.
/// The solution that gives higher score is then picked to compute the final Local Difference
/// Distance Score for the whole model.
///
/// If a string is provided as an argument to the function, residue-per-residue statistics are stored as 
/// residue properties. Specifically, the local residue-based lddt score is stored in a float property named
/// as the provided string, while the residue-based number of conserved and total distances are saved in two 
/// int properties named [string]_conserved and [string]_total.
Real DLLEXPORT_OST_MOL_ALG LocalDistDiffTest(const EntityView& mdl,
                                         const EntityView& target,
                                         Real cutoff, 
                                         Real max_dist,
                                         const String& local_ldt_property_string="");
/// \brief Wrapper around LocalDistDiffTest
Real DLLEXPORT_OST_MOL_ALG LocalDistDiffTest(const EntityView& v,
                       std::vector<EntityView>& ref_list,
                       const GlobalRDMap& glob_dist_list,
                       lDDTSettings& settings);

/// \brief Calculates the Local Distance Difference Test score for a given model starting from an alignment between a reference structure and the model. 
///
/// Calculates the Local Distance Difference Test score given an alignment between a model and a taget structure.
/// Requires a threshold on which to calculate the score and an inclusion radius to determine the interatiomic 
/// distances to check. Obviously, the strucvtures of the model and the reference must be attached to the alignment.
/// By default the first structure in the alignment is considered the reference and the second is considered the
/// model, but this can be changed by passing to the function the indexes of the two structures in the ref_index
/// and mdl_index parameters.
/// BEWARE: This algorithm uses the old version of the Local Distance Difference Test 
/// (multiple cycles, single threshold, etc. ) and will give a slightly different result than the other functions
Real DLLEXPORT_OST_MOL_ALG LocalDistDiffTest(const ost::seq::AlignmentHandle& aln,
                                         Real cutoff, Real max_dist,
                                         int ref_index=0, int mdl_index=1);

/// \brief Computes the Local Distance Difference High-Accuracy Test given a list of distances to check
///
/// Computes the Local Distance Difference High-Accuracy Test (with threshold 0.5,1,2 and 4 Angstrom)
/// Requires a list of distances to check and a model for which the score is computed
///
/// A sequence separation parameter can be passed to the function. If this happens, only distances between residues
/// whose separation is higher than the provided parameter are considered when computing the score. 
Real DLLEXPORT_OST_MOL_ALG LDDTHA(EntityView& v, const GlobalRDMap& global_dist_list, int sequence_separation=0);

/// \brief Creates a list of distances to check during a Local Difference Distance Test
///
/// Requires a reference structure and an inclusion radius (max_dist)
GlobalRDMap DLLEXPORT_OST_MOL_ALG CreateDistanceList(const EntityView& ref,Real max_dist);

/// \brief Creates a list of distances to check during a Local Difference Distance Test starting from multiple reference structures
///
/// Requires a list of reference structure and an inclusion radius (max_dist).
/// 
/// The structures in the list have to be properly prepared before being passed 
/// to the function. Corresponding residues in the structures must have the same residue number, the same chain name,
/// etc. Gaps are allowed and automatically dealt with: if information about a distance is present in at least one of
/// the structures, it will be considered.  
///
///
/// If a distance between two atoms is shorter than the inclusion radius in all structures in which the two atoms are 
/// present, it is included in the list. However, if the distance is longer than the inclusion radius in at least one 
/// of the structures, it is not be considered a local interaction and is excluded from the list
///
/// The function takes care of residues with ambigous symmetric sidechains. To decide which naming convention to use, the functions
/// computes a local distance score of each reference structure with the first reference structure in the list, using only non ambigously-named atoms. 
/// It picks then the naming convention that gives the highest score, guaranteeing that all references are processed with the correct atom names.
///
/// The cutoff list that will later be used to compute the Local Distance Difference Test score and the  sequence separation parameter 
/// must be passed to the function. These parameters do not influence the output distance list, which always includes all distances
/// within the provided max_dist (to make it consistent with the single-reference corresponding function). However, the parameters are used when
/// dealing with the naming convention of residues with ambiguous nomenclature. 
GlobalRDMap DLLEXPORT_OST_MOL_ALG CreateDistanceListFromMultipleReferences(const std::vector<EntityView>& ref_list,std::vector<Real>& cutoff_list, int sequence_separation, Real max_dist);

/// \brief Prints all distances in a global distance list to standard output
void DLLEXPORT_OST_MOL_ALG PrintGlobalRDMap(const GlobalRDMap& glob_dist_list);

/// \brief Prints all distances in a residue distance list to standard output
void DLLEXPORT_OST_MOL_ALG PrintResidueRDMap(const ResidueRDMap& res_dist_list);

// required by some helper function. Cannot reuse similar functions in other modules without creating
// circular dependencies
bool DLLEXPORT_OST_MOL_ALG IsStandardResidue(String rn);

void DLLEXPORT_OST_MOL_ALG CleanlDDTReferences(std::vector<EntityView>& ref_list);

// Prepare GlobalRDMap from reference list
GlobalRDMap DLLEXPORT_OST_MOL_ALG PreparelDDTGlobalRDMap(
    const std::vector<EntityView>& ref_list,
    std::vector<Real>& cutoff_list,
    int sequence_separation,
    Real max_dist);

void DLLEXPORT_OST_MOL_ALG CheckStructure(EntityView& ent,
                                          StereoChemicalParams& bond_table,
                                          StereoChemicalParams& angle_table,
                                          ClashingDistances& nonbonded_table,
                                          Real bond_tolerance,
                                          Real angle_tolerance);

std::vector<lDDTLocalScore> DLLEXPORT_OST_MOL_ALG GetlDDTPerResidueStats(EntityView& model,
                                              GlobalRDMap& glob_dist_list,
                                              bool structural_checks,
                                              String label);

void DLLEXPORT_OST_MOL_ALG PrintlDDTPerResidueStats(std::vector<lDDTLocalScore>& scores,
                                                    bool structural_checks,
                                                    int cutoffs_length);

}}}

#endif
