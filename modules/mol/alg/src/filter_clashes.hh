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
#ifndef OST_MOL_ALG_FILTER_CLASHES_HH
#define OST_MOL_ALG_FILTER_CLASHES_HH

#include <ost/mol/entity_view.hh>
#include <ost/mol/alg/module_config.hh>
#include <ost/mol/alg/distance_test_common.hh>
namespace ost { namespace mol { namespace alg {

class DLLEXPORT_OST_MOL_ALG BondLengthInfo
{
public:
  BondLengthInfo(): avg_length_(0),avg_zscore_(0),count_(0) {}
  BondLengthInfo(Real avg_length,Real avg_zscore, int count):
    avg_length_(avg_length),
    avg_zscore_(avg_zscore),
    count_(count){}
  Real GetAvgLength() {return avg_length_;}
  Real GetAvgZscore() {return avg_zscore_;}
  int GetCount() {return count_;}

private:
  Real avg_length_;
  Real avg_zscore_;
  int count_;
};

class DLLEXPORT_OST_MOL_ALG ClashEvent
{
public:
  ClashEvent(): atom1_(UniqueAtomIdentifier()),atom2_(UniqueAtomIdentifier()),mdl_dist_(0.0),adjusted_ref_dist_(0.0) {}
  ClashEvent(const UniqueAtomIdentifier& atom1,const UniqueAtomIdentifier& atom2, Real mdl_dist, Real adjusted_ref_dist ):
    atom1_(atom1),atom2_(atom2),mdl_dist_(mdl_dist),adjusted_ref_dist_(adjusted_ref_dist) {}
  UniqueAtomIdentifier GetFirstAtom() const { return atom1_; }
  UniqueAtomIdentifier GetSecondAtom() const { return atom2_; }
  Real GetModelDistance() const { return mdl_dist_; }
  Real GetAdjustedReferenceDistance() const { return adjusted_ref_dist_; }
private:
  UniqueAtomIdentifier atom1_;
  UniqueAtomIdentifier atom2_;
  Real mdl_dist_;
  Real adjusted_ref_dist_;
};

class DLLEXPORT_OST_MOL_ALG ClashingInfo
{

public:
  ClashingInfo(): clash_count_(0), average_offset_ (0), clash_list_(std::vector<ClashEvent>()) {}
  ClashingInfo (int clash_count, Real average_offset, const std::vector<ClashEvent>& clash_list):
    clash_count_(clash_count), average_offset_ (average_offset),clash_list_(clash_list) {}
  int GetClashCount() const {return clash_count_/2.0;}
  Real GetAverageOffset() const {return average_offset_;}
  std::vector<ClashEvent> GetClashList() const;

private:
  int clash_count_;
  Real average_offset_;
  std::vector<ClashEvent> clash_list_;
};

class DLLEXPORT_OST_MOL_ALG StereoChemicalBondViolation
{
public:
  StereoChemicalBondViolation():
    atom1_(UniqueAtomIdentifier()),atom2_(UniqueAtomIdentifier()),mdl_value_(0.0),allowed_range_(std::pair<Real,Real>(0.0,0.0)) {}
  StereoChemicalBondViolation(const UniqueAtomIdentifier& atom1,
                              const UniqueAtomIdentifier& atom2,
                              Real mdl_value, std::pair<Real,Real> allowed_range ):
    atom1_(atom1),atom2_(atom2),mdl_value_(mdl_value),allowed_range_(allowed_range) {}
  UniqueAtomIdentifier GetFirstAtom() const { return atom1_; }
  UniqueAtomIdentifier GetSecondAtom() const { return atom2_; }
  Real GetBondLength() const { return mdl_value_; }
  std::pair<Real,Real> GetAllowedRange() const { return allowed_range_; }
private:
  UniqueAtomIdentifier atom1_;
  UniqueAtomIdentifier atom2_;
  Real mdl_value_;
  std::pair<Real,Real> allowed_range_;
};

class DLLEXPORT_OST_MOL_ALG StereoChemicalAngleViolation
{
public:
  StereoChemicalAngleViolation():
    atom1_(UniqueAtomIdentifier()),atom2_(UniqueAtomIdentifier()),atom3_(UniqueAtomIdentifier()),mdl_value_(0.0),allowed_range_(std::pair<Real,Real>(0.0,0.0)) {}
  StereoChemicalAngleViolation(const UniqueAtomIdentifier& atom1,
                               const UniqueAtomIdentifier& atom2,
                               const UniqueAtomIdentifier& atom3,
                               Real mdl_value, std::pair<Real,Real> allowed_range ):
    atom1_(atom1),atom2_(atom2),atom3_(atom3),mdl_value_(mdl_value),allowed_range_(allowed_range) {}
  UniqueAtomIdentifier GetFirstAtom() const { return atom1_; }
  UniqueAtomIdentifier GetSecondAtom() const { return atom2_; }
  UniqueAtomIdentifier GetThirdAtom() const { return atom3_; }
  Real GetAngleWidth() const { return mdl_value_; }
  std::pair<Real,Real> GetAllowedRange() const { return allowed_range_; }
private:
  UniqueAtomIdentifier atom1_;
  UniqueAtomIdentifier atom2_;
  UniqueAtomIdentifier atom3_;
  Real mdl_value_;
  std::pair<Real,Real> allowed_range_;
};

class DLLEXPORT_OST_MOL_ALG StereoChemistryInfo
{
public:
  StereoChemistryInfo():
     avg_zscore_bonds_(0),
     bad_bond_count_(0),
     bond_count_(0),
     avg_zscore_angles_(0),
     bad_angle_count_(0),
     angle_count_(0),
     avg_bond_length_info_(std::map<String,BondLengthInfo>()),
     bond_violation_list_(std::vector<StereoChemicalBondViolation>()),
     angle_violation_list_(std::vector<StereoChemicalAngleViolation>()) {}
  StereoChemistryInfo(Real avg_zscore_bonds, int bad_bond_count, int bond_count,
                      Real avg_zscore_angles, int bad_angle_count, int angle_count,
                      const std::map<String,BondLengthInfo>& avg_bond_length_info,
                      const std::vector<StereoChemicalBondViolation>& bond_violation_list,
                      const std::vector<StereoChemicalAngleViolation>& angle_violation_list):
    avg_zscore_bonds_(avg_zscore_bonds),
    bad_bond_count_(bad_bond_count),
    bond_count_(bond_count),
    avg_zscore_angles_(avg_zscore_angles),
    bad_angle_count_(bad_angle_count),
    angle_count_(angle_count),
    avg_bond_length_info_(avg_bond_length_info),
    bond_violation_list_(bond_violation_list),
    angle_violation_list_(angle_violation_list) {}
  Real GetAvgZscoreBonds() const {return avg_zscore_bonds_;}
  int GetBadBondCount() const {return bad_bond_count_;}
  int GetBondCount() const {return bond_count_;}
  Real GetAvgZscoreAngles() const {return avg_zscore_angles_;}
  int GetBadAngleCount() const {return bad_angle_count_;}
  int GetAngleCount() const {return angle_count_;}
  std::map<String,BondLengthInfo> GetAvgBondLengthInfo() {return avg_bond_length_info_;}
  std::vector<StereoChemicalBondViolation> GetBondViolationList() { return bond_violation_list_; }
  std::vector<StereoChemicalAngleViolation> GetAngleViolationList() { return angle_violation_list_; }


private:
  Real avg_zscore_bonds_;
  int bad_bond_count_;
  int bond_count_;
  Real avg_zscore_angles_;
  int bad_angle_count_;
  int angle_count_;
  std::map<String,BondLengthInfo> avg_bond_length_info_;
  std::vector<StereoChemicalBondViolation> bond_violation_list_;
  std::vector<StereoChemicalAngleViolation> angle_violation_list_;
};


/// \brief List of reference atom-atom distances to detect clashes between non-bonded atoms 
class DLLEXPORT_OST_MOL_ALG ClashingDistances
{

public:
  /// \brief Default constructor (creates an empty list)
  ClashingDistances(): valid_flag_(true) {}

  /// \brief Adds or replaces an entry 
  ///
  /// Requires an atom-atom distance and a tolerance threshold
  void SetClashingDistance(const String& ele1,const String& ele2, Real min_distance, Real tolerance);

  /// \brief Recovers a reference distance and a tolerance threshold (respectively) from the list
  std::pair<Real,Real> GetClashingDistance(const String& ele1,const String& ele2) const;

  /// \brief Recovers a reference distance already adjusted by the tolerance threshold from the list
  Real GetAdjustedClashingDistance(const String& ele1,const String& ele2) const;

  /// \brief Recovers the longest distance in the list, corrected by tolerance
  Real GetMaxAdjustedDistance() const;

  /// \brief Returns true if the list is empty (i.e. in an invalid, useless state)
  bool IsEmpty() const;
  
  /// \brief Prints all distances in the list to standard output
  void PrintAllDistances() const;
  
private:

  std::map <String,std::pair<Real,Real> > min_distance_;
  Real default_min_distance_;
  Real default_min_distance_tolerance_;
  bool valid_flag_;
  
};

/// \brief List of stereo chemical parameters (Bonds and angles)
///
/// For each item (bond or angle in a specific residue), stores the mean and standard deviation 
class DLLEXPORT_OST_MOL_ALG StereoChemicalParams
{

public:
  /// \brief Adds or replaces an entry 
  void SetParam(const String& param, const String& residue, Real value, Real st_dev);

  /// \brief Recovers mean and standard deviation (respectively) of a stereo⁻chemical item (bond or angle) from the list
  ///
  /// Item format: Bond: X-Y, Angle:X-Y-Z
  std::pair<Real,Real> GetParam(const String& element,const String& residue) const;

  /// \brief Checks if the list contains an entry for a specific stereo-chemical item (a bond or atom in a specific residue)
  ///
  /// Item format: Bond: X-Y, Angle:X-Y-Z
  bool ContainsParam(const String& param,const String& residue) const;

  /// \brief Returns true if the list is empty (i.e. in an invalid, useless state)
  ///
  /// Item format: Bond: X-Y, Angle:X-Y-Z
  bool IsEmpty() const;
  
  /// \brief Prints all distances in the list to standard output
  void PrintAllParameters() const;
  
private:

  std::map<std::pair<String,String>,std::pair<Real,Real> >  params_;
 
};

/// \brief Fills a list of reference clashing distances from the content of a parameter file
///
/// Requires a list of strings holding the contents of a parameter file, one line per string
ClashingDistances DLLEXPORT_OST_MOL_ALG FillClashingDistances(std::vector<String>& stereo_chemical_props_file, bool check=false);

/// \brief Fills a list of stereo-chemical statistics from the content of a parameter file
///
/// Requires a list of strings holding the contents of a parameter file, one line per string
/// The header can be 'Bonds' to read bond statistics or 'Angles' to read angle statistics
StereoChemicalParams DLLEXPORT_OST_MOL_ALG FillStereoChemicalParams(const String& header, std::vector<String>& stereo_chemical_props_file, bool check=false);  

/// \brief Filters a structure based on detected clashes between non bonded atoms. Entity version
///
/// If a clash between two atoms (distance shorter than reference clashing distance - tolerance threshold)
/// is detected in a residue's side-chain, all atoms in the side chain are removed from the structure
/// If a clash is detected in the backbone, all atoms in the residue are removed. This behavior is changed 
/// by the always_remove_bb flag: when the flag is set to true all atoms in the residue are removed even if
/// a clash is just detected in the side-chain
std::pair<EntityView,ClashingInfo> DLLEXPORT_OST_MOL_ALG FilterClashes(const EntityView& ent,
                                               const ClashingDistances& min_distances, bool always_remove_bb=false);

/// \brief Filters a structure based on detected clashes between non bonded atoms. Handle version
///
/// If a clash between two atoms (distance shorter than reference clashing distance - tolerance threshold)
/// is detected in a residue's side-chain, all atoms in the side chain are removed from the structure
/// If a clash is detected in the backbone, all atoms in the residue are removed. This behavior is changed 
/// by the always_remove_bb flag: when the flag is set to true all atoms in the residue are removed even if
/// a clash is just detected in the side-chain
std::pair<EntityView,ClashingInfo> DLLEXPORT_OST_MOL_ALG FilterClashes(const EntityHandle& ent,
                                               const ClashingDistances& min_distances, bool always_remove_bb=false);

/// \brief Filters a structure based on detected stereo-chemical violations. Entity version
///
/// If a stereo-chemical violation (i.e., a bond or an angle with a value outside the range defined by
/// the mean value, the standard deviation and the tolerance parameter) is detected in a residue's side-chain, 
/// all atoms in the side chain are removed from the structure. If a violation is detected in the backbone, all 
/// atoms in the residue are removed. This behavior is changed by the always_remove_bb flag: when the flag is 
/// set to true all atoms in the residue are removed even if a violation is just detected in the side-chain
std::pair<EntityView,StereoChemistryInfo> DLLEXPORT_OST_MOL_ALG CheckStereoChemistry(const EntityView& ent,
                                                      const StereoChemicalParams& bond_table, 
                                                      const StereoChemicalParams& angle_table,
                                                      Real bond_tolerance,
                                                      Real angle_tolerance,
                                                      bool always_remove_bb=false);

/// \brief Filters a structure based on detected stereo-chemical violations. Handle version
///
/// If a stereo-chemical violation (i.e., a bond or an angle with a value outside the range defined by
/// the mean value, the standard deviation and the tolerance parameter) is detected in a residue's side-chain, 
/// all atoms in the side chain are removed from the structure. If a violation is detected in the backbone, all 
/// atoms in the residue are removed. This behavior is changed by the always_remove_bb flag: when the flag is 
/// set to true all atoms in the residue are removed even if a violation is just detected in the side-chain
std::pair<EntityView,StereoChemistryInfo> DLLEXPORT_OST_MOL_ALG CheckStereoChemistry(const EntityHandle& ent,
                                                      const StereoChemicalParams& bond_table, 
                                                      const StereoChemicalParams& angle_table,
                                                      Real bond_tolerance,
                                                      Real angle_tolerance,
                                                      bool always_remove_bb=false);


}}}


#endif
