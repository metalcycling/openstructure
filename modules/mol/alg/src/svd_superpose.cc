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

#include <stdexcept>
#include <iostream>

#include <Eigen/SVD>
#include <Eigen/Geometry>

#include <ost/base.hh>
#include <ost/geom/vec3.hh>
#include <ost/geom/mat4.hh>
#include <ost/mol/alg/svd_superpose.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/view_op.hh>
#include <ost/mol/atom_view.hh>

namespace ost { namespace mol { namespace alg {


typedef Eigen::Matrix<Real,3,1> EVec3;
typedef Eigen::Matrix<Real,3,3> EMat3;
typedef Eigen::Matrix<Real,4,4> EMat4;
typedef Eigen::Matrix<Real,1,3> ERVec3;
typedef Eigen::Matrix<Real,Eigen::Dynamic,Eigen::Dynamic> EMatX;
typedef Eigen::Matrix<Real,1,Eigen::Dynamic> ERVecX;





Real calc_rmsd_for_point_lists(const std::vector<geom::Vec3>& points1,
                              const std::vector<geom::Vec3>& points2,
                              const geom::Mat4& transformation)
{
  
  Real rmsd=0.0;
  std::vector<geom::Vec3>::const_iterator a_ev1=points1.begin();
  std::vector<geom::Vec3>::const_iterator a_ev1_end=points1.end();
  std::vector<geom::Vec3>::const_iterator a_ev2=points2.begin();
  for (; a_ev1_end!=a_ev1; ++a_ev1, ++a_ev2) {
    Real val=geom::Length2(geom::Vec3(transformation*geom::Vec4(*a_ev1))-*a_ev2);
    rmsd+=val;
  }
  rmsd=sqrt(rmsd/points1.size());
  return rmsd;
}

Real calc_rmsd_for_atom_lists(const mol::AtomViewList& atoms1,
                              const mol::AtomViewList& atoms2,
                              const geom::Mat4& transformation)
{
  Real rmsd=0.0;
  mol::AtomViewList::const_iterator a_ev1=atoms1.begin();
  mol::AtomViewList::const_iterator a_ev1_end=atoms1.end();
  mol::AtomViewList::const_iterator a_ev2=atoms2.begin();
  for (; a_ev1_end!=a_ev1; ++a_ev1, ++a_ev2) {
    mol::AtomView av1=*a_ev1;
    geom::Vec3 apos=av1.GetPos();
    mol::AtomView av2=*a_ev2;
    geom::Vec3 bpos=av2.GetPos();
    Real val=geom::Length2(geom::Vec3(transformation*geom::Vec4(apos))-bpos);
    rmsd+=val;
  }
  rmsd=sqrt(rmsd/atoms1.size());
  return rmsd;
}

Real calc_rmsd_for_ematx(const EMatX& atoms1,
                         const EMatX& atoms2,
                         const EMat4& transformation)
{
  EMatX transformed = EMatX::Zero(atoms1.rows(), 3);

  EMatX vector = EMatX::Zero(4,1);
  EMatX transformed_vector = EMatX::Zero(4,1);
  vector(3,0)=1;

  for(int i=0;i<atoms1.rows();++i){
    vector.block<3,1>(0,0)=atoms1.block<1,3>(i,0).transpose();
    transformed_vector = transformation*vector;
    transformed.block<1,3>(i,0)=transformed_vector.block<3,1>(0,0).transpose();
  }

  transformed = transformed - atoms2;
  transformed = transformed.array()*transformed.array();
  return sqrt(transformed.rowwise().sum().sum()/atoms1.rows());
}

Real CalculateRMSD(const mol::EntityView& ev1,
                   const mol::EntityView& ev2,
                   const geom::Mat4& transformation) {

  if(ev1.GetAtomCount() != ev2.GetAtomCount()) {
    std::stringstream ss;
    ss << "Number of atoms in input views must be consistent! ";
    ss << ev1.GetAtomCount() << " vs. " << ev2.GetAtomCount();
    throw Error(ss.str());
  }

  return calc_rmsd_for_atom_lists(ev1.GetAtomList(), ev2.GetAtomList(), 
                                  transformation);
}

geom::Vec3 EigenVec3ToVec3(const EVec3 &vec)
{
  return geom::Vec3(vec.data());
}

geom::Mat3 EigenMat3ToMat3(const EMat3 &mat)
{
  geom::Mat3 return_mat;
  for(int i=0;i<3;++i){
    for(int j=0;j<3;++j){
      return_mat(i,j) = mat(i,j);
    }
  }
  return return_mat;
}

geom::Mat4 EigenMat4ToMat4(const EMat4 &mat)
{ 
  geom::Mat4 return_mat;
  for(int i=0;i<4;++i){
    for(int j=0;j<4;++j){
      return_mat(i,j) = mat(i,j);
    }
  }
  return return_mat;
}

EMatX Mat4ToEigenMat4(const geom::Mat4 &mat){
  EMat4 res = EMat4::Zero();
  for(int i=0;i<4;++i){
    for(int j=0;j<4;++j){
      res(i,j)=mat.At(i,j);
    }
  }
  return res;
}

EVec3 Vec3ToEigenRVec(const geom::Vec3 &vec)
{
  return EVec3(&vec[0]);
}

EVec3 Vec3ToEigenVec(const geom::Vec3 &vec)
{
  return EVec3(&vec[0]);
}

EMatX MatrixShiftedBy(EMatX mat, ERVecX vec)
{
  EMatX result = mat;
  for (int row=0;  row<mat.rows();++row) {
    result.row(row) -= vec;
  }
  return result;
}


class MeanSquareMinimizerImpl {
public:
  MeanSquareMinimizerImpl(int n_atoms, bool alloc_atoms):
    n_atoms_(n_atoms), alloc_atoms_(alloc_atoms) 
  {
    if (alloc_atoms_) {
      atoms1_=EMatX::Zero(n_atoms_, 3);
      atoms2_=EMatX::Zero(n_atoms_, 3);
    }
  }

  void SetRefPos(size_t index, const geom::Vec3& pos) {
    atoms2_.row(index) = ERVec3(Vec3ToEigenVec(pos));
  }

  void SetPos(size_t index, const geom::Vec3& pos) {
    atoms1_.row(index) = ERVec3(Vec3ToEigenVec(pos));
  }

  SuperpositionResult Minimize(const EMatX& atoms, const EMatX& atoms_ref) const;

  EMatX TransformEMatX(const EMatX& mat, const EMat4& transformation) const;

  void CreateMatchingSubsets(std::vector<int>& matching_subset, const EMatX& atoms, 
                             const EMatX& atoms_ref, Real distance_threshold) const;

  SuperpositionResult IterativeMinimize(int ncycles, Real distance_threshold) const;

  SuperpositionResult MinimizeOnce() const;

private:
  int   n_atoms_;
  bool  alloc_atoms_;
  EMatX atoms1_;
  EMatX atoms2_;
};


MeanSquareMinimizer MeanSquareMinimizer::FromAtomLists(const AtomViewList& atoms,
                                                       const AtomViewList& atoms_ref) 
{
  int n_atoms = atoms.size();
  int n_atoms_ref = atoms_ref.size();
  if (n_atoms != n_atoms_ref) {
    throw Error("atom counts do not match");
  }
  if (n_atoms<3) {
    throw Error("at least 3 atoms are required for superposition");
  }
  MeanSquareMinimizer msm;
  msm.impl_ = new MeanSquareMinimizerImpl(n_atoms, true);

  for (size_t i = 0; i < atoms.size(); ++i ) {
    msm.impl_->SetRefPos(i, atoms_ref[i].GetPos());
    msm.impl_->SetPos(i, atoms[i].GetPos());
  }
  return msm;
}

MeanSquareMinimizer MeanSquareMinimizer::FromPointLists(const std::vector<geom::Vec3>& points,
                                                        const std::vector<geom::Vec3>& points_ref) 
{
  int n_points = points.size();
  int n_points_ref = points_ref.size();
  if (n_points != n_points_ref) {
    throw Error("point counts do not match");
  }
  if (n_points<3) {
    throw Error("at least 3 points are required for superposition");
  }
  MeanSquareMinimizer msm;
  msm.impl_ = new MeanSquareMinimizerImpl(n_points, true);

  for (size_t i = 0; i < points.size(); ++i ) {
    msm.impl_->SetRefPos(i, points_ref[i]);
    msm.impl_->SetPos(i, points[i]);
  }
  return msm;
}

SuperpositionResult MeanSquareMinimizerImpl::MinimizeOnce() const{
  return this->Minimize(atoms1_,atoms2_);
}

SuperpositionResult MeanSquareMinimizerImpl::IterativeMinimize(int max_cycles, Real distance_threshold) const{

  EMat4 transformation_matrix_old;
  EMat4 transformation_matrix;
  EMatX transformed_atoms;
  EMatX subset1 = atoms1_;
  EMatX subset2 = atoms2_;
  SuperpositionResult res;
  EMat4 diff;
  std::vector<int> matching_indices;

  //do initial superposition

  res = this->Minimize(atoms1_, atoms2_);
  transformation_matrix = Mat4ToEigenMat4(res.transformation);
  transformation_matrix_old = transformation_matrix;

  //note, that the initial superposition is the first cycle...
  int cycles=1;

  for(;cycles<max_cycles;++cycles){

    transformed_atoms = this->TransformEMatX(atoms1_, transformation_matrix);
    this->CreateMatchingSubsets(matching_indices, transformed_atoms, atoms2_, distance_threshold);

    if(matching_indices.size()<3){
      std::stringstream ss;
      ss<<"Failed to Converge, less than 3 superposable atoms in iteration ";
      ss<<cycles<<"!";
      throw Error(ss.str());
    }

    subset1 = EMatX::Zero(matching_indices.size(),3);
    subset2 = EMatX::Zero(matching_indices.size(),3);

    int i = 0;
    std::vector<int>::iterator it = matching_indices.begin();

    for(; it!=matching_indices.end(); ++it, ++i){
      subset1.row(i) = atoms1_.row(*it);
      subset2.row(i) = atoms2_.row(*it);
    }

    res = this->Minimize(subset1,subset2);
    transformation_matrix = Mat4ToEigenMat4(res.transformation);

    diff = transformation_matrix_old-transformation_matrix;

    if(diff.array().abs().sum()<0.0001){
      ++cycles;
      break;
    }

    transformation_matrix_old = transformation_matrix;

  }
  
  res.rmsd_superposed_atoms = calc_rmsd_for_ematx(subset1, subset2, transformation_matrix);
  res.fraction_superposed = float(subset1.rows())/atoms1_.rows();

  res.transformation = EigenMat4ToMat4(transformation_matrix);
  res.ncycles=cycles;

  return res;

}

EMatX MeanSquareMinimizerImpl::TransformEMatX(const EMatX& mat, const EMat4& transformation) const {

  EMatX transformed_mat = EMatX::Zero(mat.rows(), 3);

  EMatX vector = EMatX::Zero(4,1);
  EMatX transformed_vector = EMatX::Zero(4,1);
  vector(3,0)=1;

  for(int i=0;i<mat.rows();++i){
    vector.block<3,1>(0,0)=mat.block<1,3>(i,0).transpose();
    transformed_vector = transformation*vector;
    transformed_mat.block<1,3>(i,0)=transformed_vector.block<3,1>(0,0).transpose();
  }

  return transformed_mat;
}

void MeanSquareMinimizerImpl::CreateMatchingSubsets(std::vector<int>& matching_subset, 
                                                    const EMatX& atoms, const EMatX& atoms_ref, 
                                                    Real distance_threshold) const{


  matching_subset.clear();
  EMatX diff = EMatX::Zero(atoms.rows(),atoms.cols());
  EMatX dist = EMatX::Zero(atoms.rows(),1);

  diff = atoms-atoms_ref;
  dist = (diff.array()*diff.array()).rowwise().sum(); 

  Real squared_dist_threshold = distance_threshold * distance_threshold;
  
  for(int i=0; i<dist.rows(); ++i){
    if(dist(i,0) <= squared_dist_threshold){
      matching_subset.push_back(i);
    }
  }
}


SuperpositionResult MeanSquareMinimizerImpl::Minimize(const EMatX& atoms, const EMatX& atoms_ref) const {
  ERVec3 avg = atoms.colwise().sum()/atoms.rows();
  ERVec3 avg_ref = atoms_ref.colwise().sum()/atoms_ref.rows();

  // SVD only determines the rotational component of the superposition
  // we need to manually shift the centers of the two point sets on onto 
  // origin

  EMatX atoms_shifted = MatrixShiftedBy(atoms, avg);
  EMatX atoms_ref_shifted = MatrixShiftedBy(atoms_ref, avg_ref).transpose();

  // determine rotational component
  Eigen::JacobiSVD<EMat3> svd(atoms_ref_shifted*atoms_shifted,Eigen::ComputeFullU | Eigen::ComputeFullV);
  EMatX matrixVT=svd.matrixV().transpose();

  //determine rotation
  Real detv=matrixVT.determinant();
  Real dett=svd.matrixU().determinant();
  Real det=detv*dett;
  EMat3 rotation;
  if (det<0) {
    EMat3 tmat=EMat3::Identity();
    tmat(2,2)=-1;
    rotation = (svd.matrixU()*tmat)*matrixVT;
  } else {
    rotation = svd.matrixU()*matrixVT;
  }
  
  SuperpositionResult res;

  geom::Vec3 shift = EigenVec3ToVec3(avg_ref);
  geom::Vec3 com_vec = -EigenVec3ToVec3(avg);
  geom::Mat3 rot = EigenMat3ToMat3(rotation);
  geom::Mat4 mat4_com, mat4_rot, mat4_shift;
  mat4_rot.PasteRotation(rot);
  mat4_shift.PasteTranslation(shift);
  mat4_com.PasteTranslation(com_vec);
  //save whole transformation
  res.transformation = geom::Mat4(mat4_shift*mat4_rot*mat4_com);
  return res;
}

MeanSquareMinimizer& MeanSquareMinimizer::operator=(const MeanSquareMinimizer& rhs) {
  MeanSquareMinimizer tmp(rhs);
  this->swap(tmp);
  return *this;
}

MeanSquareMinimizer::MeanSquareMinimizer(const MeanSquareMinimizer& rhs) {
  if (rhs.impl_)
    impl_ = new MeanSquareMinimizerImpl(*rhs.impl_);
  else
    impl_ = NULL;
}


MeanSquareMinimizer::~MeanSquareMinimizer() {
  if (impl_) delete impl_;
}

SuperpositionResult MeanSquareMinimizer::MinimizeOnce() const {
  return impl_->MinimizeOnce();
}

SuperpositionResult MeanSquareMinimizer::IterativeMinimize(int ncycles, Real distance_threshold) const {
  return impl_->IterativeMinimize(ncycles, distance_threshold);
}

SuperpositionResult SuperposeAtoms(const mol::AtomViewList& atoms1,
                                   const mol::AtomViewList& atoms2,
                                   bool apply_transform)
{
  MeanSquareMinimizer msm = MeanSquareMinimizer::FromAtomLists(atoms1, atoms2);
  SuperpositionResult result = msm.MinimizeOnce();
  
  result.ncycles=1;
  result.rmsd = calc_rmsd_for_atom_lists(atoms1, atoms2, result.transformation);
  result.rmsd_superposed_atoms = result.rmsd;
  result.fraction_superposed = 1.0;
  if (apply_transform) {
    mol::AtomView jv=atoms1.front();
    mol::XCSEditor ed=jv.GetEntity().GetHandle().EditXCS();
    ed.ApplyTransform(result.transformation);
  }
  return result;

}

SuperpositionResult SuperposeSVD(const mol::EntityView& ev1,
                                 const mol::EntityView& ev2,
                                 bool apply_transform) {
  AtomViewList atoms1 = ev1.GetAtomList();
  AtomViewList atoms2 = ev2.GetAtomList();

  SuperpositionResult result = SuperposeAtoms(atoms1, atoms2, apply_transform);
  result.entity_view1 = ev1;
  result.entity_view2 = ev2;
  return result;
}

SuperpositionResult SuperposeSVD(const geom::Vec3List& pl1,
                                 const geom::Vec3List& pl2)
{
  MeanSquareMinimizer msm = MeanSquareMinimizer::FromPointLists(pl1, pl2);
  SuperpositionResult result = msm.MinimizeOnce();
  
  result.ncycles=1;
  result.rmsd = calc_rmsd_for_point_lists(pl1, pl2, result.transformation);
  return result;
}

SuperpositionResult IterativeSuperposeSVD(const mol::EntityView& ev,
                                          const mol::EntityView& ev_ref,
                                          int max_cycles,
                                          Real distance_threshold,
                                          bool apply_transform){

  AtomViewList atoms = ev.GetAtomList();
  AtomViewList atoms_ref = ev_ref.GetAtomList();

  MeanSquareMinimizer msm = MeanSquareMinimizer::FromAtomLists(atoms, atoms_ref);

  SuperpositionResult result = msm.IterativeMinimize(max_cycles, distance_threshold);

  result.rmsd = calc_rmsd_for_atom_lists(atoms, atoms_ref, result.transformation);
  if (apply_transform) {
    mol::AtomView jv=atoms.front();
    mol::XCSEditor ed=jv.GetEntity().GetHandle().EditXCS();
    ed.ApplyTransform(result.transformation);
  }

  result.entity_view1 = ev;
  result.entity_view2 = ev_ref;
  return result;
}

SuperpositionResult IterativeSuperposeSVD(const geom::Vec3List& pl1,
                                          const geom::Vec3List& pl2,
                                          int max_cycles,
                                          Real distance_threshold){

  MeanSquareMinimizer msm = MeanSquareMinimizer::FromPointLists(pl1, pl2);
  SuperpositionResult result = msm.IterativeMinimize(max_cycles, distance_threshold);
  result.rmsd = calc_rmsd_for_point_lists(pl1, pl2, result.transformation);
  return result;
}


}}} //ns
