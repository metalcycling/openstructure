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

/*
  Author: Marco Biasini
 */

#include <limits>
#include <Eigen/Eigenvalues>
#include <ost/mol/mol.hh>
#include <ost/mol/bounding_box.hh>
namespace ost { namespace mol {

namespace {

typedef Eigen::Matrix3d EMat3;

class PrincipalAxisCalc : public EntityVisitor {
public:
  PrincipalAxisCalc(const EntityHandle& ent): 
    natoms_(ent.GetAtomCount())
  {
    mean_=ent.GetCenterOfAtoms();
    cov_.setZero();
  }
  PrincipalAxisCalc(const AtomHandleList& atoms):
    natoms_(atoms.size())
  { 
    cov_.setZero();    
    mean_=geom::Vec3();
    if (natoms_>0) {
      for (AtomHandleList::const_iterator i=atoms.begin(),
           e=atoms.end(); i!=e; ++i) {
        mean_+=i->GetPos();       
      }      
    }
    mean_/=natoms_;
  }
  PrincipalAxisCalc(const EntityView& ent): 
    natoms_(ent.GetAtomCount())
  {
    mean_=ent.GetCenterOfAtoms();
    cov_.setZero();
  }  
  virtual bool VisitAtom(const AtomHandle& atom)
  {
    this->UpdateCov(atom);
    return false;
  } 
  void UpdateCov(const AtomHandle& atom)
  {
    geom::Vec3 p=atom.GetPos();
    for (size_t i=0; i<3; ++i) {
      for (size_t j=0; j<3; ++j) {
        Real var=p[i]*p[j];
        cov_(i, j)+=var;
      }
    }
  }
  
  void CalcPrincipalAxis()
  {
    for (size_t i=0; i<3; ++i) {
      for (size_t j=0; j<3; ++j) {
        cov_(i, j)=cov_(i, j)/natoms_-mean_[i]*mean_[j];
      }
    }
    typedef Eigen::SelfAdjointEigenSolver<EMat3> EigenSolver;
    EigenSolver solver(cov_);
    EigenSolver::MatrixType evalues=solver.eigenvectors();
    for (size_t i=0; i<3; ++i) {
      for (size_t j=0; j<3; ++j) {
        principal_axes_(i, j)=evalues(i, j);
      }
    }    
  }
  virtual void OnExit() 
  {
    this->CalcPrincipalAxis();
  }
  
  const geom::Mat3& GetPrincipalAxes() const 
  {
    return principal_axes_;
  }
private:
  EMat3 cov_;
  geom::Vec3  mean_;
  geom::Mat3  principal_axes_;
  int natoms_;
}; 

class  Projector : public EntityVisitor 
{
public:
  Projector(const geom::Mat3& principal_axes): axes_(principal_axes)
  {
    for (size_t i=0; i<3; ++i) {
      intervals_[i].first=std::numeric_limits<Real>::max();
      intervals_[i].second=-std::numeric_limits<Real>::max();   
    }
  }
  typedef std::pair<Real, Real> Interval;
  virtual bool VisitAtom(const AtomHandle& atom)
  {
    this->UpdateInterval(atom);
    return false;
  }
  void UpdateInterval(const AtomHandle& atom)
  {
    geom::Vec3 tp=atom.GetPos();
    for (size_t i=0; i<3; ++i) {
      Real proj=geom::Dot(axes_.GetCol(i), tp);
      intervals_[i].first=std::min(intervals_[i].first, proj);
      intervals_[i].second=std::max(intervals_[i].second, proj);
    }    
  }
  const Interval& GetInterval(int i)
  {
    return intervals_[i];
  }  
  
 geom::Cuboid GetBoundingBox() const
  {
    Real a=(intervals_[0].first+intervals_[0].second)*.5;
    Real b=(intervals_[1].first+intervals_[1].second)*.5;
    Real c=(intervals_[2].first+intervals_[2].second)*.5;        
    geom::Vec3 center=a*axes_.GetCol(0)+b*axes_.GetCol(1)+c*axes_.GetCol(2);
    geom::CuboidAxis axis_a=geom::CuboidAxis(axes_.GetCol(0), 
                           (intervals_[0].second-intervals_[0].first)*.5);
    geom::CuboidAxis axis_b=geom::CuboidAxis(axes_.GetCol(1), 
                           (intervals_[1].second-intervals_[1].first)*.5);
    geom::CuboidAxis axis_c=geom::CuboidAxis(axes_.GetCol(2), 
                           (intervals_[2].second-intervals_[2].first)*.5);
    return geom::Cuboid(center, axis_a, axis_b, axis_c);
  }
private:
  geom::Mat3 axes_;
  Interval intervals_[3];
  
};

template <typename E>
geom::Cuboid do_bounding_box_create(E ent)
{
  PrincipalAxisCalc calc(ent);
  ent.Apply(calc);
  const geom::Mat3& p_axes=calc.GetPrincipalAxes();
  Projector projector(p_axes);
  ent.Apply(projector);
  return projector.GetBoundingBox();
}

}

geom::Cuboid BoundingBoxFromEntity(const EntityHandle& ent)
{
  return do_bounding_box_create<EntityHandle>(ent);
}



geom::Cuboid BoundingBoxFromEntity(const EntityView& ent)
{
  return do_bounding_box_create<EntityView>(ent);  
}

geom::Cuboid BoundingBoxFromAtoms(const AtomHandleList& atoms)
{
  PrincipalAxisCalc calc(atoms);
  for (AtomHandleList::const_iterator i=atoms.begin(),
       e=atoms.end(); i!=e; ++i) {
    calc.UpdateCov(*i);
  }
  const geom::Mat3& p_axes=calc.GetPrincipalAxes();
  Projector projector(p_axes);
  for (AtomHandleList::const_iterator i=atoms.begin(), 
       e=atoms.end(); i!=e; ++i) {
    projector.UpdateInterval(*i);
  }  
  return projector.GetBoundingBox();
}



}}
