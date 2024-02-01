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
 * Author Niklaus Johner
 */
#include <stdexcept>
#include <ost/base.hh>
#include <ost/geom/vec3.hh>
#include <ost/base.hh>
#include <ost/geom/geom.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/view_op.hh>
#include "trajectory_analysis.hh"

namespace ost { namespace mol { namespace alg {

geom::Vec3List AnalyzeAtomPos(const CoordGroupHandle& traj, const AtomHandle& a1, unsigned int stride)
//This function extracts the position of an atom from a trajectory and returns it as a vector of geom::Vec3
{
  CheckHandleValidity(traj);
  geom::Vec3List pos;
  pos.reserve(ceil(traj.GetFrameCount()/float(stride)));
  int i1=a1.GetIndex();
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    pos.push_back(frame->GetAtomPos(i1));
  }
  return pos;
}

std::vector<Real> AnalyzeDistanceBetwAtoms(const CoordGroupHandle& traj, const AtomHandle& a1, const AtomHandle& a2, 
                                  unsigned int stride)
//This function extracts the distance between two atoms from a trajectory and returns it as a vector
{
  CheckHandleValidity(traj);
  std::vector<Real> dist;
  dist.reserve(ceil(traj.GetFrameCount()/float(stride)));
  int i1=a1.GetIndex();
  int i2=a2.GetIndex();
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    dist.push_back(frame->GetDistanceBetwAtoms(i1,i2));
  }
  return dist;
} 

std::vector<Real> AnalyzeAngle(const CoordGroupHandle& traj, const AtomHandle& a1, const AtomHandle& a2, 
                               const AtomHandle& a3, unsigned int stride)
//This function extracts the angle between three atoms from a trajectory and returns it as a vector
{
  CheckHandleValidity(traj);
  std::vector<Real> ang;
  ang.reserve(ceil(traj.GetFrameCount()/float(stride)));
  int i1=a1.GetIndex(),i2=a2.GetIndex(),i3=a3.GetIndex();
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    ang.push_back(frame->GetAngle(i1,i2,i3));
  }
  return ang;
}

std::vector<Real> AnalyzeDihedralAngle(const CoordGroupHandle& traj, const AtomHandle& a1, const AtomHandle& a2, 
                                  const AtomHandle& a3, const AtomHandle& a4, unsigned int stride)
//This function extracts the dihedral angle between four atoms from a trajectory and returns it as a vector
{
  CheckHandleValidity(traj);
  std::vector<Real> ang;
  ang.reserve(ceil(traj.GetFrameCount()/float(stride)));
  int i1=a1.GetIndex(),i2=a2.GetIndex(),i3=a3.GetIndex(),i4=a4.GetIndex();
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    ang.push_back(frame->GetDihedralAngle(i1,i2,i3,i4));
  }
  return ang;
}

geom::Vec3List AnalyzeCenterOfMassPos(const CoordGroupHandle& traj, const EntityView& sele,unsigned int stride)
//This function extracts the position of the CenterOfMass of a selection (entity view) from a trajectory 
//and returns it as a vector. 
  {
  CheckHandleValidity(traj);
  geom::Vec3List pos;
  pos.reserve(ceil(traj.GetFrameCount()/float(stride)));
  std::vector<unsigned long> indices;
  std::vector<Real> masses;
  GetIndicesAndMasses(sele, indices, masses);
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    pos.push_back(frame->GetCenterOfMassPos(indices,masses));
  }
  return pos;
}

std::vector<Real> AnalyzeDistanceBetwCenterOfMass(const CoordGroupHandle& traj, const EntityView& sele1,
                                    const EntityView& sele2, unsigned int stride)
//This function extracts the distance between the CenterOfMass of two selection (entity views) from a trajectory 
//and returns it as a vector. 
  {
  CheckHandleValidity(traj);
  std::vector<Real> dist;
  dist.reserve(ceil(traj.GetFrameCount()/float(stride)));
  std::vector<unsigned long> indices1,indices2;
  std::vector<Real> masses1,masses2;
  GetIndicesAndMasses(sele1, indices1, masses1);
  GetIndicesAndMasses(sele2, indices2, masses2);
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    dist.push_back(frame->GetDistanceBetwCenterOfMass(indices1,masses1,indices2,masses2));
  }
  return dist;
}

std::vector<Real> AnalyzeRMSD(const CoordGroupHandle& traj, const EntityView& reference_view,
                                    const EntityView& sele_view, unsigned int stride)
// This function extracts the rmsd between two entity views and returns it as a vector
// The views don't have to be from the same entity
// If you want to compare to frame i of the trajectory t, first use t.CopyFrame(i) for example:
// eh=io.LoadPDB(...),t=io.LoadCHARMMTraj(eh,...);Sele=eh.Select(...);t.CopyFrame(0);mol.alg.AnalyzeRMSD(t,Sele,Sele)
  {
  CheckHandleValidity(traj);
  int count_ref=reference_view.GetAtomCount();
  int count_sele=sele_view.GetAtomCount();
  if (count_ref!=count_sele){
    throw Error("atom counts of the two views are not equal");
  }  
  std::vector<Real> rmsd;
  rmsd.reserve(ceil(traj.GetFrameCount()/float(stride)));
  std::vector<unsigned long> sele_indices;
  std::vector<geom::Vec3> ref_pos;
  GetIndices(sele_view, sele_indices);
  GetPositions(reference_view, ref_pos);
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    rmsd.push_back(frame->GetRMSD(ref_pos,sele_indices));
  }
  return rmsd;
}
 
std::vector<Real>  AnalyzeMinDistance(const CoordGroupHandle& traj, const EntityView& view1,
                                      const EntityView& view2,unsigned int stride)
// This function extracts the minimal distance between two sets of atoms (view1 and view2) for
// each frame in a trajectory and returns it as a vector.
  {
  CheckHandleValidity(traj);
  if (!view1.HasAtoms()){
    throw Error("first EntityView is empty");
  }
  if (!view2.HasAtoms()){
    throw Error("second EntityView is empty");
  }  
  std::vector<Real> dist;
  dist.reserve(ceil(traj.GetFrameCount()/float(stride)));
  std::vector<unsigned long> indices1,indices2;
  GetIndices(view1, indices1);
  GetIndices(view2, indices2);
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    dist.push_back(frame->GetMinDistance(indices1,indices2));
  }
    return dist;
}
  
std::vector<Real> AnalyzeMinDistanceBetwCenterOfMassAndView(const CoordGroupHandle& traj, const EntityView& view_cm,
                                            const EntityView& view_atoms,unsigned int stride)
  // This function extracts the minimal distance between a set of atoms (view_atoms) and the center of mass
  // of a second set of atoms (view_cm) for each frame in a trajectory and returns it as a vector.
  {
  CheckHandleValidity(traj);
  if (!view_cm.HasAtoms()){
    throw Error("first EntityView is empty");
  }
  if (!view_atoms.HasAtoms()){
    throw Error("second EntityView is empty");
  } 
  std::vector<Real> dist, masses_cm;
  dist.reserve(ceil(traj.GetFrameCount()/float(stride)));
  std::vector<unsigned long> indices_cm,indices_atoms;
  GetIndicesAndMasses(view_cm, indices_cm,masses_cm);
  GetIndices(view_atoms, indices_atoms);
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    dist.push_back(frame->GetMinDistBetwCenterOfMassAndView(indices_cm, masses_cm, indices_atoms));
  }
  return dist;
  }

std::vector<Real> AnalyzeAromaticRingInteraction(const CoordGroupHandle& traj, const EntityView& view_ring1,
                                                   const EntityView& view_ring2,unsigned int stride)
  // This function is a crude analysis of aromatic ring interactions. For each frame in a trajectory, it calculates
  // the minimal distance between the atoms in one view and the center of mass of the other
  // and vice versa, and returns the minimum between these two minimal distances.
  // Concretely, if the two views are the heavy atoms of two rings, then it returns the minimal
  // center of mass - heavy atom distance betweent he two rings
  {
  CheckHandleValidity(traj);
  if (!view_ring1.HasAtoms()){
    throw Error("first EntityView is empty");
  }
  if (!view_ring2.HasAtoms()){
    throw Error("second EntityView is empty");
  } 
  std::vector<Real> dist, masses_ring1,masses_ring2;
  dist.reserve(ceil(traj.GetFrameCount()/float(stride)));
  std::vector<unsigned long> indices_ring1,indices_ring2;
  Real d1,d2;
  GetIndicesAndMasses(view_ring1, indices_ring1,masses_ring1);
  GetIndicesAndMasses(view_ring2, indices_ring2,masses_ring2);
  for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
    CoordFramePtr frame=traj.GetFrame(i);
    d1=frame->GetMinDistBetwCenterOfMassAndView(indices_ring1, masses_ring1, indices_ring2);
    d2=frame->GetMinDistBetwCenterOfMassAndView(indices_ring2, masses_ring2, indices_ring1);
    dist.push_back(std::min(d1,d2));
  }
  return dist;  
  }

  void AnalyzeAlphaHelixAxis(const CoordGroupHandle& traj, const EntityView& prot_seg, geom::Vec3List& directions,
                             geom::Vec3List& centers, unsigned int stride)
  //This function calculates the best fitting cylinder to the C-alpha atoms of an EntityView and returns
  //the geometric center as well as the axis of that cylinder. We take care to have the axis point towards
  //the last residue of the selection, usually the direction of the alpha-helix
  {
    CheckHandleValidity(traj);
    if (!prot_seg.HasAtoms()){
      throw Error("EntityView is empty");
    }
    std::vector<unsigned long> indices_ca;
    std::pair<geom::Line3, Real> cylinder;
    geom::Line3 axis;
    Real sign;
    directions.reserve(ceil(traj.GetFrameCount()/float(stride)));
    centers.reserve(ceil(traj.GetFrameCount()/float(stride)));
    GetCaIndices(prot_seg, indices_ca);
    unsigned int n_atoms=indices_ca.size();
    for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
      CoordFramePtr frame=traj.GetFrame(i);
      cylinder=frame->FitCylinder(indices_ca);
      axis=cylinder.first;
      sign=geom::Dot(axis.GetDirection(),(*frame)[indices_ca[n_atoms-1]]-axis.GetOrigin());
      sign=sign/fabs(sign);
      directions.push_back(sign*axis.GetDirection());
      centers.push_back(axis.GetOrigin());
    }
    return;
  }
 
  //std::vector<geom::Line3> AnalyzeBestFitLine(const CoordGroupHandle& traj, const EntityView& prot_seg,
  //                                               unsigned int stride)
  void AnalyzeBestFitLine(const CoordGroupHandle& traj, const EntityView& prot_seg, geom::Vec3List& directions,
                          geom::Vec3List& centers, unsigned int stride)
  {
    CheckHandleValidity(traj);
    if (!prot_seg.HasAtoms()){
      throw Error("EntityView is empty");
    }
    std::vector<unsigned long> indices_ca;
    geom::Line3 axis;
    directions.reserve(ceil(traj.GetFrameCount()/float(stride)));
    centers.reserve(ceil(traj.GetFrameCount()/float(stride)));
    GetIndices(prot_seg, indices_ca);
    for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
      CoordFramePtr frame=traj.GetFrame(i);
      axis=frame->GetODRLine(indices_ca);
      directions.push_back(axis.GetDirection());
      centers.push_back(axis.GetOrigin());
    }
    return;
  }
  
  void AnalyzeBestFitPlane(const CoordGroupHandle& traj, const EntityView& prot_seg, geom::Vec3List& normals,
                          geom::Vec3List& origins, unsigned int stride)
  {
    CheckHandleValidity(traj);
    if (!prot_seg.HasAtoms()){
      throw Error("EntityView is empty");
    }
    std::vector<unsigned long> indices_ca;
    geom::Plane best_plane;
    normals.reserve(ceil(traj.GetFrameCount()/float(stride)));
    origins.reserve(ceil(traj.GetFrameCount()/float(stride)));
    GetIndices(prot_seg, indices_ca);
    for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
      CoordFramePtr frame=traj.GetFrame(i);
      best_plane=frame->GetODRPlane(indices_ca);
      normals.push_back(best_plane.GetNormal());
      origins.push_back(best_plane.GetOrigin());
    }
    return;
  }
  
  std::vector<Real> AnalyzeHelicity(const CoordGroupHandle& traj, const EntityView& prot_seg,
                                    unsigned int stride)
  {
    CheckHandleValidity(traj);
    if (!prot_seg.HasAtoms()){
      throw Error("EntityView is empty");
    }
    std::vector<unsigned long> indices_c,indices_o, indices_n, indices_ca;
    std::vector<Real> helicity;
    helicity.reserve(ceil(traj.GetFrameCount()/float(stride)));
    GetCaCONIndices(prot_seg, indices_ca, indices_c, indices_o, indices_n);
    for (size_t i=0; i<traj.GetFrameCount(); i+=stride) {
      CoordFramePtr frame=traj.GetFrame(i);
      helicity.push_back(frame->GetAlphaHelixContent(indices_ca,indices_c,indices_o,indices_n));
    }
    return helicity;
  }

  //This function constructs mean structures from a trajectory
  EntityHandle CreateMeanStructure(const CoordGroupHandle& traj, const EntityView& selection,
                            int from, int to, unsigned int stride)
  {
    CheckHandleValidity(traj);
    if (to==-1)to=traj.GetFrameCount();
    unsigned int n_atoms=selection.GetAtomCount();
    if (to<from) {
      throw Error("to smaller than from");
    }
    unsigned int n_frames=ceil((to-from)/static_cast<Real>(stride));
    if (n_atoms==0){
      throw Error("EntityView is empty");
    }
    if (n_frames<=1) {
      throw Error("number of frames is too small");
    }
    std::vector<unsigned long> indices;
    std::vector<geom::Vec3> mean_positions;
    EntityHandle eh;
    eh=CreateEntityFromView(selection,1);
    GetIndices(selection,indices);
    mean_positions.assign(n_atoms,geom::Vec3(0.0,0.0,0.0));
    for (int i=from; i<to; i+=stride) {
      CoordFramePtr frame=traj.GetFrame(i);
      for (unsigned int j=0; j<n_atoms; ++j) {
        mean_positions[j]+=(*frame)[indices[j]];
      }
    }
    mol::XCSEditor edi=eh.EditXCS(mol::BUFFERED_EDIT);
    AtomHandleList atoms=eh.GetAtomList();
    for (unsigned int j=0; j<n_atoms; ++j) {
      edi.SetAtomPos(atoms[j],mean_positions[j]/Real(n_frames));
    }
    return eh;
  }

  Real AnalyzeRMSF(const CoordGroupHandle& traj, const EntityView& selection, int from, int to, unsigned int stride)
  // This function extracts the rmsf between two entity views and assigns it 
  // The views don't have to be from the same entity
  // If you want to compare to frame i of the trajectory t, first use t.CopyFrame(i) for example:
  // eh=io.LoadPDB(...),t=io.LoadCHARMMTraj(eh,...);Sele=eh.Select(...);t.CopyFrame(0);mol.alg.AnalyzeRMSD(t,Sele,Sele)
  {
    CheckHandleValidity(traj);
    if (to==-1)to=traj.GetFrameCount();
    unsigned int n_atoms=selection.GetAtomCount();
    if (to<from) {
      throw Error("to smaller than from");
    }
    unsigned int n_frames=ceil((to-from)/static_cast<Real>(stride));
    if (n_atoms==0){
      throw Error("EntityView is empty");
    }
    if (n_frames<=1) {
      throw Error("number of frames is too small");
    }
    Real rmsf=0.0;
    geom::Vec3 v;
    std::vector<unsigned long> sele_indices;
    std::vector<geom::Vec3> ref_pos(n_atoms,geom::Vec3(0.,0.,0.));
    GetIndices(selection, sele_indices);
    for (unsigned int j=0; j<n_atoms; ++j) {
      for (int i=from; i<to; i+=stride) {
        CoordFramePtr frame=traj.GetFrame(i);
        ref_pos[j]+=frame->GetAtomPos(sele_indices[j]);
      }
      ref_pos[j]/=n_frames;
    }
    for (int i=from; i<to; i+=stride) {
      CoordFramePtr frame=traj.GetFrame(i);
      for (unsigned int j=0; j<n_atoms; ++j) {
        v=frame->GetAtomPos(sele_indices[j])-ref_pos[j];
        rmsf+=geom::Dot(v,v);
      }
    }
    return sqrt(rmsf/float(n_atoms*n_frames));
  }
  
}}} //ns
