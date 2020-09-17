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
  Authors: Niklaus Johner, Ansgar Philippsen
*/

#include <ost/invalid_handle.hh>
#include <ost/integrity_error.hh>
#include <ost/log.hh>
#include <ost/mol/in_mem_coord_source.hh>
#include <ost/mol/view_op.hh>
#include <ost/mol/mol.hh>
#include "coord_frame.hh"
#include <ost/base.hh>

namespace ost { namespace mol {

  namespace {
    /*
      Returns a score between 0 and 1 measuring the distance between
      a pair of angles and a pair of reference angles
    */
    Real Eval2AngleDist(Real phi, Real phi0, Real psi, Real psi0, Real delta)
    {
      Real d1=std::abs(phi-phi0);
      Real d2=std::abs(psi-psi0);
      if (d1>M_PI) d1=std::abs(d1-2.*M_PI);
      if (d2>M_PI) d1=std::abs(d2-2.*M_PI);
      Real d=(d1*d1+d2*d2)/(delta*delta);
      return 1.0/(1.0+d);
    }
  }

  geom::Vec3List CoordFrame::GetCellVectors() const
  {
    return geom::CalculateUnitCellVectors(ucell_size_,ucell_angles_);
  }
  
  geom::Vec3 CoordFrame::GetAtomPos(const AtomHandle& atom) const 
  {
    return this->GetAtomPos(atom.GetIndex());
  }

  geom::Vec3 CoordFrame::GetAtomPos(int i1) const 
  {
    return (*this)[i1];
  }
  
  Real CoordFrame::GetDistanceBetwAtoms(const AtomHandle& a1, const AtomHandle& a2) const
  {
    return this->GetDistanceBetwAtoms(a1.GetIndex(),a2.GetIndex());
  }

  Real CoordFrame::GetDistanceBetwAtoms(int i1, int i2) const
  {
    return geom::Distance((*this)[i1],(*this)[i2]);
  }
  
  Real CoordFrame::GetAngle(const AtomHandle& a1, const AtomHandle& a2, const AtomHandle& a3) const
  {
    return this->GetAngle(a1.GetIndex(),a2.GetIndex(),a3.GetIndex());  
  }

  Real CoordFrame::GetAngle(int i1, int i2, int i3) const 
  {
    return geom::Angle((*this)[i1]-(*this)[i2],(*this)[i3]-(*this)[i2]);
  }
  
  Real CoordFrame::GetDihedralAngle(const AtomHandle& a1, const AtomHandle& a2, 
                                    const AtomHandle& a3, const AtomHandle& a4) const
  {
    return this->GetDihedralAngle(a1.GetIndex(),a2.GetIndex(),a3.GetIndex(),a4.GetIndex());
  }

  Real CoordFrame::GetDihedralAngle(int i1, int i2, int i3, int i4) const
  {
    return geom::DihedralAngle((*this)[i1],(*this)[i2],(*this)[i3],(*this)[i4]);
  }

  
  geom::Vec3 CoordFrame::GetCenterOfMassPos(const EntityView& sele) const
  {
    std::vector<unsigned long> indices;
    std::vector<Real> masses;
    GetIndicesAndMasses(sele,indices,masses);
    return this->GetCenterOfMassPos(indices,masses);
  }
  
  geom::Vec3 CoordFrame::GetCenterOfMassPos(std::vector<unsigned long>& indices,
                                            std::vector<Real>& masses) const
  {
    geom::Vec3 v;
    for (unsigned int i=0,e=indices.size();i!=e; i++) {
      v+=masses[i]*(*this)[indices[i]];
    }
    return v;
  }
  
  Real CoordFrame::GetDistanceBetwCenterOfMass(const EntityView& sele1,const EntityView& sele2) const 
  {
    std::vector<unsigned long> indices1,indices2;
    std::vector<Real> masses1,masses2;
    GetIndicesAndMasses(sele1,indices1,masses1);
    GetIndicesAndMasses(sele2,indices2,masses2);
    return this->GetDistanceBetwCenterOfMass(indices1,masses1,indices2,masses2);
  }
  
  Real CoordFrame::GetDistanceBetwCenterOfMass(std::vector<unsigned long>& indices1,
                                               std::vector<Real>& masses1,
                                               std::vector<unsigned long>& indices2,
                                               std::vector<Real>& masses2) const
  {
    geom::Vec3 v1=this->GetCenterOfMassPos(indices1, masses1);
    geom::Vec3 v2=this->GetCenterOfMassPos(indices2, masses2);
    return geom::Distance(v1,v2);
  }
  
  Real CoordFrame::GetRMSD(const std::vector<geom::Vec3>& ref_pos,
                           const std::vector<unsigned long>& indices_sele) const
  {
    Real rmsd=0.0,val;
    for (unsigned int i1=0; i1!=indices_sele.size(); ++i1) {
      geom::Vec3 av_sele=(*this)[indices_sele[i1]];
      geom::Vec3 av_ref=ref_pos[i1];
      val=geom::Length2(av_ref-av_sele);
      rmsd+=val;
    }
    return sqrt(rmsd/indices_sele.size());
  }
  
  Real CoordFrame::GetRMSD(const EntityView& reference_view,const EntityView& sele_view) const
  {
    int count_ref=reference_view.GetAtomCount();
    int count_sele=sele_view.GetAtomCount();
    if (count_ref!=count_sele){
      throw Error("atom counts of the two views are not equal");
    }
    std::vector<unsigned long> indices_sele;
    std::vector<geom::Vec3> ref_pos;
    GetIndices(sele_view,indices_sele);
    GetPositions(reference_view,ref_pos);
    return this->GetRMSD(ref_pos,indices_sele);
  }
  
  Real CoordFrame::GetMinDistance(std::vector<unsigned long>& index_list1, 
                                  std::vector<unsigned long>& index_list2) const
  {
    geom::Vec3List pos_list1,pos_list2;
    for (std::vector<unsigned long>::const_iterator i1=index_list1.begin(),e=index_list1.end(); i1!=e; ++i1) {
      pos_list1.push_back((*this)[*i1]);
    }
    for (std::vector<unsigned long>::const_iterator i1=index_list2.begin(),e=index_list2.end(); i1!=e; ++i1) {
      pos_list2.push_back((*this)[*i1]);
    }
    return geom::MinDistance(pos_list1,pos_list2);
  }

  Real CoordFrame::GetMinDistance(const mol::EntityView& view1, const mol::EntityView& view2) const
  {
    std::vector<unsigned long> index_list1,index_list2;
    GetIndices(view1,index_list1);
    GetIndices(view2,index_list2);
    return this->GetMinDistance(index_list1,index_list2);
  }
  
  Real CoordFrame::GetMinDistBetwCenterOfMassAndView(std::vector<unsigned long>& indices_cm, 
                                                     std::vector<Real>& masses_cm,
                                                     std::vector<unsigned long>& indices_atoms) const
  {
    geom::Vec3List cm_pos,atoms_pos_list;
    cm_pos.push_back(this->GetCenterOfMassPos(indices_cm, masses_cm));
    for (std::vector<unsigned long>::const_iterator i1=indices_atoms.begin(),e=indices_atoms.end(); i1!=e; ++i1) {
      atoms_pos_list.push_back((*this)[*i1]);
    }
    return geom::MinDistance(cm_pos,atoms_pos_list);
  }

  Real CoordFrame::GetMinDistBetwCenterOfMassAndView(const mol::EntityView& view_cm, 
                                                     const mol::EntityView& view_atoms) const
  {
    std::vector<unsigned long> indices_cm,indices_atoms;
    std::vector<Real> masses_cm;
    GetIndices(view_atoms,indices_atoms);
    GetIndicesAndMasses(view_cm,indices_cm,masses_cm);
    return this->GetMinDistBetwCenterOfMassAndView(indices_cm,masses_cm,indices_atoms);
  }

  geom::Line3 CoordFrame::GetODRLine(std::vector<unsigned long>& indices_ca) const
  {
    geom::Vec3List atoms_pos_list;
    atoms_pos_list.reserve(indices_ca.size());
    for (std::vector<unsigned long>::const_iterator i1=indices_ca.begin(),e=indices_ca.end(); i1!=e; ++i1) {
      atoms_pos_list.push_back((*this)[*i1]);
    }
    return atoms_pos_list.GetODRLine();
  }

  geom::Line3 CoordFrame::GetODRLine(const mol::EntityView& view1) const {
    std::vector<unsigned long> indices;
    GetIndices(view1,indices);
    return this->GetODRLine(indices);
  }
  
  geom::Plane CoordFrame::GetODRPlane(std::vector<unsigned long>& indices_ca) const
  {
    geom::Vec3List atoms_pos_list;
    atoms_pos_list.reserve(indices_ca.size());
    for (std::vector<unsigned long>::const_iterator i1=indices_ca.begin(),e=indices_ca.end(); i1!=e; ++i1) {
      atoms_pos_list.push_back((*this)[*i1]);
    }
    return atoms_pos_list.GetODRPlane();
  }

  geom::Plane CoordFrame::GetODRPlane(const mol::EntityView& view1) const
  {
    std::vector<unsigned long> indices;
    GetIndices(view1,indices);
    return this->GetODRPlane(indices);
  }  
  
  std::pair<geom::Line3, Real> CoordFrame::FitCylinder(std::vector<unsigned long>& indices_ca) const {
    geom::Vec3List atoms_pos_list;
    int n_atoms=indices_ca.size();
    atoms_pos_list.reserve(n_atoms);
    for (std::vector<unsigned long>::const_iterator i1=indices_ca.begin(),e=indices_ca.end(); i1!=e; ++i1) {
      atoms_pos_list.push_back((*this)[*i1]);
    }
    //Initial guess
    geom::Vec3 initial_axis=geom::Vec3(0.0,0.0,0.0);
    if (n_atoms<5) {
      initial_axis=atoms_pos_list[n_atoms-1]-atoms_pos_list[0];
    }
    else {
      for (geom::Vec3List::const_iterator i=atoms_pos_list.begin(),e=atoms_pos_list.end()-4; i!=e; ++i) {
        initial_axis+=(*(i+4))-(*i);
      }
    }
    return atoms_pos_list.FitCylinder(initial_axis);
  }

  std::pair<geom::Line3, Real> CoordFrame::FitCylinder(const mol::EntityView& view1) const {
    CheckHandleValidity(view1);
    std::vector<unsigned long> indices_ca;
    GetCaIndices(view1, indices_ca);
    return CoordFrame::FitCylinder(indices_ca);
  }
  
  Real CoordFrame::GetAlphaHelixContent(const mol::EntityView& segment) const 
  {
    CheckHandleValidity(segment);
    std::vector<unsigned long> indices_c, indices_n, indices_ca, indices_o;
    GetCaCONIndices(segment, indices_ca, indices_c, indices_o, indices_n);
    return CoordFrame::GetAlphaHelixContent(indices_ca,indices_c,indices_o,indices_n);
  }

  Real CoordFrame::GetAlphaHelixContent(std::vector<unsigned long>& indices_ca, 
                                        std::vector<unsigned long>& indices_c,
                                         std::vector<unsigned long>& indices_o, 
                                        std::vector<unsigned long>& indices_n) const
  {
    //See CoordFrame::GetAlphaHelixContent(const mol::EntityView) above.
    unsigned long n_atoms=indices_ca.size();
    geom::Vec3 c_previous,n,ca,c,n_next;
    Real phi,psi,phi_0=-1.2,psi_0=-0.785,delta=0.8,d_0=3.3;
    unsigned long n_helical_res=0;
    std::vector<Real> score,dist,score2,dist2;
    score.reserve(n_atoms-2);
    score2.reserve(n_atoms-2);
    dist.reserve(n_atoms-2);
    dist2.reserve(n_atoms-2);
    if (n_atoms!=indices_n.size()||n_atoms!=indices_c.size()||n_atoms!=indices_o.size()){
      throw Error("not same numbers of CA, C, O and N atoms in the selection");
    }
    if (n_atoms<=5){
      throw Error("At least five residues are needed to calculate an alpha helix similarity");
    }
    c=(*this)[indices_c[0]];
    n_next=(*this)[indices_n[1]];
    for (unsigned long i=1; i!=n_atoms-1; ++i) {
      c_previous=c;
      n=n_next;
      ca=(*this)[indices_ca[i]];
      c=(*this)[indices_c[i]];
      n_next=(*this)[indices_n[i+1]];
      phi=geom::DihedralAngle(c_previous,n,ca,c);
      psi=geom::DihedralAngle(n,ca,c,n_next);
      score.push_back(Eval2AngleDist(phi,phi_0,psi,psi_0,delta));
    }
    score2[0]=sqrt(score[0]*score[1]*score[0]*score[1]*score[0]*score[1]);
    score2[n_atoms-3]=sqrt(score[n_atoms-3]*score[n_atoms-4]*
                           score[n_atoms-3]*score[n_atoms-4]*
                           score[n_atoms-3]*score[n_atoms-4]);
    for (unsigned long i=1; i!=n_atoms-3; ++i) {
      score2[i]=score[i-1]*score[i]*score[i+1];
    }
    if (n_atoms>=8){
      for (unsigned long i=1; i!=n_atoms-1; ++i){
        if (i<n_atoms-4){
          dist.push_back(geom::Distance((*this)[indices_o[i]],(*this)[indices_n[i+4]]));
          if (i>=5) {
            dist2.push_back(std::min(dist[i-1],dist[i-5]));
          } else {
            dist2.push_back(dist[i-1]);
          }
        } else {
          dist2.push_back(geom::Distance((*this)[indices_n[i]],(*this)[indices_o[i-4]]));
        }
      }
      for (unsigned long i=0; i!=n_atoms-2; ++i){
        if (dist2[i]<=d_0 || score2[i]>=0.3) {
          n_helical_res+=1;
        }
      }
    } else {
      for (unsigned long i=0; i!=n_atoms-2; ++i){
        if (score2[i]>=0.3) {
          n_helical_res+=1;
        }
      }
    }
    return Real(n_helical_res)/Real(n_atoms-2);
  }

  CoordFrame CreateCoordFrame(const geom::Vec3List& atom_pos, const geom::Vec3& cs, const geom::Vec3& ca)
  {
    return CoordFrame(atom_pos,cs,ca);
  }

  void GetIndices(const EntityView& sele, std::vector<unsigned long>& indices)
  {
    AtomViewList atoms=sele.GetAtomList();
    indices.reserve(atoms.size());
    for (AtomViewList::const_iterator i=atoms.begin(), 
         e=atoms.end(); i!=e; ++i) {
      indices.push_back(i->GetIndex());
    }
  }

  void GetCaIndices(const EntityView& segment, std::vector<unsigned long>& indices_ca){
    mol::AtomView ca;
    ResidueViewList residues=segment.GetResidueList();
    indices_ca.reserve(residues.size());
    //for (ResidueViewList::const_iterator res=residues.begin()+1,
    //     e=residues.end(); res!=e-1; ++res){
    //  if (!InSequence((*res).GetHandle(),(*(res+1)).GetHandle())) throw Error("Residues are not in a continuous sequence");
    //}
    for (ResidueViewList::const_iterator res=residues.begin(),
         e=residues.end(); res!=e; ++res) {
      ca=res->FindAtom("CA");
      CheckHandleValidity(ca);
      indices_ca.push_back(ca.GetIndex());
      }
  }
  
  void GetCaCONIndices(const EntityView& segment, std::vector<unsigned long>& indices_ca,
                      std::vector<unsigned long>& indices_c, std::vector<unsigned long>& indices_o, std::vector<unsigned long>& indices_n){
    mol::AtomView a,n,c,o;
    ResidueViewList residues=segment.GetResidueList();
    indices_ca.reserve(residues.size());
    indices_n.reserve(residues.size());
    indices_c.reserve(residues.size());
    indices_o.reserve(residues.size());
    for (ResidueViewList::const_iterator res=residues.begin()+1,
         e=residues.end(); res!=e-1; ++res){
      if (!InSequence((*res).GetHandle(),(*(res+1)).GetHandle())) throw Error("Residues are not in a continuous sequence");
    }
    for (ResidueViewList::const_iterator res=residues.begin(),
         e=residues.end(); res!=e; ++res) {
      a=res->FindAtom("CA");
      CheckHandleValidity(a);
      indices_ca.push_back(a.GetIndex());
      c=res->FindAtom("C");
      CheckHandleValidity(c);
      indices_c.push_back(c.GetIndex());
      o=res->FindAtom("O");
      CheckHandleValidity(o);
      indices_o.push_back(o.GetIndex());
      n=res->FindAtom("N");
      CheckHandleValidity(n);
      indices_n.push_back(n.GetIndex());
    }
  }

  void GetMasses(const EntityView& sele, std::vector<Real>& masses){
    Real mass_tot=0.0;
    AtomViewList atoms=sele.GetAtomList();
    masses.reserve(atoms.size());
    for (AtomViewList::const_iterator i=atoms.begin(), 
         e=atoms.end(); i!=e; ++i) {
      masses.push_back(i->GetMass());
      mass_tot=mass_tot+i->GetMass();
    }
    for (std::vector<Real>::iterator 
         j=masses.begin(), e2=masses.end(); j!=e2; ++j) {
      (*j)/=mass_tot;
    }  
  }
  
  
  void GetIndicesAndMasses(const EntityView& sele, std::vector<unsigned long>& indices,std::vector<Real>& masses){
    GetIndices(sele, indices);
    GetMasses(sele, masses);
  }
  
  void GetPositions(const EntityView& sele, std::vector<geom::Vec3>& ref_pos){
    ref_pos.reserve(sele.GetAtomCount());
    for (ChainViewList::const_iterator ci = sele.GetChainList().begin(),
         ce = sele.GetChainList().end(); ci != ce; ++ci) {
      for (ResidueViewList::const_iterator ri = ci->GetResidueList().begin(),
           re = ci->GetResidueList().end(); ri != re; ++ri) {
        for (AtomViewList::const_iterator ai = ri->GetAtomList().begin(),
             ae = ri->GetAtomList().end(); ai != ae; ++ai) {
          ref_pos.push_back(ai->GetPos());
        }
      }
    }
  }

}}//ns

