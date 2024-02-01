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

#include <ost/base.hh>
#include <ost/mol/mol.hh>
#include "structure_analysis.hh"

namespace ost { namespace mol { namespace alg {

geom::Vec3List GetPosListFromView(const EntityView& view){
  CheckHandleValidity(view);
  geom::Vec3List vl;
  AtomViewList atoms=view.GetAtomList();
  vl.reserve(atoms.size());
  for (AtomViewList::const_iterator i=atoms.begin(), 
       e=atoms.end(); i!=e; ++i) {
    vl.push_back(i->GetPos());
  }
  return vl;
}    

std::vector<Real> CalculateAgreementWithDensityMap(const geom::Vec3List& vl, img::MapHandle& density_map){
  CheckHandleValidity(density_map);
  std::vector<Real> v;
  v.reserve(vl.size());
  for (geom::Vec3List::const_iterator v1=vl.begin(),e=vl.end(); v1!=e; ++v1) {
    img::Point p(density_map.CoordToIndex(*v1));
    v.push_back(density_map.GetReal(p));
  }
  return v;
}
  
Real CalculateAverageAgreementWithDensityMap(const geom::Vec3List& vl, img::MapHandle& density_map){
  CheckHandleValidity(density_map);
  std::vector<Real> v=CalculateAgreementWithDensityMap(vl, density_map);
  Real sum=0.0;
  for (std::vector<Real>::const_iterator i=v.begin(),e=v.end(); i!=e; ++i) {
    sum=sum+*i;
  }
  return sum/float(vl.size());
}
  
std::vector< std::vector<Real> > PariwiseDistanceMatrix(const EntityView& view1, const EntityView& view2){
  std::vector< std::vector<Real> > dist_mat;
  unsigned int n1=view1.GetAtomCount(),n2=view2.GetAtomCount();
  dist_mat.resize(n1,std::vector<Real>(n2,0));
  AtomViewList atoms1=view1.GetAtomList(),atoms2=view2.GetAtomList();
  for (unsigned int i=0; i!=n1; ++i){
    for (unsigned int j=0; j!=n2; ++j) {
      dist_mat[i][j]=geom::Distance(atoms1[i].GetPos(),atoms2[j].GetPos());
    }
  }
  return dist_mat;
}

std::vector< std::vector<Real> > PariwiseDistanceMatrix(const EntityView& view){
  std::vector< std::vector<Real> > dist_mat;
  unsigned int n=view.GetAtomCount();
  //dist_mat.resize(n,std::vector<Real>(n,0));
  AtomViewList atoms1=view.GetAtomList(),atoms2=view.GetAtomList();
  for (unsigned int i=0; i!=n; ++i){
    dist_mat.push_back(std::vector<Real>(n-i-1,0));
    for (unsigned int j=i+1; j!=n; ++j) {
      dist_mat[i][j-i-1]=geom::Distance(atoms1[i].GetPos(),atoms2[j].GetPos());
    }
  }
  return dist_mat;
}

void DLLEXPORT_OST_MOL_ALG WrapEntityInPeriodicCell(EntityHandle eh, const geom::Vec3 cell_center, const geom::Vec3 ucell_size, \
                                                    const geom::Vec3 ucell_angles, bool group_residues,bool follow_bonds){
  mol::XCSEditor edi=eh.EditXCS(mol::BUFFERED_EDIT);
  bool orthogonal;
  if (ucell_angles==geom::Vec3()){ orthogonal=true; }
  else { orthogonal=false; }
  if (group_residues) {
    geom::Vec3 cm,wrapped_cm,shift;
    ResidueHandleList residues=eh.GetResidueList();
    unsigned int n_residues=eh.GetResidueCount();
    for (unsigned int i=0; i<n_residues; ++i) {
      ResidueHandle r=residues[i];
      AtomHandleList atoms=r.GetAtomList();
      geom::Vec3 ref_pos=atoms[0].GetPos();
      if (follow_bonds) {
        int natoms=r.GetAtomCount(),n=1,cycles=1;
        AtomHandle a=atoms[0];
        AtomHandleList wrapped_atoms;
        wrapped_atoms.reserve(natoms);
        wrapped_atoms.push_back(a);
        while (n<natoms and cycles<natoms) {
          ++cycles;
          for (AtomHandleList::iterator a=wrapped_atoms.begin(), e=wrapped_atoms.end(); a!=e; ++a) {
            AtomHandleList al=(*a).GetBondPartners();
            ref_pos=(*a).GetPos();
            for (AtomHandleList::iterator a2=al.begin(), e2=al.end(); a2!=e2; ++a2) {
              if (std::find(wrapped_atoms.begin(),wrapped_atoms.end(), *a2) !=wrapped_atoms.end()){
                continue;
              } else {
                ++n;
                wrapped_atoms.push_back(*a2);
                if (orthogonal) {
                    edi.SetAtomPos((*a2),geom::WrapVec3((*a2).GetPos(),ref_pos,ucell_size));
                } else {
                  edi.SetAtomPos((*a2),geom::WrapVec3((*a2).GetPos(),ref_pos,ucell_size,ucell_angles));
          }}}}}
      } else {
        if (orthogonal) {
          for (AtomHandleList::iterator a=atoms.begin(), e=atoms.end(); a!=e; ++a) {
            edi.SetAtomPos((*a),geom::WrapVec3((*a).GetPos(),ref_pos,ucell_size));
          }}
        else {
          for (AtomHandleList::iterator a=atoms.begin(), e=atoms.end(); a!=e; ++a) {
            edi.SetAtomPos((*a),geom::WrapVec3((*a).GetPos(),ref_pos,ucell_size,ucell_angles));
        }}}
    }
    for (unsigned int i=0; i<n_residues; ++i) {
      ResidueHandle r=residues[i];
      cm=r.GetCenterOfMass();
      if (orthogonal) {wrapped_cm=geom::WrapVec3(cm,cell_center,ucell_size);}
      else {wrapped_cm=geom::WrapVec3(cm,cell_center,ucell_size,ucell_angles);}
      if (wrapped_cm==cm) continue;
      AtomHandleList atoms=r.GetAtomList();
      unsigned int n_atoms=r.GetAtomCount();
      shift=wrapped_cm-cm;
      for (unsigned int j=0; j<n_atoms; ++j) {
        edi.SetAtomPos(atoms[j],atoms[j].GetPos()+shift);
      }
    }
  }
  else {
    AtomHandleList atoms=eh.GetAtomList();
    if (orthogonal){
      for (AtomHandleList::iterator a=atoms.begin(), e=atoms.end(); a!=e; ++a) {
        edi.SetAtomPos((*a),geom::WrapVec3((*a).GetPos(),cell_center,ucell_size));
      }}
    else {
      for (AtomHandleList::iterator a=atoms.begin(), e=atoms.end(); a!=e; ++a) {
        edi.SetAtomPos((*a),geom::WrapVec3((*a).GetPos(),cell_center,ucell_size,ucell_angles));
      }}
  }
}

}}} //ns
