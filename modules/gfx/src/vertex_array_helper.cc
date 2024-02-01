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
  Author: Ansgar Philippsen
*/

#include "vertex_array_helper.hh"

namespace ost { namespace gfx { namespace detail {

// TODO: refactor creation to geom, keep static list here
PrebuildSphereEntry& GetPrebuildSphere(unsigned int level)
{
  static std::vector<PrebuildSphereEntry> prebuild_list(VA_SPHERE_MAX_DETAIL+1);

  if(prebuild_list[level].vlist.empty()) {
    unsigned int step1 = level*4;
    unsigned int step2 = level*2;
    float fact1 = M_PI*2.0/static_cast<float>(step1);
    float fact2 = M_PI/static_cast<float>(step2);
    // generate points on the sphere
    std::vector<geom::Vec3>& vlist=prebuild_list[level].vlist;
    vlist.push_back(geom::Vec3(0.0,1.0,0.0));
    for(unsigned int i=0;i<step1;++i) {
      for(unsigned int k=1;k<step2;++k) {
        float ii = M_PI+fact1*static_cast<float>(i);
        float kk = fact2*static_cast<float>(k);
        vlist.push_back(geom::Vec3(sin(kk)*cos(ii), cos(kk), sin(kk)*sin(ii)));
      }
    }
    vlist.push_back(geom::Vec3(0.0,-1.0,0.0));
    
    // now the triangle indices for the above vertices
    std::vector<unsigned int>& ilist=prebuild_list[level].ilist;
    
    unsigned int point_count=vlist.size();
    
    // top
    unsigned int j=1;
    for(unsigned int i=0;i<step1;++i) {
      unsigned int k=i+1;
      k%=step1;
      ilist.push_back(0);
      ilist.push_back(k*(step2-1)+j);
      ilist.push_back(i*(step2-1)+j);
    }
    
    // bottom
    j=step2-1;
    for(unsigned int i=0;i<step1;++i) {
      unsigned int k=i+1;
      k%=step1;
      ilist.push_back(point_count-1);
      ilist.push_back(i*(step2-1)+j);
      ilist.push_back(k*(step2-1)+j);
    }  
    
    // middle
    for(unsigned int i=0;i<step1;++i) {
      unsigned int k=i+1;
      k%=step1;
      for(unsigned j=1;j<step2-1;++j) {
        unsigned l=j+1;
        l%=(step2);
        // add 2 triangles, counterclockwise
        unsigned int p1=i*(step2-1)+j;
        unsigned int p2=i*(step2-1)+l;
        unsigned int p3=k*(step2-1)+j;
        unsigned int p4=k*(step2-1)+l;
        ilist.push_back(p1);
        ilist.push_back(p3);
        ilist.push_back(p2);
        ilist.push_back(p2);
        ilist.push_back(p3);
        ilist.push_back(p4);
      }
    }
  }
  return prebuild_list[level];
}

std::vector<geom::Vec3> GetPrebuildCyl(unsigned int level)
{
  static std::vector<std::vector<geom::Vec3> > prebuild_list(VA_CYL_MAX_DETAIL+1);
  
  if(prebuild_list[level].empty()) {
    unsigned int divs = (level+1)*4;
    float step = 2.0*M_PI/static_cast<float>(divs);
    std::vector<geom::Vec3>& vlist=prebuild_list[level];
    for(unsigned int i=0;i<divs;++i) {
      float angle=step*static_cast<float>(i);
      vlist.push_back(geom::Vec3(std::cos(angle),std::sin(angle),0.0));
    }
  }
  return prebuild_list[level];
}


void ico_sphere_subdivide(std::vector<geom::Vec3>& vlist,
                      const geom::Vec3& v1, const geom::Vec3& v2, const geom::Vec3& v3,
                      unsigned int l)
{
  if(l==0) {
    vlist.push_back(v1);
    vlist.push_back(v2);
    vlist.push_back(v3);
  } else {
    geom::Vec3 v12=geom::Normalize(0.5*(v1+v2));
    geom::Vec3 v23=geom::Normalize(0.5*(v2+v3));
    geom::Vec3 v31=geom::Normalize(0.5*(v3+v1));
    ico_sphere_subdivide(vlist,v1,v12,v31,l-1);
    ico_sphere_subdivide(vlist,v2,v23,v12,l-1);
    ico_sphere_subdivide(vlist,v3,v31,v23,l-1);
    ico_sphere_subdivide(vlist,v12,v23,v31,l-1);
  }
}

std::vector<geom::Vec3> GetPrebuildIcoSphere(unsigned int level)
{
  static std::vector<std::vector<geom::Vec3> > prebuild_list(VA_ICO_SPHERE_MAX_DETAIL+1);

    if(prebuild_list[level].empty()) {
      static float X = .525731112119133606;
      static float Z = .850650808352039932;
      
      static float vdata[12][3] = {    
        {-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},    
        {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},    
        {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0} 
      };
      
      static unsigned int tindices[20][3] = { 
        {0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},    
        {8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},    
        {7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6}, 
        {6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11} };
      
      std::vector<geom::Vec3>& tmp_list=prebuild_list[level];
      for (int i=0; i<20; ++i) {    
        geom::Vec3 v1(vdata[tindices[i][0]][0],
                      vdata[tindices[i][0]][1],
                      vdata[tindices[i][0]][2]);
        geom::Vec3 v2(vdata[tindices[i][1]][0],
                      vdata[tindices[i][1]][1],
                      vdata[tindices[i][1]][2]);
        geom::Vec3 v3(vdata[tindices[i][2]][0],
                      vdata[tindices[i][2]][1],
                      vdata[tindices[i][2]][2]);
        ico_sphere_subdivide(tmp_list,
                             Normalize(v1),Normalize(v3),Normalize(v2),
                             level);
      }
    }
    return prebuild_list[level];
}


} // detail ns    


}} // ns
