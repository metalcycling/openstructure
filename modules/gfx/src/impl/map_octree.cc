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

#include "map_octree.hh"

namespace ost { namespace gfx { namespace impl {

MapOctree::MapOctree(const img::ImageHandle& ih):
  map_(ih),
  built_(false)
{}

void MapOctree::Initialize()
{
  built_=false;
  levels_.clear();
  this->BuildOctree();
}


void MapOctree::SetNewMap(const img::ImageHandle& ih)
{
  map_=ih;
}

bool MapOctree::IsMapManageable (const img::ImageHandle ih)
{
 bool manageable = true;
 if (ih.GetExtent().GetWidth() > 256 ||
     ih.GetExtent().GetHeight() > 256 ||
     ih.GetExtent().GetDepth() > 256 ) manageable=false;
 return manageable;
}


uint32_t MapOctree::GetNumNodesForLevel(uint8_t level) const
{
  img::Size size=map_.GetExtent().GetSize();
  OcRangeVector range_vec(size[0]-1, size[1]-1, size[2]-1);
  uint16_t mask=0;
  for (unsigned char i=0; i<level; ++i) {
    mask|=1 << i;
  }
  return (((range_vec.x & mask))+1)*(((range_vec.y) & mask)+1)*
         (((range_vec.z & mask))+1);
}

void MapOctree::BuildOctree()
{
  img::Size size=map_.GetExtent().GetSize();
  OcRangeVector range_vec(size[0]-1, size[1]-1, size[2]-1);
  img::RealSpatialImageState* p=NULL;
  p=dynamic_cast<img::RealSpatialImageState*>(map_.ImageStatePtr().get());
  assert(p && "Octree only supports real spatial images");
  OctreeNode dummy;
  this->BuildOctreeRec(range_vec, 0, p, map_.GetExtent(), dummy);
  built_=true;
}

std::pair<float,float> MapOctree::BuildOctreeRec(const OcRangeVector& range_vec,
                                                 uint16_t level,
                                                 img::RealSpatialImageState* map,
                                                 const img::Extent& ext,
                                                 OctreeNode& parent)
{
  static img::Point OFFSET[]={ 
    img::Point(0, 0, 1), img::Point(0, 1, 1), img::Point(1, 0, 1),
    img::Point(1, 1, 1), img::Point(1, 1, 0), img::Point(0, 1, 0),
    img::Point(1, 0, 0)
  };
  // determine highest-order bit in all 3 directions
  char highest_order_bit=-1;
  for (char i=15; i>=0; --i) {
    uint16_t m=1<<i;
    if (range_vec.x & m || range_vec.y & m || range_vec.z & m) {
      highest_order_bit=i;
      break;
    }
  } 

  if (highest_order_bit==-1) {
    parent.SetLeaf();
    img::Point mend=map->GetExtent().GetEnd();
    float val=map->Value(ext.GetStart());
    std::pair<float, float> minmax(val, val);
    for (int i=0; i<7; ++i) {
      img::Point p=OFFSET[i]+ext.GetStart();
      if (mend[0]>=p[0] && mend[1]>=p[1] && mend[2]>=p[2]) {
        val=map->Value(p);
        if (val<minmax.first) {
          minmax.first=val;
        } else if (minmax.second<val) {
          minmax.second=val;
        }
      }
    }
    return minmax;
  }
  assert(parent.IsLeaf()==false);
  if (levels_.size()<=level) {
    levels_.push_back(OcNodeEntryList());
    levels_.back().reserve(this->GetNumNodesForLevel(level));
    assert(level<9 && "MAP TOO BIG for Octree");
  }  
  // determine branch pattern
  uint16_t highest_order_mask=1 << highest_order_bit;
  bool branch_x=(range_vec.x & highest_order_mask)!= 0;
  bool branch_y=(range_vec.y & highest_order_mask)!= 0;
  bool branch_z=(range_vec.z & highest_order_mask)!= 0;
  int range_x[2][2];
  int range_y[2][2];
  int range_z[2][2];
  // update ranges. this is basically the same code as in BuildOctreeRec()
  if (branch_x) {
    uint16_t ems=ext.GetEnd()[0]-ext.GetStart()[0];
    range_x[0][0]=ext.GetStart()[0];
    range_x[0][1]=range_x[0][0]+highest_order_mask-1;
    range_x[1][0]=ext.GetEnd()[0]-(ems & ~highest_order_mask);
    range_x[1][1]=ext.GetEnd()[0];
  } else {
    range_x[0][0]=ext.GetStart()[0];
    range_x[0][1]=ext.GetEnd()[0];
  }
  if (branch_y) {
    uint16_t ems=ext.GetEnd()[1]-ext.GetStart()[1];
    range_y[0][0]=ext.GetStart()[1];
    range_y[0][1]=range_y[0][0]+highest_order_mask-1;
    range_y[1][0]=ext.GetEnd()[1]-(ems & ~highest_order_mask);
    range_y[1][1]=ext.GetEnd()[1];
  } else {
    range_y[0][0]=ext.GetStart()[1];
    range_y[0][1]=ext.GetEnd()[1];
  }
  if (branch_z) {
    uint16_t ems=ext.GetEnd()[2]-ext.GetStart()[2];
    range_z[0][0]=ext.GetStart()[2];
    range_z[0][1]=range_z[0][0]+highest_order_mask-1;
    range_z[1][0]=ext.GetEnd()[2]-(ems & ~highest_order_mask);
    range_z[1][1]=ext.GetEnd()[2];
  } else {
    range_z[0][0]=ext.GetStart()[2];
    range_z[0][1]=ext.GetEnd()[2];
  }
  std::pair<float, float> minmax;
  uint32_t fc=0;
  if (levels_.size()>uint32_t(level+1)) {
    fc=levels_[level+1].size();
  }
  levels_[level].push_back(OctreeNode(fc, branch_x, branch_y, branch_z));
  for (int i=0; i<=int(branch_x); ++i) {
    for (int j=0; j<=int(branch_y); ++j) {
      for (int k=0; k<=int(branch_z); ++k) {        
        std::pair<float, float> local_minmax;
        img::Extent ext(img::Point(range_x[i][0], range_y[j][0], 
                                     range_z[k][0]),
                         img::Point(range_x[i][1], range_y[j][1], 
                                     range_z[k][1]));
#if !defined(NDEBUG)                                     
        if (!map_.GetExtent().Contains(ext)) {
          std::cout << ext << " is not contained in " << map_.GetExtent() 
                    <<  "(x=" << i << ", y=" << j << ", z=" << k 
                    << ", branch_x=" << branch_x << ", branch_y=" << branch_y 
                    << ", branch_z=" << branch_z << ")\n";
          assert(0 && "Internal Octree Error");
        }          
#endif                           
        OcRangeVector range(range_x[i][1]-range_x[i][0],
                            range_y[j][1]-range_y[j][0],
                            range_z[k][1]-range_z[k][0]);
        
        local_minmax=this->BuildOctreeRec(range, level+1, map, ext, 
                                          levels_[level].back());
        if (i+j+k) {
          if (local_minmax.first<minmax.first) {
            minmax.first=local_minmax.first;
          }
          if (local_minmax.second>minmax.second) {
            minmax.second=local_minmax.second;
          }
        } else {
          minmax=local_minmax;
        }
      }
    }    
  }
  levels_[level].back().SetMin(minmax.first);
  levels_[level].back().SetMax(minmax.second);  
  return minmax;
}

}}}
