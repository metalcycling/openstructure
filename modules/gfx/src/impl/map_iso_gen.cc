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

#include <iostream>
#include <vector>

#include <ost/gfx/module_config.hh>
#include <ost/log.hh>
#include <ost/geom/geom.hh>

#include <ost/img/image_state.hh>

//#define MAP_ISO_DEBUG

#include <ost/gfx/map_iso_prop.hh>
#include "map_iso_spec.hh"
#include "map_iso_gen.hh"

namespace ost { namespace gfx { namespace impl {

// singleton implementation
IsosurfaceGenerator& IsosurfaceGenerator::Instance()
{
  static IsosurfaceGenerator inst;
  return inst;
}

IsosurfaceGenerator::IsosurfaceGenerator():
  add_lf_func_()
{
  map_iso::GenerateLFList(add_lf_func_);
}

/*
  Main isosurfacing code
  
  This is in essence the marching cubes algorithm, but with a sophisticated vertex id 
  book-keeping routine to avoid recalculation of interpolated isosurface values and
  reduce vertex count during rendering. Some other speed optimization are also included.

  Wether a vertex id is actually set/used is ignored at this high level, it is less
  expensive to blindly assign integers than to check wether they are actually used by
  the Add* functions.

  vertex id offsets per cube
       ---9---
      5.     6|
     / 8    / |
     ---1--   10
    |  ..11|..|
    0 4    2 7  
    |.     |/
     --3---

  corner offsets

      5        6
       -------
      /.     /|
     / .    / |
   1 ------ 2 |
    |  ....|..|
    | .4   | / 7
    |.     |/
     ------
   0        3

*/
namespace {

inline uint gen_vertex_0(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_1(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_2(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_3(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_4(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_5(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_6(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_7(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_8(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_9(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_A(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);
inline uint gen_vertex_B(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level);


}

void IsosurfaceGenerator::Generate(const img::MapHandle& mh, IndexedVertexArray& va, float level)
{
#ifdef MAP_ISO_DEBUG
  va.Add(geom::Vec3(0.0,0.0,0.0),geom::Vec3(),Color(1.0,1.0,1.0));
#endif
  img::ImageStateBasePtr is=mh.ImageStatePtr();
  img::RealSpatialImageState* rsis=dynamic_cast<img::RealSpatialImageState*>(is.get());
  assert(rsis);
  // TODO: use *start to *end instead 0 to *count
  int ucount=mh.GetSize()[0];
  int vcount=mh.GetSize()[1];
  int wcount=mh.GetSize()[2];
  img::Point s=mh.GetExtent().GetStart();
  uint vertex_id[12];
  
  // storage for vertex ids for the last xy plane
  // _x_ for positions 0 and 2
  std::vector<uint> vid_uvplane_x_1(ucount*(vcount-1));
  std::vector<uint> vid_uvplane_x_2(ucount*(vcount-1));
  // _y_ for positions 1 and 3
  std::vector<uint> vid_uvplane_y_1((ucount-1)*vcount);
  std::vector<uint> vid_uvplane_y_2((ucount-1)*vcount);

  // storage for vertex ids for the last y row
  // _x_ for positions 5 and 6
  std::vector<uint> vid_vrow_x_1(ucount);
  std::vector<uint> vid_vrow_x_2(ucount);
  // _z_ for position 9 (1 already covered by uvplane)
  std::vector<uint> vid_vrow_z_1(ucount-1);
  std::vector<uint> vid_vrow_z_2(ucount-1);

  // and finally for positions 4 and 8
  uint vid_ucube_1[2];

  // precalc first entries for uvplane_x
  for(int v=0;v<vcount-1;++v) {
    for(int u=0;u<ucount;++u) {
      vid_uvplane_x_1[v*ucount+u]=gen_vertex_0(rsis,va,u+s[0],
                                               v+s[1],s[2],level);
    }
  }

  // precalc first entries for uvplane_y
  for(int v=0;v<vcount;++v) {
    for(int u=0;u<ucount-1;++u) {
      vid_uvplane_y_1[v*(ucount-1)+u]=gen_vertex_3(rsis,va,u+s[0],v+s[1],
                                                   s[2],level);
    }
  }

  // main loop over planes
  for(int w=0;w<wcount-1;++w) {
    // precalc first entries for vrow_x
    for(int u=0;u<ucount;++u) {
      vid_vrow_x_1[u]=gen_vertex_4(rsis,va,u+s[0],s[1],w+s[2],level);
    }
    /*
      precalc first entries for vrow_z
      these are transfered to uvplane below
    */
    for(int u=0;u<ucount-1;++u) {
      vid_vrow_z_1[u]=gen_vertex_B(rsis,va,u+s[0],s[1],w+s[2],level);
    }

    // main loop over rows
    for(int v=0;v<vcount-1;++v) {
      /*
  precalc the first u entries
  these are transfered to vrow and uvplane below
      */
      vid_ucube_1[0]=gen_vertex_5(rsis,va,s[0],v+s[1],w+s[2],level);
      vid_ucube_1[1]=gen_vertex_8(rsis,va,s[0],v+s[1],w+s[2],level);

      // main loop over cubes
      for(int u=0;u<ucount-1;++u) {
        int us=s[0]+u, vs=s[1]+v, ws=s[2]+w;
        // determine cube bitpattern
        unsigned int pattern = 0;
        if(static_cast<float>(rsis->Value(img::Point(us,vs,ws)))<level) pattern|=1<<0;
        if(static_cast<float>(rsis->Value(img::Point(us,vs+1,ws)))<level) pattern|=1<<1;
        if(static_cast<float>(rsis->Value(img::Point(us+1,vs+1,ws)))<level) pattern|=1<<2;
        if(static_cast<float>(rsis->Value(img::Point(us+1,vs,ws)))<level) pattern|=1<<3;
        if(static_cast<float>(rsis->Value(img::Point(us,vs,ws+1)))<level) pattern|=1<<4;
        if(static_cast<float>(rsis->Value(img::Point(us,vs+1,ws+1)))<level) pattern|=1<<5;
        if(static_cast<float>(rsis->Value(img::Point(us+1,vs+1,ws+1)))<level) pattern|=1<<6;
        if(static_cast<float>(rsis->Value(img::Point(us+1,vs,ws+1)))<level) pattern|=1<<7;
        
        LOG_TRACE(u << " " << v << " " << w << ": " << pattern << " " 
                  << mh.GetReal(img::Point(u,v,w)));
        
        /*
          set up already known vertex-ids from book keeping
          arrays, in fact 9 out of 12
        */
        vertex_id[0]=vid_uvplane_x_1[v*ucount+u];
        vertex_id[1]=vid_uvplane_y_1[(v+1)*(ucount-1)+u];
        vertex_id[2]=vid_uvplane_x_1[v*ucount+u+1];
        vertex_id[3]=vid_uvplane_y_1[v*(ucount-1)+u];
        
        vertex_id[4]=vid_vrow_x_1[u];
        vertex_id[7]=vid_vrow_x_1[u+1];
        vertex_id[11]=vid_vrow_z_1[u];
        
        vertex_id[5]=vid_ucube_1[0];
        vertex_id[8]=vid_ucube_1[1];
        
        // generate remaining ones
        vertex_id[6]=gen_vertex_6(rsis,va, us, vs, ws,level);
        vertex_id[9]=gen_vertex_9(rsis,va,us, vs, ws,level);
        vertex_id[10]=gen_vertex_A(rsis,va,u+s[0],v+s[1],w+s[2],level);
        
        // add lines and faces according to the particular cube type
        add_lf_func_[pattern](va,vertex_id);
        
        /*
          store old and new vertex_ids in book keeping arrays
          duplication is on purpose, because the vids are
          re-used at different times
        */
        vid_ucube_1[0]=vertex_id[6];
        vid_ucube_1[1]=vertex_id[10];
        
        vid_vrow_x_2[u]=vertex_id[5];
        vid_vrow_x_2[u+1]=vertex_id[6];
        vid_vrow_z_2[u]=vertex_id[9];
  
        vid_uvplane_x_2[v*ucount+u]=vertex_id[8];
        vid_uvplane_y_2[(v+1)*(ucount-1)+u]=vertex_id[9];
        vid_uvplane_x_2[v*ucount+u+1]=vertex_id[10];
        vid_uvplane_y_2[v*(ucount-1)+u]=vertex_id[11];
      } // u
      std::swap(vid_vrow_x_1,vid_vrow_x_2);
      std::swap(vid_vrow_z_1,vid_vrow_z_2);
    } // v
    std::swap(vid_uvplane_x_1,vid_uvplane_x_2);
    std::swap(vid_uvplane_y_1,vid_uvplane_y_2);
  } // w
}

void IsosurfaceGenerator::AddLF(uint code, IndexedVertexArray& va,unsigned int vertex_id[12])
{
  add_lf_func_[code](va,vertex_id);
}

namespace {

// lots of similar code, for performance purposes

static Color def_color(1.0,1.0,1.0);



uint gen_vertex_0(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u,v,w))); 
  float v1=static_cast<float>(mh->Value(img::Point(u,v+1,w)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u),
                                        static_cast<float>(v)+static_cast<float>((level-v0)/(v1-v0)),
                                        static_cast<float>(w))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_1(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u,v+1,w))); 
  float v1=static_cast<float>(mh->Value(img::Point(u+1,v+1,w)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u)+static_cast<float>((level-v0)/(v1-v0)),
                                        static_cast<float>(v+1),
                                        static_cast<float>(w))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_2(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u+1,v,w))); 
  float v1=static_cast<float>(mh->Value(img::Point(u+1,v+1,w)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u+1),
                                        static_cast<float>(v)+static_cast<float>((level-v0)/(v1-v0)),
                                        static_cast<float>(w))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_3(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u,v,w))); float v1=static_cast<float>(mh->Value(img::Point(u+1,v,w)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u)+static_cast<float>((level-v0)/(v1-v0)),
                                        static_cast<float>(v),
                                        static_cast<float>(w))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_4(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u,v,w))); float v1=static_cast<float>(mh->Value(img::Point(u,v,w+1)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u),
                                        static_cast<float>(v),
                                        static_cast<float>(w)+static_cast<float>((level-v0)/(v1-v0)))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_5(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u,v+1,w))); float v1=static_cast<float>(mh->Value(img::Point(u,v+1,w+1)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u),
                                        static_cast<float>(v+1),
                                        static_cast<float>(w)+static_cast<float>((level-v0)/(v1-v0)))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_6(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u+1,v+1,w))); float v1 = static_cast<float>(mh->Value(img::Point(u+1,v+1,w+1)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u+1),
                                        static_cast<float>(v+1),
                                        static_cast<float>(w)+static_cast<float>((level-v0)/(v1-v0)))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_7(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u+1,v,w))); float v1=static_cast<float>(mh->Value(img::Point(u+1,v,w+1)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u+1),
                                        static_cast<float>(v),
                                        static_cast<float>(w)+static_cast<float>((level-v0)/(v1-v0)))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_8(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u,v,w+1))); float v1=static_cast<float>(mh->Value(img::Point(u,v+1,w+1)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u),
                                        static_cast<float>(v)+static_cast<float>((level-v0)/(v1-v0)),
                                        static_cast<float>(w+1))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_9(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u,v+1,w+1))); float v1=static_cast<float>(mh->Value(img::Point(u+1,v+1,w+1)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u)+static_cast<float>((level-v0)/(v1-v0)),
                                        static_cast<float>(v+1),
                                        static_cast<float>(w+1))),
                  geom::Vec3(),def_color);
  }
  return 0;
}

uint gen_vertex_A(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u+1,v,w+1))); float v1=static_cast<float>(mh->Value(img::Point(u+1,v+1,w+1)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u+1),
                                        static_cast<float>(v)+static_cast<float>((level-v0)/(v1-v0)),
                                        static_cast<float>(w+1))),
                  geom::Vec3(),def_color);
  }
  return 0;
}
  
uint gen_vertex_B(img::RealSpatialImageState* mh, IndexedVertexArray& va, int u, int v, int w, float level)
{
  float v0=static_cast<float>(mh->Value(img::Point(u,v,w+1))); float v1=static_cast<float>(mh->Value(img::Point(u+1,v,w+1)));
  if((level>v0 && level<=v1) || (level>v1 && level<=v0)) {
    return va.Add(mh->FractionalIndexToCoord(geom::Vec3(static_cast<float>(u)+static_cast<float>((level-v0)/(v1-v0)),
                                        static_cast<float>(v),
                                        static_cast<float>(w+1))),
                  geom::Vec3(),def_color);
  }
  return 0;
}
  
}

}}} // ns
