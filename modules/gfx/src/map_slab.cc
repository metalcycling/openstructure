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
#include <vector>

#include <ost/log.hh>
#include <ost/profile.hh>
#include <ost/img/alg/stat.hh>

#include "gl_helper.hh"

#include "scene.hh"
#include "map_slab.hh"
#include "shader.hh"

namespace {

  float InterpolatedValue(const ost::img::ImageHandle& mh, const geom::Vec3& xyz)
  {
    geom::Vec3 uvw=mh.CoordToIndex(xyz);
    float x=fmod(float(uvw[0]), float(1.0));
    float y=fmod(float(uvw[1]), float(1.0));
    float z=fmod(float(uvw[2]), float(1.0));
    uint u=static_cast<uint>(uvw[0]);
    uint v=static_cast<uint>(uvw[1]);
    uint w=static_cast<uint>(uvw[2]);

    float value=0.0;
    if(u<mh.GetSize()[0]-1 && v<mh.GetSize()[1]-1 && w<mh.GetSize()[2]-1) {
     value=mh.GetReal(ost::img::Point(u,v,w))*(1.0-x)*(1.0-y)*(1.0-z)+
       mh.GetReal(ost::img::Point(u+1,v,w))*(x)*(1.0-y)*(1.0-z)+
       mh.GetReal(ost::img::Point(u,v+1,w))*(1.0-x)*(y)*(1.0-z)+
       mh.GetReal(ost::img::Point(u,v,w+1))*(1.0-x)*(1.0-y)*(z)+
       mh.GetReal(ost::img::Point(u+1,v,w+1))*(x)*(1.0-y)*(z)+
       mh.GetReal(ost::img::Point(u,v+1,w+1))*(1.0-x)*(y)*(z)+
       mh.GetReal(ost::img::Point(u+1,v+1,w))*(x)*(y)*(1.0-z)+
       mh.GetReal(ost::img::Point(u+1,v+1,w+1))*(x)*(y)*(z);
    }
    return static_cast<float>(value);
  }
} // ns

namespace ost { 

namespace gfx {

MapSlab::MapSlab(const String& name, const img::MapHandle& mh, const geom::Plane& p):
  GfxObj(name),
  mh_(mh),
  plane_(p),
  usize_(256),
  vsize_(256),
  tex_id_(0),
  tex_data_(new float[usize_*vsize_*3]),
  rebuild_(false),
  valid_(false),
  c0_(), c1_(), c2_(), c3_(),
  grad_()  
{
  grad_.SetColorAt(0.0, Color(0,0,0));
  grad_.SetColorAt(1.0, Color(1,1,1));
  Rebuild();
  img::alg::Stat stat;
  mh_.Apply(stat);
  minv_ = stat.GetMinimum();
  maxv_ = stat.GetMaximum();
}

geom::AlignedCuboid MapSlab::GetBoundingBox(bool use_tf) const
{
  geom::Vec3 minc = plane_.At(-1.0,-1.0);
  geom::Vec3 maxc = plane_.At(1.0,1.0);
  geom::AlignedCuboid bb(minc,maxc);
  return use_tf ? transform_.Apply(bb) : bb;
}

geom::Vec3 MapSlab::GetCenter() const
{
  geom::Vec3 nrvo = plane_.At(0.0,0.0);
  return nrvo;
}
  
void MapSlab::CustomPreRenderGL(bool flag)
{
  if(tex_id_==0) {
    glGenTextures(1,&tex_id_);
  }
}

void MapSlab::CustomRenderGL(RenderPass pass)
{
  if(rebuild_) rebuild_tex();

  if(pass==STANDARD_RENDER_PASS) {
    render_slab();
  }
}

void MapSlab::CustomRenderPov(PovState& pov)
{
}


void MapSlab::Rebuild()
{
  valid_=calc_poly();
  if(valid_)  rebuild_=true;
}

geom::Plane MapSlab::GetPlane() const
{
  return plane_;
}

void MapSlab::SetPlane(const geom::Plane& p)
{
  plane_=p;
  Rebuild();
  Scene::Instance().RequestRedraw();
}

void MapSlab::ColorBy(const Gradient& g, float minv, float maxv)
{
  grad_=g;
  minv_=minv;
  maxv_=maxv;
  rebuild_=true;
  Scene::Instance().RequestRedraw();
}

void MapSlab::ColorBy(const Gradient& g)
{
  grad_=g;
  ost::img::alg::Stat stat;
  mh_.Apply(stat);
  minv_=stat.GetMinimum();
  maxv_=stat.GetMaximum();
}

void MapSlab::ColorBy(const Color& c1, const Color& c2, float minv, float maxv)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  ColorBy(g,minv,maxv);
}

void MapSlab::ColorBy(const Color& c1, const Color& c2)
{
  Gradient g;
  g.SetColorAt(0.0,c1);
  g.SetColorAt(1.0,c2);
  ColorBy(g);
}


/////////////////////////////
// private methods

namespace {

std::pair<float, float> calc_sr(float ax, float ay, 
                                  float bx, float by,
                                  float cx, float cy)
{
  float det = ax*by-ay*bx;
  if(std::abs(det)<1e-10f) {
    return std::make_pair(0.0f,0.0f);
  }
  return std::make_pair((cx*by-cy*bx)/det,
                        (ax*cy-ay*cx)/det);
}

struct PEntry {
  PEntry() {}
  PEntry(const geom::Vec3& pp, uint aa): p(pp), a(aa) {}
  geom::Vec3 p;
  uint a;
};

typedef std::vector<PEntry> PList;

typedef std::vector<std::pair<PEntry,PEntry> > LList;

} // anon ns

/*
  this routine does all the slab work, and could
  probably be written in a more compact form
*/
bool MapSlab::calc_poly()
{
  static std::vector<geom::Vec3> corner;
  if(corner.empty()) {
    corner.push_back(geom::Vec3(0.0,0.0,0.0));
    corner.push_back(geom::Vec3(0.0,1.0,0.0));
    corner.push_back(geom::Vec3(1.0,1.0,0.0));
    corner.push_back(geom::Vec3(1.0,0.0,0.0));
    corner.push_back(geom::Vec3(0.0,0.0,1.0));
    corner.push_back(geom::Vec3(0.0,1.0,1.0));
    corner.push_back(geom::Vec3(1.0,1.0,1.0));
    corner.push_back(geom::Vec3(1.0,0.0,1.0));
  }      
  static uint edge[][2]={{0,1},{1,2},{2,3},{3,0},
                         {0,4},{1,5},{2,6},{3,7},
                         {4,5},{5,6},{6,7},{7,4}};
  static uint face[][2]={{0,1},{0,4},{0,3},{0,5},
                         {5,1},{1,4},{4,3},{3,5},
                         {1,2},{4,2},{3,2},{5,2}};
  
  geom::Vec3 ori = plane_.GetOrigin();
  geom::Vec3 norm = plane_.GetNormal();
  geom::Vec3 size(mh_.GetSize()[0],mh_.GetSize()[1],mh_.GetSize()[2]);

  // gather the edge intersection points
  PList plist;
  for(uint ec=0;ec<12;++ec) {
    geom::Vec3 p0 = mh_.FractionalIndexToCoord(geom::CompMultiply(size,corner[edge[ec][0]]));
    geom::Vec3 p1 = mh_.FractionalIndexToCoord(geom::CompMultiply(size,corner[edge[ec][1]]));
    float den = geom::Dot(p1-p0,norm);
    if(std::abs(den)>1e-10) {
      float r = geom::Dot(ori-p0,norm)/den;
      if(r>=0.0 && r<=1.0) {
        plist.push_back(PEntry(p0+r*(p1-p0),ec));
      }
    }
  }

  if(plist.empty()) {
    return false;
  }

  // gather the lines
  LList llist2;
  for(uint i=0;i<plist.size()-1;++i) {
    for(uint k=i+1;k<plist.size();++k) {
      if(face[plist[i].a][0]==face[plist[k].a][0] ||
         face[plist[i].a][1]==face[plist[k].a][0] ||
         face[plist[i].a][0]==face[plist[k].a][1] ||
         face[plist[i].a][1]==face[plist[k].a][1]) {
        llist2.push_back(std::make_pair(plist[i],plist[k]));
      }
    }
  }

  // and sort them
  LList llist;
  llist.push_back(llist2[0]);
  llist2.erase(llist2.begin());
  uint maxcount=llist2.size();
  uint count=0;
  while(!llist2.empty() && count<maxcount) {
    ++count;
    for(LList::iterator lit = llist2.begin();lit!=llist2.end();++lit) {
      bool ins=false;
      for(LList::iterator nit = llist.begin();nit!=llist.end();++nit) {
        if(lit->second.a==nit->first.a) {
          llist.insert(nit,*lit);
          llist2.erase(lit);
          ins=true;
          break;
        } else if(lit->first.a==nit->first.a) {
          llist.insert(nit,std::make_pair(lit->second,lit->first));
          llist2.erase(lit);
          ins=true;
          break;
        } else if(lit->first.a==nit->second.a) {
          ++nit;
          llist.insert(nit,*lit);
          llist2.erase(lit);
          ins=true;
          break;
        } else if(lit->first.a==nit->first.a) {
          ++nit;
          llist.insert(nit,std::make_pair(lit->second,lit->first));
          llist2.erase(lit);
          ins=true;
          break;
        }
      }
      if(ins) break;
    }
  }

  // find rectangle enclosing all points
  geom::Vec3 F = geom::Normalize(plist[0].p-ori);
  geom::Vec3 G = geom::Normalize(geom::Cross(norm,F));

  uint ui,vi;
  if(norm[2]>norm[0] && norm[2]>norm[1]) {
    ui=0; vi=1;
  } else if(norm[1]>norm[0] && norm[1]>norm[2]) {
    ui=0; vi=2;
  } else {
    ui=1; vi=2;
  }

  geom::Vec3 H = plist[0].p-ori;
  std::pair<float,float> min_sr = calc_sr(F[ui],F[vi],G[ui],G[vi],H[ui],H[vi]);
  std::pair<float,float> max_sr = min_sr;

  for(uint pc=1;pc<plist.size();++pc) {
    H = plist[pc].p-ori;
    std::pair<float,float> sr = calc_sr(F[ui],F[vi],G[ui],G[vi],H[ui],H[vi]);
    min_sr.first = std::min(min_sr.first,sr.first);
    min_sr.second = std::min(min_sr.second,sr.second);
    max_sr.first = std::max(max_sr.first,sr.first);
    max_sr.second = std::max(max_sr.second,sr.second);
  }

  // these are the enclosing rectangle limits
  c0_ = ori+min_sr.first*F+min_sr.second*G;
  c1_ = ori+max_sr.first*F+min_sr.second*G;
  c2_ = ori+max_sr.first*F+max_sr.second*G;
  c3_ = ori+min_sr.first*F+max_sr.second*G;

  // now assemble to poly list from the line-list and the rect limits
  poly_list_.clear();
  F = (c1_-c0_);
  G = (c3_-c0_);
  for(LList::const_iterator lit=llist.begin();lit!=llist.end();++lit) {
    geom::Vec3 p = lit->first.p;
    H = p-c0_;
    std::pair<float,float> sr = calc_sr(F[ui],F[vi],G[ui],G[vi],H[ui],H[vi]);
    geom::Vec2 t(sr.first,sr.second);
    poly_list_.push_back(PolyListEntry(p,t));
  }

  return true;
}

void MapSlab::rebuild_tex()
{
  rebuild_=false;
  if(valid_) {
    float iu=1.0/static_cast<float>(usize_-1);
    float iv=1.0/static_cast<float>(vsize_-1);
    geom::Vec3 d0=c1_-c0_;
    geom::Vec3 d1=c3_-c0_;
    for(uint vc=0;vc<vsize_;++vc) {
      geom::Vec3 p0=c0_+static_cast<float>(vc)*iv*d1;
      for(uint uc=0;uc<usize_;++uc) {
        geom::Vec3 p1 = static_cast<float>(uc)*iu*d0+p0;
        float val = InterpolatedValue(mh_,p1);
        uint off = 3*(vc*usize_+uc);
        Color col = grad_.GetColorAt((val-minv_)/(maxv_-minv_));
        tex_data_[off+0]=col.Red();
        tex_data_[off+1]=col.Green();
        tex_data_[off+2]=col.Blue();
      }
    }
    glBindTexture(GL_TEXTURE_2D, tex_id_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);   
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, usize_,vsize_,0,
                 GL_RGB, GL_FLOAT, tex_data_.get());
  }
}

void MapSlab::render_slab()
{
  if(valid_) {
#if OST_SHADER_SUPPORT_ENABLED
    Shader::Instance().PushProgram();
    Shader::Instance().Activate("");
#endif
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, tex_id_);
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glBegin(GL_POLYGON);
    glColor3f(1.0,1.0,1.0);
    glNormal3v(plane_.GetNormal().Data());
    for(PolyList::iterator pit=poly_list_.begin();pit!=poly_list_.end();++pit) {
      glTexCoord2v(pit->t.Data());
      glVertex3v(pit->p.Data());
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
    glBegin(GL_LINE_LOOP);
    glColor3f(1.0,1.0,1.0);
    glNormal3v(plane_.GetNormal().Data());
    for(PolyList::iterator pit=poly_list_.begin();pit!=poly_list_.end();++pit) {
      glVertex3v(pit->p.Data());
    }
    glEnd();

    glBegin(GL_LINES);
    glColor3f(1.0,1.0,1.0);
    glVertex3v(plane_.GetOrigin().Data());
    glVertex3v((plane_.GetOrigin()+plane_.GetNormal()).Data());
    glEnd();

    glPopAttrib();
#if OST_SHADER_SUPPORT_ENABLED
    Shader::Instance().PopProgram();
#endif
  }
}

}} // ns
