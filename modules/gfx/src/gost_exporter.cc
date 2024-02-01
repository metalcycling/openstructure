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

#include <sstream>

#include "gost_exporter.hh"
#include "scene.hh"

namespace ost { namespace gfx {

  namespace {
    enum TYPE {
      GOST_UNKNOWN=0,
      GOST_SCENE=1,
      GOST_START_NODE=2,
      GOST_END_NODE=3,
      GOST_DATA=4,
      GOST_FRAME=5
    };

    enum NODE_TYPE {
      GOST_ROOT=1,
      GOST_GROUP=2,
      GOST_OBJ=3
    };

    enum DATA_TYPE {
      GOST_VDATA=1,
      GOST_NDATA=2,
      GOST_CDATA=3,
      GOST_TDATA=4,
      GOST_PINDEX=5,
      GOST_LINDEX=6,
      GOST_TINDEX=7,
      GOST_QINDEX=8
    };
  }

GostExporter::GostExporter(const std::string& fname):
  Exporter(),
  file_(0)
{
  file_=fopen(fname.c_str(),"w");
  if(!file_) {
    std::ostringstream m;
    m << "Could not open '" << fname << "' for writing";
    throw Error(m.str());
  }
  char headr[]={'G','O','S','T','1','0'};
  fwrite(headr,sizeof(char),6,file_);
}

GostExporter::~GostExporter()
{
  fclose(file_);
  file_=0;
}

void GostExporter::SceneStart(const Scene* scene)
{
  int type=GOST_SCENE;
  int subtype=0;
  size_t size=0;

  fwrite(&type,sizeof(int),1,file_);
  fwrite(&subtype,sizeof(int),1,file_);
  fwrite(&size,sizeof(size_t),1,file_);
}

void GostExporter::SceneEnd(const Scene* scene)
{}

void GostExporter::NodeStart(const std::string& name, NodeType t)
{
  int type=GOST_START_NODE;
  int subtype=0;
  if(t==Exporter::ROOT) subtype=GOST_ROOT;
  else if(t==Exporter::GROUP) subtype=GOST_GROUP;
  else if(t==Exporter::OBJ) subtype=GOST_OBJ;
  size_t size=name.size();

  fwrite(&type,sizeof(int),1,file_);
  fwrite(&subtype,sizeof(int),1,file_);
  fwrite(&size,sizeof(size_t),1,file_);
  fwrite(name.c_str(),sizeof(char), size,file_);
}

void GostExporter::NodeEnd(const std::string& name)
{
  int type=GOST_END_NODE;
  int subtype=0;
  size_t size=0;

  fwrite(&type,sizeof(int),1,file_);
  fwrite(&subtype,sizeof(int),1,file_);
  fwrite(&size,sizeof(size_t),1,file_);
}

void GostExporter::WriteVertexData(const float* vdata, const float* ndata, 
                                   const float* cdata, const float* tdata,
                                   size_t stride2, size_t count)
{
  std::vector<float> buffer(count*4);
  int type=GOST_DATA;
  int subtype=0;
  size_t size=0;
  size_t stride=stride2/sizeof(float);

  if(vdata) {
    // positions
    subtype=GOST_VDATA;
    size=sizeof(float)*3*count+sizeof(size_t);
    float* dest=&buffer[0];
    const float* src=vdata;
    for(size_t i=0;i<count;++i) {
      dest[0]=src[0];
      dest[1]=src[1];
      dest[2]=src[2];
      src+=stride;
      dest+=3;
    }
    fwrite(&type,sizeof(int),1,file_);
    fwrite(&subtype,sizeof(int),1,file_);
    fwrite(&size,sizeof(size_t),1,file_);
    //
    fwrite(&count,sizeof(size_t),1,file_);
    fwrite(&buffer[0],sizeof(float),3*count,file_);
  }

  if(ndata) {
    // normals
    subtype=GOST_NDATA;
    size=sizeof(float)*count*3+sizeof(size_t);
    float* dest=&buffer[0];
    const float* src=ndata;
    for(size_t i=0;i<count;++i) {
      dest[0]=src[0];
      dest[1]=src[1];
      dest[2]=src[2];
      src+=stride;
      dest+=3;
    }
    fwrite(&type,sizeof(int),1,file_);
    fwrite(&subtype,sizeof(int),1,file_);
    fwrite(&size,sizeof(size_t),1,file_);
    //
    fwrite(&count,sizeof(size_t),1,file_);
    fwrite(&buffer[0],sizeof(float),3*count,file_);
  }
  
  if(cdata) {
    // colors
    subtype=GOST_CDATA;
    size=sizeof(float)*count*4+sizeof(size_t);
    float* dest=&buffer[0];
    const float* src=cdata;
    for(size_t i=0;i<count;++i) {
      dest[0]=src[0];
      dest[1]=src[1];
      dest[2]=src[2];
      dest[3]=src[3];
      src+=stride;
      dest+=4;
    }
    fwrite(&type,sizeof(int),1,file_);
    fwrite(&subtype,sizeof(int),1,file_);
    fwrite(&size,sizeof(size_t),1,file_);
    //
    fwrite(&count,sizeof(size_t),1,file_);
    fwrite(&buffer[0],sizeof(float),4*count,file_);
  }

  if(tdata) {
    // texture coordinates
    subtype=GOST_TDATA;
    size=sizeof(float)*count*2+sizeof(size_t);
    float* dest=&buffer[0];
    const float* src=tdata;
    for(size_t i=0;i<count;++i) {
      dest[0]=src[0];
      dest[1]=src[1];
      src+=stride;
      dest+=2;
    }
    fwrite(&type,sizeof(int),1,file_);
    fwrite(&subtype,sizeof(int),1,file_);
    fwrite(&size,sizeof(size_t),1,file_);
    //
    fwrite(&count,sizeof(size_t),1,file_);
    fwrite(&buffer[0],sizeof(float),2*count,file_);
  }
}

void GostExporter::WritePointData(const unsigned int* i, size_t count)
{
  int type=GOST_DATA;
  int subtype=GOST_PINDEX;
  size_t size=sizeof(unsigned int)*count+sizeof(size_t);
  fwrite(&type,sizeof(int),1,file_);
  fwrite(&subtype,sizeof(int),1,file_);
  fwrite(&size,sizeof(size_t),1,file_);
  //
  fwrite(&count,sizeof(size_t),1,file_);
  fwrite(i,sizeof(unsigned int),count,file_);
}

void GostExporter::WriteLineData(const unsigned int* ij, size_t count)
{
  int type=GOST_DATA;
  int subtype=GOST_LINDEX;
  size_t size=sizeof(unsigned int)*2*count+sizeof(size_t);
  fwrite(&type,sizeof(int),1,file_);
  fwrite(&subtype,sizeof(int),1,file_);
  fwrite(&size,sizeof(size_t),1,file_);
  //
  fwrite(&count,sizeof(size_t),1,file_);
  fwrite(ij,sizeof(unsigned int),2*count,file_);
}

void GostExporter::WriteTriData(const unsigned int* ijk, size_t count)
{
  int type=GOST_DATA;
  int subtype=GOST_TINDEX;
  size_t size=sizeof(unsigned int)*3*count+sizeof(size_t);
  fwrite(&type,sizeof(int),1,file_);
  fwrite(&subtype,sizeof(int),1,file_);
  fwrite(&size,sizeof(size_t),1,file_);
  //
  fwrite(&count,sizeof(size_t),1,file_);
  fwrite(ijk,sizeof(unsigned int),3*count,file_);
}

void GostExporter::WriteQuadData(const unsigned int* ijkl, size_t count)
{
  int type=GOST_DATA;
  int subtype=GOST_QINDEX;
  size_t size=sizeof(unsigned int)*4*count+sizeof(size_t);
  fwrite(&type,sizeof(int),1,file_);
  fwrite(&subtype,sizeof(int),1,file_);
  fwrite(&size,sizeof(size_t),1,file_);
  //
  fwrite(&count,sizeof(size_t),1,file_);
  fwrite(ijkl,sizeof(unsigned int),4*count,file_);
}

void GostExporter::SetFrame(size_t frame)
{
  int type=GOST_FRAME;
  int subtype=0;
  size_t size=sizeof(size_t);
  fwrite(&type,sizeof(int),1,file_);
  fwrite(&subtype,sizeof(int),1,file_);
  fwrite(&size,sizeof(size_t),1,file_);
  //
  fwrite(&frame,sizeof(size_t),1,file_);
}

}} // ns
