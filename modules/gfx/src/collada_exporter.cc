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
#include <boost/format.hpp>

#include "collada_exporter.hh"
#include "scene.hh"

namespace ost { namespace gfx {

ColladaExporter::ColladaExporter(const std::string& file):
  Exporter(),
  file_(file),
  out_(file_.c_str()),
  obj_()
{
  out_ << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";
  out_ << "<COLLADA version=\"1.4.0\" xmlns=\"http://www.collada.org/2005/11/COLLADASchema\">\n";
}

ColladaExporter::~ColladaExporter()
{
  out_ << "</COLLADA>\n";
}

void ColladaExporter::SceneStart(const Scene* scene)
{
  out_ << " <library_cameras>\n";
	out_ << "  <camera id=\"Camera-Camera\" name=\"Camera-Camera\">\n";
  out_ << "   <optics>\n";
  out_ << "    <technique_common>\n";
  out_ << "      <perspective>\n";
  out_ << "       <yfov>" << scene->GetFOV() << "</yfov>\n";
  out_ << "       <znear>" << scene->GetNear() << "</znear>\n";
  out_ << "       <zfar>" << scene->GetFar() << "</zfar>\n";
  out_ << "      </perspective>\n";
  out_ << "    </technique_common>\n";
  out_ << "   </optics>\n";
	out_ << "  </camera>\n";
  out_ << " </library_cameras>\n";

  out_ << " <library_geometries>\n" << std::endl;
}

void ColladaExporter::SceneEnd(const Scene* scene)
{
  out_ << " </library_geometries>\n" << std::endl;

  out_ << " <library_visual_scenes>\n";
  out_ << "  <visual_scene id=\"Scene\" name=\"Scene\">\n";

  out_ << "   <node id=\"Camera\" name=\"Camera\">\n";
  out_ << "    <instance_camera url=\"#Camera-Camera\"/>\n";
  out_ << "   </node>\n";

  out_ << "   <node id=\"Root\" name=\"Root\">\n";
  if(GetToOrigin()) {
    out_ << "    <matrix>\n";
    geom::Vec3 cen=scene->GetTransform().GetCenter();
    out_ << "    1 0 0 " << -cen[0] << "\n";
    out_ << "    0 1 0 " << -cen[1] << "\n";
    out_ << "    0 0 1 " << -cen[2] << "\n";
    out_ << "    0 0 0 1\n";
    out_ << "    </matrix>\n";
  } else {
    out_ << "    <matrix>\n";
    geom::Mat4 tm=scene->GetTransform().GetMatrix();
    out_ << "    " << tm(0,0) << " " << tm(0,1) << " " << tm(0,2) << " " << tm(0,3) << "\n";
    out_ << "    " << tm(1,0) << " " << tm(1,1) << " " << tm(1,2) << " " << tm(1,3) << "\n";
    out_ << "    " << tm(2,0) << " " << tm(2,1) << " " << tm(2,2) << " " << tm(2,3) << "\n";
    out_ << "    " << tm(3,0) << " " << tm(3,1) << " " << tm(3,2) << " " << tm(3,3) << "\n";
    out_ << "    </matrix>\n";
  }
  for(std::vector<std::string>::const_iterator oit=obj_.begin();oit!=obj_.end();++oit) {
    out_ << "    <node id=\"" << *oit << "\" name=\"" << *oit <<"\">\n";
    out_ << "     <instance_geometry url=\"#" << *oit << "\"/>\n";
    out_ << "    </node>\n";
  }
  out_ << "   </node>\n";
  out_ << "  </visual_scene>\n";
	out_ << " </library_visual_scenes>\n";

	out_ << " <scene>\n";
  out_ << "  <instance_visual_scene url=\"#Scene\"/>\n";
	out_ << " </scene>\n";
}

void ColladaExporter::NodeStart(const std::string& name, NodeType t)
{
  if(name!="Scene") {
    obj_.push_back(name);
    out_ << "<geometry id=\"" << name << "\" name=\"" << name << "\">\n";
    out_ << " <mesh>\n";
  }
}

void ColladaExporter::NodeEnd(const std::string& name)
{
  if(name!="Scene") {
    out_ << " </mesh>\n";
    out_ << "</geometry>\n";
  }
}

void ColladaExporter::WriteVertexData(const float* vdata,
                                      const float* ndata, 
                                      const float* cdata,
                                      const float* tdata,
                                      size_t stride2, size_t count)
{
  std::string name=obj_.back();
  size_t stride=stride2/sizeof(float);

  out_ << "  <source id=\"" << name+"-Positions" << "\">\n";
  out_ << "   <float_array count=\"" << count*3 << "\" id=\"" << name+"-Positions-array" << "\">\n";
  if(vdata) {
    const float* src=vdata;
    float tmpv[3];
    for(unsigned int i=0;i<count;++i) {
      tmpv[0]=src[0]; tmpv[1]=src[1]; tmpv[2]=src[2];
      TransformPosition(tmpv);
      out_ << tmpv[0] << " ";
      out_ << tmpv[1] << " ";
      out_ << tmpv[2] << " ";
      src+=stride;
    }
  } else {
    for(unsigned int i=0;i<count;++i) {
      out_ << "0.0 0.0 0.0 ";
    }
  }
  out_ << "   </float_array>\n";
  out_ << "   <technique_common>\n";
  out_ << "    <accessor count=\"" << count << "\" source=\"#" << name + "-Positions-array" << "\" stride=\"3\">\n";
  out_ << "     <param name=\"X\" type=\"float\"/>\n";
  out_ << "     <param name=\"Y\" type=\"float\"/>\n";
  out_ << "     <param name=\"Z\" type=\"float\"/>\n";
  out_ << "    </accessor>\n";
  out_ << "   </technique_common>\n";
  out_ << "  </source>\n";

  // normals, lots of code duplication for now
  out_ << "  <source id=\"" << name+"-Normals" << "\">\n";
  out_ << "   <float_array count=\"" << count*3 << "\" id=\"" << name+"-Normals-array" << "\">\n";
  if(ndata) {
    const float* src=ndata;
    float tmpn[3];
    for(unsigned int i=0;i<count;++i) {
      tmpn[0]=src[0]; tmpn[1]=src[1]; tmpn[2]=src[2];
      TransformNormal(tmpn);
      out_ << tmpn[0] << " ";
      out_ << tmpn[1] << " ";
      out_ << tmpn[2] << " ";
      src+=stride;
    }
  } else {
    for(unsigned int i=0;i<count;++i) {
      out_ << "0.0 0.0 0.0 ";
    }
  }
  out_ << "   </float_array>\n";
  out_ << "   <technique_common>\n";
  out_ << "    <accessor count=\"" << count << "\" source=\"#" << name + "-Normals-array" << "\" stride=\"3\">\n";
  out_ << "     <param name=\"X\" type=\"float\"/>\n";
  out_ << "     <param name=\"Y\" type=\"float\"/>\n";
  out_ << "     <param name=\"Z\" type=\"float\"/>\n";
  out_ << "    </accessor>\n";
  out_ << "   </technique_common>\n";
  out_ << "  </source>\n";

  // colors, again lots of code duplication
  out_ << "  <source id=\"" << name+"-Colors" << "\">\n";
  //out_ << "   <float_array count=\"" << count*4 << "\" id=\"" << name+"-Colors-array" << "\">\n";
  out_ << "   <float_array count=\"" << count*3 << "\" id=\"" << name+"-Colors-array" << "\">\n";
  if(cdata) {
    const float* src=cdata;
    for(unsigned int i=0;i<count;++i) {
      out_ << src[0] << " ";
      out_ << src[1] << " ";
      out_ << src[2] << " ";
      //out_ << src[3] << " ";
      src+=stride;
    }
  } else {
    for(unsigned int i=0;i<count;++i) {
      out_ << "0.0 0.0 0.0 ";
    }
  }
  out_ << "   </float_array>\n";
  out_ << "   <technique_common>\n";
  //out_ << "    <accessor count=\"" << count << "\" source=\"#" << name + "-Colors-array" << "\" stride=\"4\">\n";
  out_ << "    <accessor count=\"" << count << "\" source=\"#" << name + "-Colors-array" << "\" stride=\"3\">\n";
  out_ << "     <param name=\"R\" type=\"float\"/>\n";
  out_ << "     <param name=\"G\" type=\"float\"/>\n";
  out_ << "     <param name=\"B\" type=\"float\"/>\n";
  //out_ << "     <param name=\"A\" type=\"float\"/>\n";
  out_ << "    </accessor>\n";
  out_ << "   </technique_common>\n";
  out_ << "  </source>\n";

  out_ << "  <vertices id=\"" << name+"-Vertex" << "\">\n";
  out_ << "   <input semantic=\"POSITION\" source =\"#" << name+"-Positions" << "\"/>\n";
  //out_ << "   <input semantic=\"NORMAL\" source =\"#" << name+"-Normals" << "\"/>\n";
  //out_ << "   <input semantic=\"COLOR\" source =\"#" << name+"-Colors" << "\"/>\n";
  out_ << "  </vertices>\n";
}

void ColladaExporter::WritePointData(const unsigned int* i, size_t count)
{
}

void ColladaExporter::WriteLineData(const unsigned int* ij, size_t count)
{
}

void ColladaExporter::WriteTriData(const unsigned int* ijk, size_t count)
{
  std::string name=obj_.back();

  out_ << "  <triangles count=\"" << count << "\">\n";
  out_ << "   <input offset=\"0\" semantic=\"VERTEX\" source =\"#" << name+"-Vertex" << "\"/>\n";
  out_ << "   <input offset=\"1\" semantic=\"NORMAL\" source =\"#" << name+"-Normals" << "\"/>\n";
  out_ << "   <input offset=\"2\" semantic=\"COLOR\" source =\"#" << name+"-Colors" << "\"/>\n";
  out_ << "   <p>\n";
  for(unsigned int c=0;c<count*3;c+=3) {
    out_ << ijk[c+0] << " ";
    out_ << ijk[c+0] << " ";
    out_ << ijk[c+0] << " ";
    out_ << ijk[c+1] << " ";
    out_ << ijk[c+1] << " ";
    out_ << ijk[c+1] << " ";
    out_ << ijk[c+2] << " ";
    out_ << ijk[c+2] << " ";
    out_ << ijk[c+2] << " ";
  }
  out_ << "   </p>\n";
  out_ << "  </triangles>\n";
}

void ColladaExporter::WriteQuadData(const unsigned int* ijkl, size_t count)
{
  std::string name=obj_.back();

  out_ << "  <triangles count=\"" << count << "\">\n";
  out_ << "   <input offset=\"0\" semantic=\"VERTEX\" source =\"#" << name+"-Vertex" << "\"/>\n";
  out_ << "   <input offset=\"1\" semantic=\"NORMAL\" source =\"#" << name+"-Normals" << "\"/>\n";
  out_ << "   <input offset=\"2\" semantic=\"COLOR\" source =\"#" << name+"-Colors" << "\"/>\n";
  out_ << "   <p>\n";
  for(unsigned int c=0;c<count*4;c+=4) {
    out_ << ijkl[c+0] << " ";
    out_ << ijkl[c+0] << " ";
    out_ << ijkl[c+0] << " ";
    out_ << ijkl[c+1] << " ";
    out_ << ijkl[c+1] << " ";
    out_ << ijkl[c+1] << " ";
    out_ << ijkl[c+2] << " ";
    out_ << ijkl[c+2] << " ";
    out_ << ijkl[c+2] << " ";
    //
    out_ << ijkl[c+0] << " ";
    out_ << ijkl[c+0] << " ";
    out_ << ijkl[c+0] << " ";
    out_ << ijkl[c+2] << " ";
    out_ << ijkl[c+2] << " ";
    out_ << ijkl[c+2] << " ";
    out_ << ijkl[c+3] << " ";
    out_ << ijkl[c+3] << " ";
    out_ << ijkl[c+3] << " ";
  }
  out_ << "   </p>\n";
  out_ << "  </triangles>\n";
}

}} // ns
