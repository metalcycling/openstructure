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
#ifndef OST_GFX_EXPORTER_HH
#define OST_GFX_EXPORTER_HH

#include <ost/geom/vec3.hh>
#include <ost/geom/vec4.hh>
#include <ost/geom/mat3.hh>
#include <ost/geom/mat4.hh>

#include <ost/gfx/module_config.hh>

namespace ost { namespace gfx {

class Scene;

class DLLEXPORT_OST_GFX Exporter
{
public:
  enum NodeType {
    ROOT=1,
    GROUP=2,
    OBJ=3
  };

  Exporter() :
    scale_(1.0),
    to_origin_(true)
  {}
  virtual ~Exporter() {}
  virtual void SceneStart(const Scene* scene) {}
  virtual void SceneEnd(const Scene* scene) {}

  virtual void NodeStart(const std::string& name, NodeType t) {}
  virtual void NodeEnd(const std::string& name) {}
  
  // this indicates beginning of new data, including a reset of the indices
  // may occur more than once for a given node
  virtual void WriteVertexData(const float* v, const float* n, const float* c, const float* t,
                               size_t stride, size_t count) {}
  virtual void WritePointData(const unsigned int* i, size_t count) {}
  virtual void WriteLineData(const unsigned int* ij, size_t count) {}
  virtual void WriteTriData(const unsigned int* ijk, size_t count) {}
  virtual void WriteQuadData(const unsigned int* ijkl, size_t count) {}

  // scale positions for absolute data formats (like dae)
  void SetScale(float s) {scale_=s;}
  float GetScale() const {return scale_;}
  // if true (default), re-orient so that center is at (0,0,0)
  // and viewing direction is along z
  // (basically apply modelview matrix)
  void SetToOrigin(bool b) {to_origin_=b;}
  bool GetToOrigin() const {return to_origin_;}

  // used by Scene::Export
  void SetupTransform(const Scene* scene);
  // used by WriteVertexData in derived classes
  // modifies input arg!!
  void TransformPosition(float* p) const;
  void TransformNormal(float* n) const;

private:
  float scale_;
  bool to_origin_;
  geom::Mat4 vertex_tf_;
  geom::Mat3 normal_tf_;
};



}} // ns

#endif
