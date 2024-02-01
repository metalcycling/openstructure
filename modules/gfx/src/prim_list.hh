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
#ifndef OST_GFX_PRIM_LIST_HH
#define OST_GFX_PRIM_LIST_HH

/*
  Author: Ansgar Philippsen
*/

#include <boost/shared_ptr.hpp>

#include <ost/geom/geom.hh>

#include "gfx_object.hh"
#include "gfx_prim.hh"

namespace ost { namespace gfx {

/// \see gfx::Cuboid, \ref primitives.py "Displaying Lines and Quads",
///      gfx::Quad, \ref gradient.py "Gradient Example",
///      \ref random_lines.py "Random Lines"
class DLLEXPORT_OST_GFX PrimList: public GfxObj 
{
  struct PointEntry {
    PointEntry(const geom::Vec3& p, float r, const Color& c):
      pos(p), rad(r), col(c) {}
    geom::Vec3 pos;
    float rad;
    Color col;
  };

  typedef std::vector<PointEntry> PointEntryList;

  struct LineEntry {
    LineEntry(const geom::Vec3& p1, const geom::Vec3& p2, float r1, float r2, const Color& c1, const Color& c2):
      pos1(p1), pos2(p2), radius_1(r1), radius_2(r2), col1(c1), col2(c2) {}
    geom::Vec3 pos1, pos2;
    float radius_1, radius_2;
    Color col1, col2;
  };

  typedef std::vector<LineEntry> LineEntryList;

 public:
  /// \brief create new prim list
  PrimList(const String& name);

  virtual geom::AlignedCuboid GetBoundingBox(bool use_tf=true) const;

  virtual void ProcessLimits(geom::Vec3& minc, geom::Vec3& maxc, 
                             const geom::Transform& tf) const;
  /// \brief get center      
  virtual geom::Vec3 GetCenter() const;

  virtual void CustomRenderPov(PovState& pov);

  virtual void CustomRenderGL(RenderPass pass);

  virtual void OnRenderModeChange();

  virtual void SetLineWidth(float w);
  /// \brief clear all prims
  void Clear();

  /*!
     \brief add point

     python interface:

     PrimList.AddPoint(pos, color=gfx.Color())
  */
  void AddPoint(const geom::Vec3& p, const Color& col);

  /*!
    \brief add line

    Python interface:

    PrimList.AddLine(pos1,pos2,color=gfx.WHITE,color1=gfx.WHITE,color2=gfx.WHITE)
  */
  void AddLine(const geom::Vec3& p1, const geom::Vec3& p2, const Color& col1, const Color& col2);

  /*!
    \brief add sphere

    Python interface:

    PrimList.AddSphere(cen,radius=1.0,color=gfx.WHITE)
  */
  void AddSphere(const geom::Vec3& cen, float rad, const Color& col);

  /*!
    \brief add cylinder

    Python interface:

    PrimList.AddCyl(pos1,pos2,radius=1.0,radius1=1.0,radius2=1.0,color=gfx.WHITE,color1=gfx.WHITE,color2=gfx.WHITE)
  */
  void AddCyl(const geom::Vec3& p0, const geom::Vec3& p1, float r1, float r2, const Color& col1, const Color& col2);

  /*!
    \brief add text

    Python interface:

    PrimList.AddText(text,pos,color=gfx.WHITE,point_size=1.0)
  */
  void AddText(const std::string& text, const geom::Vec3& pos, const Color& col, float point_size);

  /// defunct
  void SetDiameter(float d);
  
  /// defunct
  void SetRadius(float r);

  /// \brief set global prims color, overriding individual ones
  void SetColor(const Color& c);

  void SetSphereDetail(unsigned int d);
  unsigned int GetSphereDetail() const {return sphere_detail_;}

  void SetArcDetail(unsigned int d);
  unsigned int GetArcDetail() const {return arc_detail_;}

  // TODO: add point and line pixel width

  /*!
    \brief add triangle mesh

    v : pointer to nv*3 floats for the positions (mandatory)
    n : pointer to nv*3 floats for the normals (may be NULL)
    c : pointer to nv*4 floats for the colors (may be NULL)
    nv: number of vertices, normals, and colors
    i : pointer to ni*3 vertex indices
    ni: number of index triplets

    Python interface, using numpy arrays:

      AddMesh(vertex_array, normal_array, color_array, index_array)

    where normal_array and color_array may be None
  */
  void AddMesh(float* v, float* n, float* c, size_t nv, unsigned int* i, size_t ni);

 protected:
  virtual void CustomPreRenderGL(bool flag);

 private:
  SpherePrimList points_;
  CylinderPrimList lines_;
  SpherePrimList spheres_;
  CylinderPrimList cyls_;
  TextPrimList texts_;
  unsigned int sphere_detail_;
  unsigned int arc_detail_;

  IndexedVertexArray simple_va_;

  std::vector<IndexedVertexArray> vas_;
  
  void prep_simple_va();
  void prep_va();
  void render_text();
};

/// \example primitives.py
/// 
/// Shows how to create grapical objects to render lines and rectangles in 3D 
/// space.
/// \sa \ref bbox.py "Bounding Box Example", \ref random_lines.py "Random Lines"

/// \example random_lines.py
/// 
///  Draw lines with random orientation in 3D space
/// \sa \ref bbox.py "Bounding Box Example", \ref primitives.py "Lines And Quads" 
}} // ns

#endif
