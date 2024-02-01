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
#ifndef OST_GFX_PRIMITIVES_HH
#define OST_GFX_PRIMITIVES_HH

#include <ost/geom/geom.hh>
#include <ost/gfx/module_config.hh>
#include <ost/gfx/color.hh>
#include <ost/gfx/gfx_node.hh>

/*
  Author: Marco Biasini
 */
namespace ost { namespace gfx {

/// \brief base class for geometric primitives such as cuboids, and circles
///
/// The geometric primitives support rendering in both outline and solid mode.
/// These properties are controllable with SetFill() and SetOutline() as well
/// as SetFillColor() and SetOutlineColor(). The line width can be adjusted
/// with SetLineWidth().
///
/// \sa gfx::Cuboid, gfx::Quad
/// \todo transparent primitives would require the gfx objects to be 
///     ordered back to front. This has to be properly implemented in the
///     rendering system.
class DLLEXPORT_OST_GFX Primitive : 
public GfxNode {
protected:
  Primitive(const String& name);
public:
  void SetFill(const Color& color, bool fill=true);
  void SetOutline(const Color& color, bool line=true, float width=1.0);
  void SetFillColor(const Color& color);
  void SetOutlineColor(const Color& color);
  const Color& GetFillColor() const;
  const Color& GetOutlineColor() const;
  bool HasFill() const;
  bool HasOutline() const;
  float GetLineWidth() const;
  
  void SetLineWidth(float width);
private:
  bool  fill_;
  bool  line_;
  Color outline_color_;
  Color fill_color_;
  float line_width_;
};
/// \brief a rendered cuboid in 3D space

class DLLEXPORT_OST_GFX Cuboid : public Primitive {
public:
  Cuboid(const String& name, const geom::Cuboid& cuboid);
  virtual void RenderGL(RenderPass pass);  
private:
  geom::Cuboid cuboid_;
};

/// \brief a simple rendered quad in 3D space. 
/// \todo this class should be part of a more sophisticated polygon rendering
///    class.
class DLLEXPORT_OST_GFX Quad : public Primitive {
public:
  Quad(const String& name, const geom::Vec3& a, const geom::Vec3& b,
       const geom::Vec3& c, const geom::Vec3& d);
  virtual void RenderGL(RenderPass pass);
private:
  geom::Vec3 corner_points_ [4];
};

/// \example primitives.py
///
/// Shows how to display quads and lines in the graphical scene

}}

#endif
