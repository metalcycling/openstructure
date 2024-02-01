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
#ifndef OST_GFX_OBJ_BASE_HH
#define OST_GFX_OBJ_BASE_HH

/*
  Author: Ansgar Philippsen
*/

#include <ost/geom/geom.hh>
#include <ost/mol/entity_handle.hh>
#include <ost/mol/entity_view.hh>
#include <ost/img/map.hh>

#include <ost/config.hh>
#include "module_config.hh"

#include "gfx_node.hh"
#include "color.hh"
#include "render_mode.hh"
#include "gradient.hh"
#include "material.hh"
#include "color_ops/color_op.hh"

namespace ost { namespace gfx {

class Scene; // fw decl

/// \brief main interface for all graphic objects, both in C++ and Python
/*
  this abstract base class is for both the actual GfxObj as well as a 
  shortly implemented re-design of the GfxNode, and hence temporary in
  its current form
*/
class DLLEXPORT_OST_GFX GfxObjBase: public GfxNode
{
 public:
  GfxObjBase(const String& name): GfxNode(name) {}
  
  /// \brief material ambient color
  virtual void SetMatAmb(const Color& c) = 0;
  /// \brief material diffuse color
  virtual void SetMatDiff(const Color& c) = 0;
  /// \brief material specular color
  virtual void SetMatSpec(const Color& c) = 0;
  /// \brief material shininess
  virtual void SetMatShin(float s) = 0;
  /// \brief material emmissive color
  virtual void SetMatEmm(const Color& c) = 0;

  /// \brief rendering context switch, some things need to be regenerated
  virtual void ContextSwitch() = 0;

  /// \brief change render mode
  virtual void SetRenderMode(RenderMode::Type m) = 0;
  /// \brief current render mode
  virtual RenderMode::Type GetRenderMode() const = 0;
  
  /// \brief get geometric center
  virtual geom::Vec3 GetCenter() const = 0;

  /// \brief set line width
  virtual void SetLineWidth(float w) = 0;
  
  /// \brief set polygon mode
  virtual void SetPolyMode(unsigned int m) = 0;

  /// \brief turn on sophisticated line anti-aliasing, requires shader
  virtual void SetAALines(bool f) = 0;

  /// \brief turn on line halo of given strength
  virtual void SetLineHalo(float f) = 0;

  /// \brief turn outline rendering on or off
  virtual void SetOutline(bool f) = 0;
  /// \brief get state of outline rendering
  virtual bool GetOutline() const = 0;
  /// \brief set outline mode, 1, 2 or 3
  virtual void SetOutlineMode(int m) = 0;
  /// \brief get current outline mode
  virtual int GetOutlineMode() const = 0;
  /// \brief set outline width in pixels (modes 1 + 2)
  /// this does not scale with resolution
  virtual void SetOutlineWidth(float f) = 0;
  /// \brief get current outline width
  virtual float GetOutlineWidth() const = 0;
  /// \brief set outline expansion factor (mode 3)
  /// this scales with resolution
  virtual void SetOutlineExpandFactor(float f) = 0;
  /// \brief get current outline expand factor (mode 3)
  virtual float GetOutlineExpandFactor() const = 0;
  /// \brief set outline color (mode 3)
  virtual void SetOutlineExpandColor(const Color& c) = 0;
  /// \brief get current outline color (mode 3)
  virtual Color GetOutlineExpandColor() const = 0;

  /// \brief set opacity (1 = no transparency)
  virtual void SetOpacity(float f) = 0;
  /// \brief returns a value smaller than 1.0 if transparency is used in this object
  virtual float GetOpacity() const = 0;

  virtual void SetSolid(bool f) = 0;
  virtual bool GetSolid() const = 0;
  virtual void SetSolidColor(const Color& c) = 0;
  virtual Color GetSolidColor() const = 0;

  virtual void SetClip(bool f) = 0;
  virtual bool GetClip() const = 0;
  virtual void SetClipPlane(const geom::Vec4&) = 0;
  virtual geom::Vec4 GetClipPlane() const = 0;
  virtual void SetClipOffset(float f) = 0;
  virtual float GetClipOffset() const = 0;
  
  /// \brief color each component based on the gradient-mapped property of 
  ///    the given entity
  virtual void ColorBy(const mol::EntityView& ev, 
                       const String& prop,
                       const Gradient& g, float minv, float maxv) = 0;

  /// \brief color based on image
  virtual void ColorBy(const img::MapHandle& mh,
                       const String& prop,
                       const Gradient& g,float minv, float maxv) = 0;    
};

}} // ns

#endif
