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
#ifndef OST_GFX_OBJ_HH
#define OST_GFX_OBJ_HH

/*
  Author: Ansgar Philippsen
*/
// I have really no clue who defines this macro, but apparently  its there.
#ifdef check
#undef check
#endif
#ifndef Q_MOC_RUN
#include <boost/ptr_container/ptr_vector.hpp>

#include <vector>

#include <ost/config.hh>
#include <ost/gfx/module_config.hh>

#include <ost/geom/transform.hh>

#include "gfx_object_fw.hh"
#include "gfx_object_base.hh"
#include "gfx_prim.hh"
#include "vertex_array.hh"
#include "input.hh"
#include "exporter_fw.hh"
#endif
namespace ost { namespace gfx {

class Scene; // fw decl

/// \brief main class for all graphic objects
class DLLEXPORT_OST_GFX GfxObj: public GfxObjBase
{
public:
  GfxObj(const String& name);

  // gfx node interface
  virtual GfxNodeP Copy() const;
  virtual void DeepSwap(GfxObj& go);
  virtual void RenderGL(RenderPass pass);
  virtual void RenderPov(PovState& pov);
  virtual void Export(Exporter* ex);
  virtual void Apply(GfxNodeVisitor& v, GfxNodeVisitor::Stack st);
  virtual int GetType() const;
  //

  // gfx obj base interface
  virtual void SetMatAmb(const Color& c);
  virtual void SetMatDiff(const Color& c);
  virtual void SetMatSpec(const Color& c);
  virtual void SetMatShin(float s);
  virtual void SetMatEmm(const Color& c);
  virtual void ContextSwitch();
  virtual void SetRenderMode(RenderMode::Type m);
  virtual RenderMode::Type GetRenderMode() const;
  virtual geom::Vec3 GetCenter() const;
  virtual void SetLineWidth(float w);
  virtual void SetPolyMode(unsigned int m);
  virtual void SetAALines(bool f);
  virtual void SetLineHalo(float f);
  virtual void SetOutline(bool f);
  virtual bool GetOutline() const {return outline_flag_;};
  virtual void SetOutlineMode(int m);
  virtual int GetOutlineMode() const {return outline_mode_;}
  virtual void SetOutlineWidth(float f);
  virtual float GetOutlineWidth() const;
  virtual void SetOutlineExpandFactor(float f);
  virtual float GetOutlineExpandFactor() const;
  virtual void SetOutlineExpandColor(const Color& c);
  virtual Color GetOutlineExpandColor() const;
  virtual void SetOpacity(float f);
  virtual float GetOpacity() const {return opacity_;}
  virtual void SetSolid(bool f);
  virtual bool GetSolid() const {return solid_;}
  virtual void SetSolidColor(const Color& c);
  virtual Color GetSolidColor() const {return solid_color_;}
  
  virtual void SetClip(bool f);
  virtual bool GetClip() const {return clip_flag_;}
  virtual void SetClipPlane(const geom::Vec4&);
  virtual geom::Vec4 GetClipPlane() const {return clip_plane_;}
  virtual void SetClipOffset(float f);
  virtual float GetClipOffset() const {return clip_offset_;}

  virtual void ColorBy(const mol::EntityView& ev, 
                       const String& prop,
                       const Gradient& g, float minv, float maxv);
  virtual void ColorBy(const img::MapHandle& mh,
                       const String& prop,
                       const Gradient& g,float minv, float maxv);

  // end of gfx obj base interface

  // new gfx obj virtual interface starts here

  /*!
    \brief returns the bounding box of this object

    The bounding box, i.e. the left-bottom-front and the right-top-back corner
    of the object, is obtained with this method. The single boolean parameter
    denotes whether to return local coordinates or global scene coordinates; 
    global scene coordinates are local coordinates with the object's transformation
    (if present) applied. The default is false, i.e. return local coordinates.
  */
  virtual geom::AlignedCuboid GetBoundingBox(bool use_tf=false) const;

  /*!
    \brief adjust minimum and maximum extent based on graphical object

    this routine will adjust the provided minimum and maximum points
    based on the vertices of the underlying graphical representation,
    combining the given Transform with the object transform.

    If derived classes do not implement this method then the limits
    will be adjusted based on the Cuboid returned by GetBoundingBox(true)
  */
  virtual void ProcessLimits(geom::Vec3& minc, geom::Vec3& maxc, 
                             const geom::Transform& tf) const;

  // implemented in derived classes for the actual GL rendering
  /*
    note: only STANDARD_RENDER_PASS and GLOW_RENDER_PASS need to 
    be implemented, the rest is taken care of by GfxObj::RenderGL
  */
  virtual void CustomRenderGL(RenderPass pass);

  // implemented in derived classes to deal with initialization etc
  // called just before CustomRenderGL is called
  // the boolean flag indicated that a re-build was requested
  virtual void CustomPreRenderGL(bool rebuild);

  // implemented in derived classes for first GL initialization
  // which should be done here, not in the ctor
  virtual void InitGL();

  // implemented in derived classes for the actual POVray export
  virtual void CustomRenderPov(PovState& pov);

  // handle pick (select) event
  /*
    should return the position of the picked object in order for the caller
    to determine the closest one, and indicate with the boolean return wether
    the selection was successful
    the boolean flag indicates whether a true pick or just a feedback is requested
    the given zlim is the projected z of the best selection so far, and should be
    used for both picking as well as feedback if the selection is in front of the
    currently best one
  */
  virtual bool OnSelect(const geom::Line3& l, geom::Vec3& result, float zlim, bool pick_flag);

  // input event entry point
  virtual void OnInput(const InputEvent& e);

  // informs derived class that render mode has changes
  virtual void OnRenderModeChange();
  virtual void OnGLCleanup();
  //

  /// \brief removes all graphical elements
  void Clear();
  
  /// \brief get transform
  const geom::Transform& GetTF() const;
  /// \brief set transform
  void SetTF(const geom::Transform& tf);
  
  // add a label at the given position
  void AddLabel(const String& s, const geom::Vec3& pos, const Color& col, float psize);
  // convenience method
  void AddLabel(const String& s, const geom::Vec3& pos, const Color& col);
  // convenience method
  void AddLabel(const String& s, const geom::Vec3& pos, float psize);
  // convenience method
  void AddLabel(const String& s, const geom::Vec3& pos);
  // remove all labels
  void ClearLabels();

  void FlagRebuild();
  void FlagRefresh();

  bool GetAALines() const {return aalines_flag_;}
  float GetLineWidth() const {return line_width_;}
  float GetLineHalo() const {return line_halo_;}

  void SetNormalSmoothFactor(float smoothf);
  float GetNormalSmoothFactor() const;

  Material GetMaterial() const;
  void SetMaterial(const Material& m);

  // experimental, don't use
  void SmoothVertices(float smoothf);
 
  void GLCleanup();

  // this really should not be here
  static Color Ele2Color(const String& ele);

  void Debug(unsigned int flags);

  IndexedVertexArray& GetVA() {return va_;}
  const IndexedVertexArray& GetVA() const {return va_;}

 protected:
  
  void PreRenderGL(bool flag);

 private:
  GfxObj(const GfxObj& o);
  GfxObj& operator=(const GfxObj&);

 protected:
  void AppendColorOp(gfx::ColorOp* op);
  void CleanColorOps();
  void ReapplyColorOps();
  void render_labels() const;
  void render_depth_only();

  IndexedVertexArray va_;
  RenderMode::Type render_mode_;
  unsigned int debug_flags_;
 
  geom::Transform transform_;
  bool rebuild_;
  bool refresh_;
  float line_width_;
  unsigned int poly_mode_;
  bool aalines_flag_;
  float line_halo_;

  Material mat_;
  unsigned int mat_dlist_;
  bool mat_update_;

  float opacity_;
  float smoothf_;
  bool outline_flag_;
  int outline_mode_;

  bool solid_;
  Color solid_color_;

  bool clip_flag_;
  geom::Vec4 clip_plane_;
  float clip_offset_;

  boost::ptr_vector<gfx::ColorOp> c_ops_;

  TextPrimList labels_;

  bool use_occlusion_;
};

}} //ns

#endif
