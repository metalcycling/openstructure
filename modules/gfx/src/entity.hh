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
#ifndef OST_GFX_ENTITY_HH
#define OST_GFX_ENTITY_HH

/*
  Author: Ansgar Philippsen, Marco Biasini
*/
#include <vector>

#include <boost/ptr_container/ptr_map.hpp>

#include <ost/geom/geom.hh>

#include <ost/mol/query_view_wrapper.hh>
#include <ost/gfx/render_options/render_options.hh>
#include <ost/gfx/color_ops/color_op.hh>
#include <ost/gfx/color_ops/by_element_color_op.hh>
#include <ost/gfx/color_ops/by_chain_color_op.hh>
#include <ost/gfx/color_ops/uniform_color_op.hh>
#include <ost/gfx/color_ops/gradient_level_color_op.hh>
#include <ost/gfx/color_ops/entity_view_color_op.hh>
#include <ost/gfx/color_ops/map_handle_color_op.hh>

#include <ost/gfx/impl/entity_renderer.hh>

#include "gradient.hh"
#include "entity_fw.hh"
#include "impl/backbone_trace.hh"
#include "impl/entity_detail.hh"

namespace ost { namespace gfx {

typedef std::vector<RenderMode::Type> RenderModeTypes;

/// \brief graphical rendering of \ref mol::EntityHandle entites
///
/// Entity is responsible for rendering of 
/// \ref mol::EntityHandle "molecular entities". A bunch of different 
/// \ref RenderMode::Type render modes are supported. To switch the render mode
/// of the full entity, Entity::SetRenderMode(RenderMode::Type) can be used. To 
/// only change one part to a different rendering style, use 
/// Entity::SetRenderMode(RenderMode::Type, mol::EntityView, bool)
class DLLEXPORT_OST_GFX Entity: public GfxObj {

public:
  /// \brief Initialize with an object name, an mol::Entity handle, and 
  ///    optionally with a mol::Query
  /// 
  /// If the latter is ommitted, the full view is displayed. A later call to 
  /// Rebuild() will re-apply the mol::Query to the given mol::EntityHandle
  Entity(const String& name,
         const mol::EntityHandle& eh,
         const mol::Query& q=mol::Query(),
         mol::QueryFlags f=0);

  /// \brief variant with explicit graphics mode instead of the default
  Entity(const String& name,
         RenderMode::Type m,
         const mol::EntityHandle& eh,
         const mol::Query& q=mol::Query(),
         mol::QueryFlags f=0);

  /// \brief Initialize with an object name and an explicit mol::Entity view; 
  ///    later calls to Rebuild always use this mol::Entity view
  Entity(const String& name,
         const mol::EntityView& ev);

  /// \brief ctor variant with explicit graphics mode instead of the default
  Entity(const String& name,
         RenderMode::Type m,
         const mol::EntityView& ev);
  virtual ~Entity();

  virtual geom::AlignedCuboid GetBoundingBox(bool use_tf=false) const;

  // ProcessLimits uses the default implementation of bounding box
  
  /// internal routine
  virtual void CustomRenderGL(RenderPass pass);

  virtual void RefreshVA();
  
  virtual bool OnSelect(const geom::Line3& l, geom::Vec3& result, float zlim,
                        bool pick_flag);

  /// \brief pick atom
  /// 
  /// In case the line intersects several atoms, the atom closest to the 
  /// near clipping plane is returned. Returns an invalid handle in case no 
  /// atom was close to the line
  /// \todo honour object transformation
  mol::AtomHandle PickAtom(const geom::Line3& line, Real line_width=0.5);
  

  /// \brief pick bond
  /// 
  /// In case the line intersects several bonds, the bond closest to the 
  /// near clipping plane is returned. Returns an invalid handle if no bond was 
  /// close to the line.
  /// \todo honour object transformation
  mol::BondHandle PickBond(const geom::Line3& line, Real line_width=0.5);


  virtual void OnRenderModeChange();
  
  const String GetRenderModeName(RenderMode::Type mode);

  void SetEnableRenderMode(RenderMode::Type mode, bool enable);

  bool IsRenderModeEnabled(RenderMode::Type mode);

  RenderModeTypes GetNotEmptyRenderModes();

  void SetRenderMode(RenderMode::Type mode, const mol::EntityView& view, 
                     bool keep=false);
  void SetRenderMode(RenderMode::Type mode, const String& selection, 
                     bool keep=false);
  virtual void SetRenderMode(RenderMode::Type mode);  
  
  mol::EntityView GetRenderView(RenderMode::Type mode);

  virtual void SetVisible(const mol::EntityView& view, bool visible);

  virtual void SetVisible(const String& sel, bool visible);
  virtual void OptionsChanged(RenderMode::Type mode);

  virtual void SetOpacity(float f);
  virtual float GetOpacity() const {return opacity_;}
  virtual void SetOutlineWidth(float f);
  virtual void SetOutlineExpandFactor(float f);
  virtual void SetOutlineExpandColor(const Color& c);
  virtual void SetClipOffset(float f);

  /// \brief resets used entity handle
  /// replaces underlying entity, keeps query and flags intact
  void Reset(const mol::EntityHandle& eh);
  /// \brief resets used entity handle and query
  /// replaces underlying entity and query, keeps flags intact
  void Reset(const mol::EntityHandle& eh, const mol::Query& q);
  /// \brief resets used entity handle, query and flags
  /// this has the same effect as the ctor call with the same parameters
  void Reset(const mol::EntityHandle& eh, const mol::Query& q, mol::QueryFlags flags);
  /// \brief resets entity view
  /// this as the same effect as the ctor call with the same parameters
  void Reset(const mol::EntityView& ev);
  /// \brief rebuild graphical object (see ctor comments)
  /*
    the naming here is misleading - this method WON'T be called upon FlagRebuild
  */
  void Rebuild();

  /// \brief only grab updated positions, dont rebuild the whole thing
  /// views won't be regenerated from stored queries
  void UpdatePositions();

  /// \brief forces all views to be regenerated from stored queries
  void UpdateView();

  /// \brief set color for selection
  void SetColor(const Color& col, const String& selection=String(""));

  // \brief detail coloring
  void SetDetailColor(const Color& col, const String& selection=String(""));
  
  /// \brief set color for specific atom
  void SetColorForAtom(const Color& col, 
                       const mol::AtomHandle& atom);

  /// \brief color by element
  void ColorByElement();
  
  /// \brief color by element for a specific selection
  void ColorByElement(const String& selection);

  /// \brief color by chain
  void ColorByChain();

  /// \brief color by chain for a specific selection
  void ColorByChain(const String& selection);

  /// \brief get view
  mol::EntityView GetView() const;

  /// \brief set a new query to use (deprecated)
  /// this will re-create the object based on the given selection
  void SetQuery(const mol::Query& q);

  /// return internally used query view
  mol::QueryViewWrapper GetQueryView() const;
  /// set new query view, rebuilding object
  void SetQueryView(const mol::QueryViewWrapper& qv);

  /// return underlying entity
  mol::EntityHandle GetEntity() const;

  // turn blur on or off (experimental feature)
  void SetBlur(bool f);
  // set atom positions as n-1 for blur (experimental feature)
  void BlurSnapshot();
  // blur transparency falloffs (experimental feature)
  void SetBlurFactors(float bf1,float bf2);

  /// \brief set selection
  /// \sa gfx_ent
  void SetSelection(const mol::EntityView& view);
  
  /// \brief get selection
  /// \sa gfx_ent
  mol::EntityView GetSelection() const;
  
  // GfxObj property interface
  virtual void ColorBy(const mol::EntityView& ev, 
                       const String& prop,
                       const Gradient& g, float minv, float maxv);

  // GfxObj property interface
  virtual void ColorBy(const img::MapHandle& mh,
                       const String& prop,
                       const Gradient& g,float minv, float maxv);

  // map property to color gradient from minv to maxv
  void ColorBy(const String& prop, 
               const Gradient& gradient,
               float minv,float maxv,
               mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  // temporarily here, will be moved to py interface
  void ColorBy(const String& prop, 
               const Gradient& gradient,
               float minv,float maxv,
               bool clamp);

  // temporary, should be incorporated with ColorBy
  void DetailColorBy(const String& prop, 
                     const Gradient& gradient,
                     float minv,float maxv,
                     mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  // convenience
  void ColorBy(const String& prop, 
               const Gradient& gradient,
               mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  // convenience
  void ColorBy(const String& prop,
               const Gradient& gradient,
               const String& selection);

  // convenience
  void ColorBy(const String& prop, 
               const Color& c1, const Color& c2, 
               float min, float max,
               mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  // convenience
  void ColorBy(const String& prop, 
               const Color& c1, const Color& c2,
               mol::Prop::Level hint=mol::Prop::UNSPECIFIED);


  void RadiusBy(const String& prop,
                float rmin, float rmax, 
                float vmin, float vmax,
                mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  void RadiusBy(const String& prop,
                float rmin, float rmax,
                mol::Prop::Level hint=mol::Prop::UNSPECIFIED);

  void ResetRadiusBy();

  void Apply(const gfx::ByElementColorOp& op, bool store=true);
  void Apply(const gfx::ByChainColorOp& op, bool store=true);
  void Apply(const gfx::UniformColorOp& op, bool store=true);
  void Apply(const gfx::GradientLevelColorOp& op, bool store=true);
  void Apply(const gfx::EntityViewColorOp& op, bool store=true);
  void Apply(const gfx::MapHandleColorOp& op, bool store=true);

  void CleanColorOps();
  void ReapplyColorOps();
  
  /// \brief Get render options for given render mode
  /// 
  /// In Python, also available as the properties:
  /// \li \c sline_options
  /// \li \c simple_options
  /// \li \c tube_options
  /// \li \c cartoon_options
  /// \li \c cpk_options
  /// \li \c trace_options
  /// \li \c custom_options
  RenderOptionsPtr GetOptions(RenderMode::Type render_mode);
  void SetOptions(RenderMode::Type render_mode, 
                  RenderOptionsPtr& render_options);
  void ApplyOptions(RenderMode::Type render_mode,
                          RenderOptionsPtr& render_options);
  bool HasSelection() const;

  void SetSeqHack(bool b);
  bool GetSeqHack() const;
  
  virtual void Export(Exporter* ex);

protected:

  virtual void CustomPreRenderGL(bool flag);
  virtual void CustomRenderPov(PovState& pov);
  void UpdateSelection();
  bool UpdateIfNeeded() const;
  void CacheBoundingBox() const;
  impl::EntityRenderer* GetOrCreateRenderer(RenderMode::Type);
private:
  mol::QueryViewWrapper         qv_;
                                
  mutable mol::EntityView       cached_view_;
  mutable bool                  update_view_;
  
  mutable geom::AlignedCuboid   bbox_;
  mol::EntityView               sel_;
  bool                          sel_update_;
  mutable impl::BackboneTrace   trace_;
                           
  void init(RenderMode::Type);

  void set_static_max_rad();
  void do_update_view() const;
  
  typedef boost::ptr_map<RenderMode::Type, impl::EntityRenderer> RendererMap;
  mutable RendererMap renderer_;

  float opacity_;
  bool blur_;
  float blurf1_;
  float blurf2_;
  mutable bool needs_update_;
};


/// \example load_and_display.py
/// 
/// Shows how to  display one \ref ost::gfx::Entity "entity" with several render modes
/// at once. The sidechains  are displayed simple mode, whereas the backbone is 
/// displayed with smooth lines.

///  \example rendermodes.py
/// 
/// Shows how to switch between different \ref ost::gfx::RenderMode "render modes" and
/// explains some of the rendermode parameters.
/// \sa \ref load_and_display.py "Loading and Displaying an Entity"

///  \example gfx_selection.py
/// 
/// Graphical entities have an active selection. In essence this selection is a 
/// subset of atoms, residues, chains and bonds. The active selection is 
/// displayed with a green halo around the structure. The selection can be set
/// either programmatically by using Entity::SetSelection() or interactively by 
/// using the selection tool.
/// 
/// \sa \ref load_and_display.py "Loading and Displaying an Entity"


/// \example color_by_property.py
/// 
/// Color \ref ost::gfx::Entity "graphical entity" by property using a gradient
/// \sa \ref gradient.py "Gradient Example"
}} // ns

#endif
