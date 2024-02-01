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

#ifndef OST_GFX_CARTOON_RENDERER_HH
#define OST_GFX_CARTOON_RENDERER_HH

#include <ost/geom/geom.hh>

#include <ost/gfx/entity_fw.hh>

#include <ost/gfx/impl/trace_renderer_base.hh>
#include <ost/gfx/impl/entity_detail.hh>
#include <ost/gfx/render_options/cartoon_render_options.hh>
#include <ost/gfx/impl/backbone_trace.hh>

namespace ost { namespace gfx { namespace impl {

/// \brief display entity in cartoon mode
/// \internal
class DLLEXPORT_OST_GFX CartoonRenderer: public TraceRendererBase {
public:
  CartoonRenderer(BackboneTrace* trace, bool force_tube=false);

  virtual geom::AlignedCuboid GetBoundingBox() const;  
  
  virtual void PrepareRendering();

  virtual bool CanSetOptions(RenderOptionsPtr& render_options);
  virtual void SetOptions(RenderOptionsPtr& render_options);
  virtual RenderOptionsPtr GetOptions();

  virtual void SetForceTube(bool force_tube);
  
private:
  void prepare_rendering(const BackboneTrace&,
                         IndexedVertexArray&, 
                         SplineEntryListList&);

  void fudge_spline_obj(SplineEntryListList&);

  void rebuild_spline_obj(IndexedVertexArray&,
                          const SplineEntryListList&,
                          bool);
  
  void cap_profile(const impl::TraceProfile&, 
                   const impl::SplineEntry&, 
                   bool,
                   IndexedVertexArray&);
                  
  void assemble_profile(const TraceProfile&,
                        const TraceProfile&, 
                        IndexedVertexArray&,
                        size_t);

  TraceProfile transform_and_add_profile(const std::vector<TraceProfile>&, 
                                         const SplineEntry&, 
                                         IndexedVertexArray&);

  TraceProfile get_circ_profile(unsigned int detail, float rx, float ry, unsigned int type, float ecc);

  bool force_tube_;
  CartoonRenderOptionsPtr options_;
  SplineEntryListList    spline_list_list_;
  SplineEntryListList    sel_spline_list_list_;
};

}}}

#endif // OST_GFX_CARTOON_RENDERER_HH
