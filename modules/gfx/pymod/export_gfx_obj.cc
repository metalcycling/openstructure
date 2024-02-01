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
#include <boost/python.hpp>
using namespace boost::python;

#include <ost/gfx/gfx_object.hh>
using namespace ost;
using namespace ost::gfx;

#include "color_by_def.hh"

  // convenience for python
  void set_mat_amb2(GfxObjBase* b, float c) {b->SetMatAmb(Color(c,c,c,1.0));}
  void set_mat_diff2(GfxObjBase* b, float c) {b->SetMatDiff(Color(c,c,c,1.0));}
  void set_mat_spec2(GfxObjBase* b, float c) {b->SetMatSpec(Color(c,c,c,1.0));}
  void set_mat_emm2(GfxObjBase* b, float c) {b->SetMatEmm(Color(c,c,c,1.0));}
  void set_mat1(GfxObjBase* b, float a, float d, float s, float p)
  {
    set_mat_amb2(b,a);
    set_mat_diff2(b,d);
    set_mat_spec2(b,s);
    b->SetMatShin(p);
    set_mat_emm2(b,0.0);
  }
  void set_mat2(GfxObjBase* b, Color a, Color d, Color s, float p)
  {
    b->SetMatAmb(a);
    b->SetMatDiff(d);
    b->SetMatSpec(s);
    b->SetMatShin(p);
    set_mat_emm2(b,0.0);
  }

  void set_outline(GfxObjBase* b, bool f)
  {
    LOG_INFO("Outline(bool) is deprecated, use SetOutline(bool) instead");
    b->SetOutline(f);
  }
  void set_aalines(GfxObjBase* b, bool f)
  {
    LOG_INFO("AALines(bool) is deprecated, use SetAALines(bool) instead");
    b->SetAALines(f);
  }

  class GfxObjWrap: public GfxObj, public wrapper<GfxObj>
  {
  public:
    GfxObjWrap(const std::string& name):
      GfxObj(name)
    {}

    virtual geom::AlignedCuboid GetBoundingBox(bool return_global=true) const
    {
      if(override f = this->get_override("GetBoundingBox")) {
        return f(return_global);
      } else {
        return GfxObj::GetBoundingBox(return_global);
      }
    }

    geom::AlignedCuboid default_GetBoundingBox(bool return_global) const {
      return GfxObj::GetBoundingBox(return_global);
    }

    virtual void CustomRenderGL(RenderPass pass) {
      if(override f = this->get_override("_CustomRenderGL")) {
        f(pass);
      } else {
        GfxObj::CustomRenderGL(pass);
      }
    }

    void default_CustomRenderGL(RenderPass pass) {
        GfxObj::CustomRenderGL(pass);
    }

    virtual void CustomPreRenderGL(bool rebuild) {
      if(override f = this->get_override("_CustomPreRenderGL")) {
        f(rebuild);
      } else {
        GfxObj::CustomPreRenderGL(rebuild);
      }
    }

    void default_CustomPreRenderGL(bool rebuild) {
        GfxObj::CustomPreRenderGL(rebuild);
    }

    virtual void InitGL() {
      if(override f = this->get_override("_InitGL")) {
        f();
      } else {
        GfxObj::InitGL();
      }
    }

    void default_InitGL() {
        GfxObj::InitGL();
    }
  };

void export_GfxObj()
{
  class_<GfxObjBase, boost::shared_ptr<GfxObjBase>, bases<GfxNode>, boost::noncopyable>("GfxObjBase",no_init)
    .def("SetMatAmb",&GfxObjBase::SetMatAmb)
    .def("SetMatAmb",set_mat_amb2)
    .def("SetMatDiff",&GfxObjBase::SetMatDiff)
    .def("SetMatDiff",set_mat_diff2)
    .def("SetMatSpec",&GfxObjBase::SetMatSpec)
    .def("SetMatSpec",set_mat_spec2)
    .def("SetMatEmm",&GfxObjBase::SetMatEmm)
    .def("SetMatEmm",set_mat_emm2)
    .def("SetMatShin",&GfxObjBase::SetMatShin)
    .def("SetMat",set_mat1)
    .def("SetMat",set_mat2)
    .def("ContextSwitch", &GfxObjBase::ContextSwitch)
    .def("SetRenderMode", &GfxObjBase::SetRenderMode)
    .def("GetRenderMode", &GfxObjBase::GetRenderMode)
    .def("GetCenter",&GfxObjBase::GetCenter) 
    .add_property("center", &GfxObjBase::GetCenter)
    .def("SetLineWidth", &GfxObjBase::SetLineWidth)
    .def("SetPolyMode",&GfxObjBase::SetPolyMode)
    .def("AALines",set_aalines) /* deprecated */
    .def("SetAALines",&GfxObjBase::SetAALines)
    .def("SetLineHalo",&GfxObjBase::SetLineHalo)
    .def("Outline",set_outline) /* deprecated */
    .def("SetOutline",&GfxObjBase::SetOutline)
    .def("GetOutline",&GfxObjBase::GetOutline)
    .add_property("outline",&GfxObjBase::GetOutline,&GfxObjBase::SetOutline)
    .def("SetOutlineMode",&GfxObjBase::SetOutlineMode)
    .add_property("outline_mode",&GfxObjBase::GetOutlineMode,&GfxObjBase::SetOutlineMode)
    .def("SetOutlineWidth",&GfxObjBase::SetOutlineWidth)
    .add_property("outline_width",&GfxObjBase::GetOutlineWidth,&GfxObjBase::SetOutlineWidth)
    .def("SetOutlineExpandFactor",&GfxObjBase::SetOutlineExpandFactor)
    .add_property("outline_expand_factor",&GfxObjBase::GetOutlineExpandFactor,&GfxObjBase::SetOutlineExpandFactor)
    .def("SetOutlineExpandColor",&GfxObjBase::SetOutlineExpandColor)
    .add_property("outline_expand_color",&GfxObjBase::GetOutlineExpandColor,&GfxObjBase::SetOutlineExpandColor)
    .add_property("outline_color",&GfxObjBase::GetOutlineExpandColor,&GfxObjBase::SetOutlineExpandColor)
    .def("SetOpacity",&GfxObjBase::SetOpacity)
    .def("GetOpacity",&GfxObjBase::GetOpacity)
    .add_property("opacity",&GfxObjBase::GetOpacity,&GfxObjBase::SetOpacity)
    .add_property("solid",&GfxObjBase::GetSolid,&GfxObjBase::SetSolid)
    .add_property("solid_color",&GfxObjBase::GetSolidColor,&GfxObjBase::SetSolidColor)
    .add_property("clip",&GfxObjBase::GetClip,&GfxObjBase::SetClip)
    .add_property("clip_plane",&GfxObjBase::GetClipPlane,&GfxObjBase::SetClipPlane)
    .add_property("clip_offset",&GfxObjBase::GetClipOffset,&GfxObjBase::SetClipOffset)
    COLOR_BY_DEF()
   ;

  enum_<RenderPass>("RenderPass")
    .value("STANDARD_RENDER_PASS",STANDARD_RENDER_PASS)
    .value("TRANSPARENT_RENDER_PASS",TRANSPARENT_RENDER_PASS)
    ;        

  class_<GfxObjWrap, boost::shared_ptr<GfxObj>, bases<GfxObjBase>, boost::noncopyable>("GfxObj",init<const std::string&>())
    .def("GetTF", &GfxObj::GetTF, return_value_policy<copy_const_reference>())
    .def("SetTF", &GfxObj::SetTF)
    .def("FlagRebuild",&GfxObj::FlagRebuild)
    .def("FlagRefresh",&GfxObj::FlagRefresh)
    .def("SetNormalSmoothFactor",&GfxObj::SetNormalSmoothFactor)
    .def("GetNormalSmoothFactor",&GfxObj::GetNormalSmoothFactor)
    .def("SmoothVertices",&GfxObj::SmoothVertices)
    .def("Debug",&GfxObj::Debug)
    .def("GetAALines",&GfxObj::GetAALines)
    .def("GetLineWidth",&GfxObj::GetLineWidth)
    .def("GetLineHalo",&GfxObj::GetLineHalo)
    .def("GetBoundingBox",&GfxObj::GetBoundingBox,&GfxObjWrap::default_GetBoundingBox)
    .def("_CustomRenderGL",&GfxObj::CustomRenderGL, &GfxObjWrap::default_CustomRenderGL)
    .def("_CustomPreRenderGL",&GfxObj::CustomPreRenderGL, &GfxObjWrap::default_CustomPreRenderGL)
    .def("_InitGL",&GfxObj::InitGL, &GfxObjWrap::default_InitGL)
    ;    

}
