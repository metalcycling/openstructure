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

#include <ost/gfx/render_options/render_options.hh>
#include <ost/gfx/render_options/cpk_render_options.hh>
#include <ost/gfx/render_options/custom_render_options.hh>
#include <ost/gfx/render_options/cartoon_render_options.hh>
#include <ost/gfx/render_options/line_render_options.hh>
#include <ost/gfx/render_options/simple_render_options.hh>
#include <ost/gfx/render_options/line_trace_render_options.hh>
#include <ost/gfx/render_options/sline_render_options.hh>
#include <ost/gfx/render_options/trace_render_options.hh>
using namespace boost::python;

using namespace ost;
using namespace ost::gfx;

void export_RenderOptions()
{
  class_<RenderOptions, boost::shared_ptr<RenderOptions>, boost::noncopyable>("RenderOptions", no_init)
    .def("GetRenderMode",&RenderOptions::GetRenderMode)
    .def("NotifyStateChange",&RenderOptions::NotifyStateChange)
    .def("AddObserver",&RenderOptions::AddObserver)
    .def("RemoveObserver",&RenderOptions::RemoveObserver)
  ;

  class_<LineRenderOptions, boost::shared_ptr<LineRenderOptions>, bases<RenderOptions>, boost::noncopyable>("LineRenderOptions", no_init)
    .def("SetLineWidth",&LineRenderOptions::SetLineWidth)
    .def("GetLineWidth",&LineRenderOptions::GetLineWidth)
    .add_property("line_width", &LineRenderOptions::GetLineWidth, 
                  &LineRenderOptions::SetLineWidth)    
    .def("SetAALines",&LineRenderOptions::SetAALines)
    .def("GetAALines",&LineRenderOptions::GetAALines)
    .add_property("aa_lines", &LineRenderOptions::GetAALines, 
                  &LineRenderOptions::SetAALines)
  ;

  class_<SimpleRenderOptions, boost::shared_ptr<SimpleRenderOptions>, bases<LineRenderOptions>, boost::noncopyable>("SimpleRenderOptions")
    .def(init<const SimpleRenderOptions&>())
    .def("GetBondOrderFlag", &SimpleRenderOptions::GetBondOrderFlag)
    .def("SetBondOrderFlag", &SimpleRenderOptions::SetBondOrderFlag) 
    .add_property("draw_bond_order", &SimpleRenderOptions::GetBondOrderFlag, 
                  &SimpleRenderOptions::SetBondOrderFlag)
    .add_property("bond_order_flag", &SimpleRenderOptions::GetBondOrderFlag, 
                  &SimpleRenderOptions::SetBondOrderFlag)
    .def("GetBondOrderDistance", &SimpleRenderOptions::GetBondOrderDistance)
    .def("SetBondOrderDistance", &SimpleRenderOptions::SetBondOrderDistance)
    .add_property("bond_order_distance", &SimpleRenderOptions::GetBondOrderDistance, 
                  &SimpleRenderOptions::SetBondOrderDistance)    
  ;

  class_<LineTraceRenderOptions, boost::shared_ptr<LineTraceRenderOptions>, bases<LineRenderOptions>, boost::noncopyable>("LineTraceRenderOptions")
    .def(init<const LineTraceRenderOptions&>())
  ;

  class_<SlineRenderOptions, boost::shared_ptr<SlineRenderOptions>, bases<LineRenderOptions>, boost::noncopyable>("SlineRenderOptions")
    .def(init<const SlineRenderOptions&>())
    .def("SetSplineDetail", &SlineRenderOptions::SetSplineDetail)
    .def("GetSplineDetail", &SlineRenderOptions::GetSplineDetail)
    .add_property("spline_detail",&SlineRenderOptions::GetSplineDetail,&SlineRenderOptions::SetSplineDetail)
  ;

  class_<CPKRenderOptions, boost::shared_ptr<CPKRenderOptions>, bases<RenderOptions>, boost::noncopyable>("CPKRenderOptions")
    .def(init<const CPKRenderOptions&>())
    .def("SetSphereDetail", &CPKRenderOptions::SetSphereDetail)
    .def("GetSphereDetail", &CPKRenderOptions::GetSphereDetail)
    .add_property("sphere_detail", &CPKRenderOptions::GetSphereDetail, 
                  &CPKRenderOptions::SetSphereDetail)
    .def("SetSphereMode", &CPKRenderOptions::SetSphereMode)
    .def("GetSphereMode", &CPKRenderOptions::GetSphereMode)
    .add_property("sphere_mode", &CPKRenderOptions::GetSphereMode, 
                  &CPKRenderOptions::SetSphereMode)
    .def("SetCPKMode", &CPKRenderOptions::SetSphereMode)
    .def("GetCPKMode", &CPKRenderOptions::GetSphereMode)
    .add_property("cpk_mode", &CPKRenderOptions::GetSphereMode, 
                  &CPKRenderOptions::SetSphereMode)
    .add_property("radius_mult",&CPKRenderOptions::GetRadiusMult,&CPKRenderOptions::SetRadiusMult)
  ;

  class_<CustomRenderOptions, boost::shared_ptr<CustomRenderOptions>, bases<RenderOptions>, boost::noncopyable>("CustomRenderOptions")
    .def(init<const CustomRenderOptions&>())
    .def("SetSphereDetail", &CustomRenderOptions::SetSphereDetail)
    .def("GetSphereDetail", &CustomRenderOptions::GetSphereDetail)
    .add_property("sphere_detail", &CustomRenderOptions::GetSphereDetail, &CustomRenderOptions::SetSphereDetail)
    .def("SetArcDetail", &CustomRenderOptions::SetArcDetail)
    .def("GetArcDetail", &CustomRenderOptions::GetArcDetail)
    .add_property("arc_detail", &CustomRenderOptions::GetArcDetail, 
                  &CustomRenderOptions::SetArcDetail)
    .def("SetSphereRad", &CustomRenderOptions::SetSphereRad)
    .def("GetSphereRad", &CustomRenderOptions::GetSphereRad)
    .add_property("sphere_rad", &CustomRenderOptions::GetSphereRad, &CustomRenderOptions::SetSphereRad)    
    .add_property("sphere_radius", &CustomRenderOptions::GetSphereRad, &CustomRenderOptions::SetSphereRad)    
    .def("SetBondRad", &CustomRenderOptions::SetBondRad)
    .def("GetBondRad", &CustomRenderOptions::GetBondRad)
    .add_property("bond_rad", &CustomRenderOptions::GetBondRad, &CustomRenderOptions::SetBondRad)
    .add_property("bond_radius", &CustomRenderOptions::GetBondRad, &CustomRenderOptions::SetBondRad)
  ;

  class_<CartoonRenderOptions, boost::shared_ptr<CartoonRenderOptions>, bases<RenderOptions>, boost::noncopyable>("CartoonRenderOptions")
    .def(init<const CartoonRenderOptions&>())
    .def(init<optional<bool> >())
    .def("SetSplineDetail", &CartoonRenderOptions::SetSplineDetail)
    .def("GetSplineDetail", &CartoonRenderOptions::GetSplineDetail)
    .add_property("spline_detail", &CartoonRenderOptions::GetSplineDetail, 
                  &CartoonRenderOptions::SetSplineDetail)
    .def("SetPolyMode", &CartoonRenderOptions::SetPolyMode)
    .def("GetPolyMode", &CartoonRenderOptions::GetPolyMode)
    .add_property("poly_mode", &CartoonRenderOptions::GetPolyMode,&CartoonRenderOptions::SetPolyMode)
    .def("SetArcDetail", &CartoonRenderOptions::SetArcDetail)
    .def("GetArcDetail", &CartoonRenderOptions::GetArcDetail)
    .add_property("arc_detail", &CartoonRenderOptions::GetArcDetail, 
                  &CartoonRenderOptions::SetArcDetail)
    .def("SetNormalSmoothFactor", &CartoonRenderOptions::SetNormalSmoothFactor)
    .def("GetNormalSmoothFactor", &CartoonRenderOptions::GetNormalSmoothFactor)
    .add_property("normal_smooth_factor",&CartoonRenderOptions::GetNormalSmoothFactor,
                  &CartoonRenderOptions::SetNormalSmoothFactor)
    .def("SetTubeRadius", &CartoonRenderOptions::SetTubeRadius)
    .def("GetTubeRadius", &CartoonRenderOptions::GetTubeRadius)
    .add_property("tube_radius",&CartoonRenderOptions::GetTubeRadius,
                  &CartoonRenderOptions::SetTubeRadius)
    .def("SetTubeRatio", &CartoonRenderOptions::SetTubeRatio)
    .def("GetTubeRatio", &CartoonRenderOptions::GetTubeRatio)
    .add_property("tube_ratio",&CartoonRenderOptions::GetTubeRatio,
                  &CartoonRenderOptions::SetTubeRatio)
    .def("SetTubeProfileType", &CartoonRenderOptions::SetTubeProfileType)
    .def("GetTubeProfileType", &CartoonRenderOptions::GetTubeProfileType)
    .add_property("tube_profile_type",&CartoonRenderOptions::GetTubeProfileType,
                  &CartoonRenderOptions::SetTubeProfileType)
    .def("SetHelixWidth", &CartoonRenderOptions::SetHelixWidth)
    .def("GetHelixWidth", &CartoonRenderOptions::GetHelixWidth)
    .add_property("helix_width",&CartoonRenderOptions::GetHelixWidth,
                  &CartoonRenderOptions::SetHelixWidth)
    .def("SetHelixThickness", &CartoonRenderOptions::SetHelixThickness)
    .def("GetHelixThickness", &CartoonRenderOptions::GetHelixThickness)
    .add_property("helix_thickness",&CartoonRenderOptions::GetHelixThickness,
                  &CartoonRenderOptions::SetHelixThickness)
    .def("SetHelixEcc", &CartoonRenderOptions::SetHelixEcc)
    .def("GetHelixEcc", &CartoonRenderOptions::GetHelixEcc)
    .add_property("helix_ecc",&CartoonRenderOptions::GetHelixEcc,
                  &CartoonRenderOptions::SetHelixEcc)
    .def("SetHelixProfileType", &CartoonRenderOptions::SetHelixProfileType)
    .def("GetHelixProfileType", &CartoonRenderOptions::GetHelixProfileType)
    .add_property("helix_profile_type",&CartoonRenderOptions::GetHelixProfileType,
                  &CartoonRenderOptions::SetHelixProfileType)
    .def("SetHelixMode", &CartoonRenderOptions::SetHelixMode)
    .def("GetHelixMode", &CartoonRenderOptions::GetHelixMode)
    .add_property("helix_mode",&CartoonRenderOptions::GetHelixMode,
                  &CartoonRenderOptions::SetHelixMode)
    .def("SetStrandWidth", &CartoonRenderOptions::SetStrandWidth)
    .def("GetStrandWidth", &CartoonRenderOptions::GetStrandWidth)
    .add_property("strand_width",&CartoonRenderOptions::GetStrandWidth,
                  &CartoonRenderOptions::SetStrandWidth)
    .def("SetStrandThickness", &CartoonRenderOptions::SetStrandThickness)
    .def("GetStrandThickness", &CartoonRenderOptions::GetStrandThickness)
    .add_property("strand_thickness",&CartoonRenderOptions::GetStrandThickness,
                  &CartoonRenderOptions::SetStrandThickness)
    .def("SetStrandEcc", &CartoonRenderOptions::SetStrandEcc)
    .def("GetStrandEcc", &CartoonRenderOptions::GetStrandEcc)
    .add_property("strand_ecc",&CartoonRenderOptions::GetStrandEcc,
                  &CartoonRenderOptions::SetStrandEcc)
    .def("SetStrandProfileType", &CartoonRenderOptions::SetStrandProfileType)
    .def("GetStrandProfileType", &CartoonRenderOptions::GetStrandProfileType)
    .add_property("strand_profile_type",&CartoonRenderOptions::GetStrandProfileType,
                  &CartoonRenderOptions::SetStrandProfileType)
    .def("SetStrandMode", &CartoonRenderOptions::SetStrandMode)
    .def("GetStrandMode", &CartoonRenderOptions::GetStrandMode)
    .add_property("strand_mode",&CartoonRenderOptions::GetStrandMode,
                  &CartoonRenderOptions::SetStrandMode)
    .def("SetColorBlendMode", &CartoonRenderOptions::SetColorBlendMode)
    .def("GetColorBlendMode", &CartoonRenderOptions::GetColorBlendMode)
    .add_property("color_blend_mode",&CartoonRenderOptions::GetColorBlendMode,
                  &CartoonRenderOptions::SetColorBlendMode)
    .add_property("twist_hack",&CartoonRenderOptions::GetTwistHack,
                  &CartoonRenderOptions::SetTwistHack)
  ;
  
  class_<TraceRenderOptions, boost::shared_ptr<TraceRenderOptions>, bases<RenderOptions>, boost::noncopyable>("TraceRenderOptions")
    .def(init<const TraceRenderOptions&>())
    .def("SetArcDetail", &TraceRenderOptions::SetArcDetail)
    .def("GetArcDetail", &TraceRenderOptions::GetArcDetail)
    .add_property("arc_detail", &TraceRenderOptions::GetArcDetail, 
                  &TraceRenderOptions::SetArcDetail)
    .def("SetNormalSmoothFactor", &TraceRenderOptions::SetNormalSmoothFactor)
    .def("GetNormalSmoothFactor", &TraceRenderOptions::GetNormalSmoothFactor)
    .add_property("normal_smooth_factor", 
                  &TraceRenderOptions::GetNormalSmoothFactor, 
                  &TraceRenderOptions::SetNormalSmoothFactor)    
    .def("SetTubeRadius", &TraceRenderOptions::SetTubeRadius)
    .def("GetTubeRadius", &TraceRenderOptions::GetTubeRadius)
    .add_property("tube_radius", &TraceRenderOptions::GetTubeRadius, 
                  &TraceRenderOptions::SetTubeRadius)
  ;
}



