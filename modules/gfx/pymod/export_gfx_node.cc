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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;

#include <ost/geom/export_helper/vector.hh>
#include <ost/gfx/gfx_node.hh>


using namespace ost;
using namespace ost::gfx;


const GfxNodeVector& (GfxNode::*get_children)() const=&GfxNode::GetChildren;

void export_GfxNode()
{
  void (GfxNode::* node_add1)(GfxObjP) = &GfxNode::Add;
  void (GfxNode::* node_rem1)(GfxObjP) = &GfxNode::Remove;
  void (GfxNode::* node_add2)(GfxNodeP) = &GfxNode::Add;
  void (GfxNode::* node_rem2)(GfxNodeP) = &GfxNode::Remove;
  void (GfxNode::* node_rem3)(const String&) = &GfxNode::Remove;

  class_<GfxNodeVector>("GfxNodeList", init<>())
    .def(vector_indexing_suite<GfxNodeVector, true>())
    .def(geom::VectorAdditions<GfxNodeVector>())
  ;
  class_<GfxNode, GfxNodeP, 
         boost::noncopyable>("GfxNode", init<const String&>())
    .def("GetName",&GfxNode::GetName)
    .def("Hide",&GfxNode::Hide)
    .def("Show",&GfxNode::Show)
    .def("IsVisible",&GfxNode::IsVisible)
    .def("Add",node_add1)
    .def("Remove",node_rem1)
    .def("Add",node_add2)
    .add_property("name", &GfxNode::GetName)    
    .def("Remove",node_rem2)
    .def("Remove",node_rem3)
    .add_property("parent", &GfxNode::GetParent)
    .add_property("children", make_function(get_children, 
                                            return_value_policy<copy_const_reference>()))
    ;
}
