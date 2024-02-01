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

#include <ost/gfx/exporter.hh>
#include <ost/gfx/gost_exporter.hh>
#include <ost/gfx/collada_exporter.hh>

using namespace ost;
using namespace ost::gfx;

void export_Exporter()
{
  class_<Exporter, boost::noncopyable>("Exporter", no_init)
    .add_property("scale",&Exporter::GetScale,&Exporter::SetScale)
    .add_property("to_origin",&Exporter::GetToOrigin,&Exporter::SetToOrigin)
  ;

  // internal class, factory function in __init__.py
  class_<GostExporter, bases<Exporter>, boost::noncopyable>("GostExporter_", init<const std::string&>())
    ;

  // internal class, factory function in __init__.py
  class_<ColladaExporter, bases<Exporter>, boost::noncopyable>("ColladaExporter_", init<const std::string&>())
    ;
}
