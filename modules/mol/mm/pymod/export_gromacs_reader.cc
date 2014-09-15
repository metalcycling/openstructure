//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2011 by the OpenStructure authors
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
#include <ost/mol/mm/gromacs_reader.hh>
#include <ost/mol/residue_handle.hh>

using namespace boost::python;
using namespace ost::mol::mm;


void export_GromacsReader()
{
  class_<ost::mol::mm::GromacsReader>("GromacsReader",init<String>())
    .def("ReadGromacsForcefield",&ost::mol::mm::GromacsReader::ReadGromacsForcefield)
    .def("SetPreprocessorDefinition",&ost::mol::mm::GromacsReader::SetPreprocessorDefinition)
    .def("GetForcefield",&ost::mol::mm::GromacsReader::GetForcefield)
    .def("SetForcefield",&ost::mol::mm::GromacsReader::SetForcefield)
    .def("ReadResidueDatabase",&ost::mol::mm::GromacsReader::ReadResidueDatabase)
    .def("ReadITP",&ost::mol::mm::GromacsReader::ReadITP)   
    .def("ReadCHARMMPRM",&ost::mol::mm::GromacsReader::ReadCHARMMPRM) 
    .def("ReadCHARMMRTF",&ost::mol::mm::GromacsReader::ReadCHARMMRTF)                                                                      
  ;

  boost::python::register_ptr_to_python<ost::mol::mm::GromacsReaderPtr>();
}