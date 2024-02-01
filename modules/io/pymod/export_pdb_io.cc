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

#include <ost/io/mol/io_profile.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/io/mol/pdb_writer.hh>
#include <ost/io/mol/pdb_str.hh>

using namespace ost;
using namespace ost::io;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(X_import,
                                       PDBReader::Import, 1, 2)

BOOST_PYTHON_FUNCTION_OVERLOADS(load_PDB_ov, LoadPDB, 1, 2)

void (PDBWriter::*write_a)(const mol::EntityHandle&)=&PDBWriter::Write;
void (PDBWriter::*write_b)(const mol::EntityView&)=&PDBWriter::Write;

String (*pdb_str_a)(const mol::EntityHandle&, const IOProfile&)=&EntityToPDBString;
String (*pdb_str_b)(const mol::EntityView&, const IOProfile&)=&EntityToPDBString;

void remove_profiles() {
  IOProfileRegistry::RemoveProfiles();
}

void export_pdb_io()
{
  class_<IOProfile>("IOProfile",
         init<String,bool,bool,bool,bool,bool,bool,
              conop::ProcessorPtr>((arg("dialect")="PDB",
                                    arg("quack_mode")=false,
                                    arg("fault_tolerant")=false,
                                    arg("join_spread_atom_records")=false,
                                    arg("no_hetatms")=false,
                                    arg("calpha_only")=false,
                                    arg("read_conect")=false,
                                    arg("processor")=conop::ProcessorPtr())))
    .def(init<const IOProfile&>())
    .def_readwrite("dialect", &IOProfile::dialect)
    .def_readwrite("fault_tolerant", &IOProfile::fault_tolerant)
    .def_readwrite("quack_mode", &IOProfile::quack_mode)
    .def_readwrite("no_hetatms", &IOProfile::no_hetatms)
    .def_readwrite("calpha_only", &IOProfile::calpha_only)
    .def_readwrite("join_spread_atom_records", &IOProfile::join_spread_atom_records)
    .def_readwrite("read_conect", &IOProfile::read_conect)
    .def_readwrite("processor", &IOProfile::processor)
    .def("Copy", &IOProfile::Copy)
    .def(self_ns::str(self))
  ;
  class_<IOProfileRegistry>("IOProfileRegistry", no_init)
    .def("Get", &IOProfileRegistry::Get,  
         return_internal_reference<>())
    .def("Set", &IOProfileRegistry::Set)
    .def("Instance", &IOProfileRegistry::Instance,
         //return_internal_reference).staticmethod("Instance")
         return_value_policy<reference_existing_object>()).staticmethod("Instance")
    .def("GetDefault", &IOProfileRegistry::GetDefault,
         return_internal_reference<>())
  ;
  class_<PDBReader, boost::noncopyable>("PDBReader", init<String, const IOProfile&>())
    .def("HasNext", &PDBReader::HasNext)
    .def("Import", &PDBReader::Import, 
         X_import(args("entity", "restrict_chains")))
    .add_property("read_seqres", &PDBReader::GetReadSeqRes, 
                  &PDBReader::SetReadSeqRes)
    .add_property("seqres", &PDBReader::GetSeqRes)
  ;
  
  class_<PDBWriter, boost::noncopyable>("PDBWriter", init<String, const IOProfile&>())
    .def("Write", write_a)
    .def("SetIsPQR", &PDBWriter::SetIsPQR)
    .def("IsPQR", &PDBWriter::IsPQR)
    .def("SetWriteMultiModel", &PDBWriter::SetWriteMultiModel)
    .def("GetWriteMultiModel", &PDBWriter::GetWriteMultiModel)
    .add_property("write_multi_model", &PDBWriter::GetWriteMultiModel, 
                  &PDBWriter::SetWriteMultiModel)
    .def("Write", write_b)    
  ;
  def("EntityToPDBStr", pdb_str_a,
      (arg("entity"), arg("profile")=IOProfile()));
  def("EntityToPDBStr", pdb_str_b,
      (arg("entity"), arg("profile")=IOProfile()));

  def("PDBStrToEntity", &PDBStringToEntity, (arg("pdb_string"),
                                             arg("profile")=IOProfile(),
                                             arg("process")=false));

  // we need to make sure there are no pending references to Python objects
  // tied to the IOProfileRegistry singleton. The destructor of 
  // IOProfileRegistry may be called after Python is shutdown which results
  // in a segfault.
  scope().attr("__dict__")["atexit"]=handle<>(PyImport_ImportModule("atexit"));

  def("_remove_profiles", &remove_profiles);
  object r=scope().attr("_remove_profiles");
  scope().attr("atexit").attr("register")(r);
}
