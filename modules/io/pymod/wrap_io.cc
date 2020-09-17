//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <ost/io/io_manager.hh>
#include <ost/io/mol/load_entity.hh>
#include <ost/io/mol/load_surface.hh>
#include <ost/io/mol/save_entity.hh>
#include <ost/io/seq/load.hh>
#include <ost/io/seq/save.hh>
#include <ost/io/mol/entity_io_pdb_handler.hh>
#include <ost/io/mol/entity_io_crd_handler.hh>
#include <ost/io/mol/entity_io_pqr_handler.hh>
#include <ost/io/mol/entity_io_mae_handler.hh>
#include <ost/io/mol/entity_io_sdf_handler.hh>
#include <ost/io/mol/pdb_reader.hh>
#include <ost/io/mol/dcd_io.hh>
#include <ost/io/stereochemical_params_reader.hh>
using namespace ost;
using namespace ost::io;

namespace {

BOOST_PYTHON_FUNCTION_OVERLOADS(load_entity_ov,LoadEntity,1,2);
BOOST_PYTHON_FUNCTION_OVERLOADS(load_mentity_ov,LoadManagedEntity,2,3);

BOOST_PYTHON_FUNCTION_OVERLOADS(load_surface_ov,LoadSurface,1,2);
BOOST_PYTHON_FUNCTION_OVERLOADS(load_msurface_ov,LoadManagedSurface,2,3);


void save_ent_view(const mol::EntityView& en, const String& filename,
                   const String& format="auto") {
  SaveEntity(en, filename, format);
}

void save_ent_handle(const mol::EntityHandle& en, const String& filename,
                     const String& format="auto") {
  SaveEntity(en, filename, format);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(save_entity_handle_ov,
                                save_ent_handle, 2, 3)

BOOST_PYTHON_FUNCTION_OVERLOADS(save_entity_view_ov,
                                save_ent_view, 2, 3)

ost::mol::alg::StereoChemicalProps (*read_props_a)(String filename, bool check) = &ReadStereoChemicalPropsFile;
ost::mol::alg::StereoChemicalProps (*read_props_b)(bool check) = &ReadStereoChemicalPropsFile;

}

void export_pdb_io();
void export_mmcif_io();
#if OST_IMG_ENABLED
void export_map_io();
#endif
BOOST_PYTHON_MODULE(_ost_io)
{
  class_<IOManager, boost::noncopyable>("IOManager", no_init)
  ;

  // global functions
  def("Instance",IOManager::Instance,
      return_value_policy<reference_existing_object>() );

  def("LoadEntity", LoadEntity,
      load_entity_ov(args("filename", "format")));
  def("SaveEntity", &save_ent_view,
      save_entity_view_ov(args("entity", "filename", "format")));
  def("SaveEntity", &save_ent_handle,
      save_entity_handle_ov(args("entity", "filename", "format")));

  def("LoadAlignment", &LoadAlignment,
      (arg("filename"), arg("format")="auto"));
  def("AlignmentFromString", &AlignmentFromString);
  def("AlignmentFromStream", &AlignmentFromStream);
  def("AlignmentToString", &AlignmentToString);
  def("LoadSequenceList", &LoadSequenceList,
      (arg("filename"), arg("format")="auto"));
  def("LoadSequence", &LoadSequence,
      (arg("filename"), arg("format")="auto"));
  def("SequenceListFromString", &SequenceListFromString);
  def("SequenceFromString", &SequenceFromString);  
  def("SaveAlignment", &SaveAlignment,
      (arg("aln"), arg("filename"), arg("format")="auto"));
  
  def("LoadSequenceProfile", &LoadSequenceProfile,
      (arg("filename"), arg("format")="auto"));
  def("SequenceProfileFromString", &SequenceProfileFromString,
      (arg("data"), arg("format")));

  def("LoadSurface",LoadSurface,load_surface_ov());
  def("LoadManagedSurface",LoadManagedSurface,load_msurface_ov());

  def("SequenceToString", &SequenceToString);
  def("SequenceListToString", &SequenceListToString); 
  def("SaveSequenceList", &SaveSequenceList,
      (arg("seq_list"), arg("filename"), arg("format")="auto"));
  def("SaveSequence", &SaveSequence,
      (arg("sequence"), arg("filename"), arg("format")="auto"));
  def("LoadSDF", &LoadSDF);

  def("LoadCRD", &LoadCRD);
  def("LoadCHARMMTraj_", &LoadCHARMMTraj, (arg("ent"), arg("trj_filename"), 
                                           arg("stride")=1, arg("lazy_load")=false,
                                           arg("detect_swap")=true,arg("swap_bytes")=false))
;
  def("LoadMAE", &LoadMAE);
  def("LoadPQR", &LoadPQR);

  def("ReadStereoChemicalPropsFile", read_props_a,
    (arg("filename"), arg("check")=true)); 
  def("ReadStereoChemicalPropsFile", read_props_b,
    (arg("check")=true));
  def("GetStereoChemicalPropsFile", &GetStereoChemicalPropsFile);

  export_pdb_io();
  export_mmcif_io();
#if OST_IMG_ENABLED  
  export_map_io();
#endif
  def("SaveCHARMMTraj", &SaveCHARMMTraj, 
      (arg("traj"), arg("pdb_filename"), arg("dcd_filename"), arg("stride")=1, 
       arg("profile")=IOProfile()));
}
