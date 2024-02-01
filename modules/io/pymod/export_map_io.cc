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
#include  <ost/io/img/map_io_df3_handler.hh>
#include  <ost/io/img/map_io_dx_handler.hh>
#include  <ost/io/img/map_io_spi_handler.hh>
#include  <ost/io/img/map_io_mrc_handler.hh>
#include  <ost/io/img/map_io_dm3_handler.hh>
#include  <ost/io/img/map_io_situs_handler.hh>
#include  <ost/io/img/map_io_tiff_handler.hh>
#include  <ost/io/img/map_io_png_handler.hh>
#include  <ost/io/img/map_io_dat_handler.hh>
#include  <ost/io/img/map_io_jpk_handler.hh>
#include  <ost/io/img/map_io_nanoscope_handler.hh>
#include  <ost/io/img/map_io_ipl_handler.hh>
#include  <ost/io/img/image_format.hh>
#include  <ost/io/img/load_map.hh>

using namespace boost::python;
using namespace ost;
using namespace ost::io;

void  save_image1(const img::ImageHandle& image,const boost::filesystem::path& loc)
{
  SaveImage(image,loc);
}

void  save_image2(const img::ImageHandle& image,const boost::filesystem::path& loc, const ImageFormatBase& formatstruct)
{
  SaveImage(image,loc,formatstruct);
}

img::ImageHandle  load_image1(const String& loc)
{
  return LoadImage(loc);
}

img::ImageHandle  load_image2(const String& loc, const ImageFormatBase& formatstruct)
{
  return LoadImage(loc,formatstruct);
}

void export_map_io()
{
  class_<boost::logic::tribool>("tribool", init<boost::logic::tribool>())
    .def(init<bool>())
  ;
  implicitly_convertible<bool,boost::logic::tribool>();

  enum_<Endianess>("Endianess")
        .value("OST_BIG_ENDIAN", OST_BIG_ENDIAN)
        .value("OST_LITTLE_ENDIAN", OST_LITTLE_ENDIAN)
        .value("OST_VAX_DATA", OST_VAX_DATA)
        .value("OST_LOCAL_ENDIAN", OST_LOCAL_ENDIAN)
        .export_values()
  ;

  enum_<Format>("Format")
        .value("OST_BIT8_FORMAT", OST_BIT8_FORMAT)
        .value("OST_BIT16_FORMAT", OST_BIT16_FORMAT)
        .value("OST_BIT32_FORMAT", OST_BIT32_FORMAT)
        .value("OST_FLOAT_FORMAT", OST_FLOAT_FORMAT)
        .value("OST_DOUBLE_FORMAT", OST_DOUBLE_FORMAT)
        .value("OST_DEFAULT_FORMAT", OST_DEFAULT_FORMAT)
        .export_values()
  ;

  enum_<Subformat>("Subformat")
        .value("MRC_NEW_FORMAT", MRC_NEW_FORMAT)
        .value("MRC_OLD_FORMAT", MRC_OLD_FORMAT)
        .value("MRC_AUTO_FORMAT", MRC_AUTO_FORMAT)
        .export_values()
  ;


  class_<ImageFormatBase>("ImageFormatBase",no_init)
    .def("GetMaximum", &ImageFormatBase::GetMaximum)
    .def("GetMinimum", &ImageFormatBase::GetMinimum)
  ;
  
  class_<DF3, bases<ImageFormatBase> >("DF3", init<bool>(arg("normalize_on_save") = false))
    .def("SetNormalizeOnSave", &DF3::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &DF3::GetNormalizeOnSave)
  ;

  class_<DX, bases<ImageFormatBase> >("DX", init<bool>(arg("normalize_on_save") = false))
    .def("SetNormalizeOnSave", &DX::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &DX::GetNormalizeOnSave)
  ;

  class_<Spider, bases<ImageFormatBase> >("Spider", init<bool,Endianess>
               ((arg("normalize_on_save") = false, arg("endianess_on_save")=OST_LOCAL_ENDIAN)))
    .def("SetNormalizeOnSave", &Spider::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &Spider::GetNormalizeOnSave)
    .def("SetEndianessOnSave", &Spider::SetEndianessOnSave)
    .def("GetEndianessOnSave", &Spider::GetEndianessOnSave)
  ;

  class_<MRC, bases<ImageFormatBase> >("MRC", init<bool,Subformat,Endianess,Format>
           ((arg("normalize_on_save") = false,arg("subformat")=MRC_AUTO_FORMAT,arg("endianess_on_save")=OST_LOCAL_ENDIAN,arg("format")=OST_DEFAULT_FORMAT)))
    .def("SetBitDepth", &MRC::SetBitDepth)
    .def("GetBitDepth", &MRC::GetBitDepth)
    .def("SetNormalizeOnSave", &MRC::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &MRC::GetNormalizeOnSave)
    .def("SetSubformat", &MRC::SetSubformat)
    .def("GetSubformat", &MRC::GetSubformat)
    .def("SetEndianessOnSave", &MRC::SetEndianessOnSave)
    .def("GetEndianessOnSave", &MRC::GetEndianessOnSave)
  ;

  class_<CCP4, bases<MRC> >("CCP4", init<bool,Endianess,Format>
           ((arg("normalize_on_save") = false,arg("endianess_on_save")=OST_LOCAL_ENDIAN,arg("format")=OST_DEFAULT_FORMAT)))
  ;

  class_<DM3, bases<ImageFormatBase> >("DM3", init<>())
  ;

  class_<Situs, bases<ImageFormatBase> >("Situs", init<bool>(arg("normalize_on_save") = false))
    .def("SetNormalizeOnSave", &Situs::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &Situs::GetNormalizeOnSave)
  ;

  class_<TIF, bases<ImageFormatBase> >("TIF", init<boost::logic::tribool,Format,bool,bool,int>
				       ((arg("normalize_on_save") =  boost::logic::tribool(boost::logic::indeterminate),arg("format")=OST_DEFAULT_FORMAT,arg("signed")=false,arg("phasecolor")=false,arg("subimage") = -1)))
    .def("SetNormalizeOnSave", &TIF::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &TIF::GetNormalizeOnSave)
    .def("SetFormat", &TIF::SetFormat)
    .def("GetFormat", &TIF::GetFormat)
    .def("SetSigned", &TIF::SetSigned)
    .def("GetSigned", &TIF::GetSigned)
    .def("SetBitDepth", &TIF::SetBitDepth)
    .def("GetBitDepth", &TIF::GetBitDepth)
    .def("GetPhasecolor", &TIF::GetPhasecolor)
    .def("SetPhasecolor", &TIF::SetPhasecolor)
    .def("SetSubimage", &TIF::SetSubimage)
    .def("GetSubimage", &TIF::GetSubimage)
  ;

  class_<PNG, bases<ImageFormatBase> >("PNG", init<bool>(arg("normalize_on_save") = false))
    .def("SetNormalizeOnSave", &PNG::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &PNG::GetNormalizeOnSave)
  ;

  class_<DAT, bases<ImageFormatBase> >("DAT", init<bool,Format,bool,Endianess>((arg("normalize_on_save") = false,arg("format")=OST_DEFAULT_FORMAT,arg("signed")=false,arg("endianess") = OST_LOCAL_ENDIAN)))
    .def("SetNormalizeOnSave", &DAT::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &DAT::GetNormalizeOnSave)
    .def("SetSigned", &DAT::SetSigned)
    .def("GetSigned", &DAT::GetSigned)
    .def("SetBitDepth", &DAT::SetBitDepth)
    .def("GetBitDepth", &DAT::GetBitDepth)
  ;

  class_<IPL, bases<ImageFormatBase> >("IPL", init<bool,Format>((arg("normalize_on_save") = true,arg("format")=OST_DEFAULT_FORMAT)))
    .def("SetNormalizeOnSave", &IPL::SetNormalizeOnSave)
    .def("GetNormalizeOnSave", &IPL::GetNormalizeOnSave)
    .def("SetBitDepth", &IPL::SetBitDepth)
    .def("GetBitDepth", &IPL::GetBitDepth)
  ;

  class_<JPK, bases<TIF> >("JPK", init<boost::logic::tribool,Format,bool,bool,int>
           ((arg("normalize_on_save") =  boost::logic::tribool(boost::logic::indeterminate),arg("format")=OST_DEFAULT_FORMAT,arg("signed")=false,arg("phasecolor")=false,arg("subimage") = -1)))
  ;

  class_<Nanoscope, bases<ImageFormatBase> >("Nanoscope", init <int>(arg("subimage")=-1))
  ;

  def("SaveImage",save_image1);
  def("SaveImage",save_image2);
  def("LoadImage",load_image1);
  def("LoadImage",load_image2);

}
