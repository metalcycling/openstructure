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

#include <limits>
#include <string>
#include <ost/units.hh>
#include <ost/img/alg/transform.hh>
#include <ost/units.hh>
#include <ost/img/alg/transform.hh>
#include <ost/img/image_handle.hh>
#include <ost/img/progress.hh>
#include <ost/geom/vec3.hh>
#include <ost/geom/mat4.hh>
#include <ost/log.hh>
#include <ost/io/io_exception.hh>

#include <boost/function.hpp>

#include "jpk_util.hh"
#include "tiff_util.hh"
#include "map_io_jpk_handler.hh"

namespace ost { namespace io {

String JPK::FORMAT_STRING = "defined_jpk";

JPK::JPK(boost::logic::tribool normalize_on_save, Format bit_depth, bool sign, bool phasecolor, int subimage):
  TIF(FORMAT_STRING,normalize_on_save,bit_depth,sign,phasecolor,subimage)
{}


void MapIOJpkHandler::Import(img::ImageHandle& image, const boost::filesystem::path& location, const ImageFormatBase& formatstruct)
{

  JPK form;
  JPK& formatjpk = form;
  if (formatstruct.GetFormatString()==JPK::FORMAT_STRING) {
    formatjpk = formatstruct.As<JPK>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  String tmptype;
  JPK newformat=formatjpk;
  if(formatjpk.GetSubimage()==-1){
    newformat.SetSubimage(1);
  }

  TIFF* tfile=open_subimage_file(location.string(),newformat);
  load_image_data(tfile,image,newformat);
  int count;
  Real* scale;
  Real* base;
  int* num_slots;
  char* def_slot_name;
  char* slot_name;

  if (TIFFGetField(tfile, jpk_tags::SLOT_COUNT, &count, &num_slots)) {
    if (TIFFGetField(tfile, jpk_tags::DEFAULT_SLOT, &count, &def_slot_name)) {
      for (int i=0;i<(*num_slots);++i) {
        if (TIFFGetField(tfile, jpk_tags::SlotName(i), &count, &slot_name)){
          if(String(def_slot_name).compare(String(slot_name)) == 0){
            if (TIFFGetField(tfile, jpk_tags::ScalingMultiply(i), &count, &scale)
                &&TIFFGetField(tfile, jpk_tags::ScalingOffset(i), &count, &base)) {
              Real units_scale=1.0;
              char *encunit;
              if (TIFFGetField(tfile, jpk_tags::EncoderUnit(i), &count, &encunit)){
                if(String(encunit).compare("m") == 0){
                  units_scale=Units::m;
                }
              }
              image*=(*scale*units_scale);
              image+=(*base*units_scale);
              LOG_INFO("jpk_plugin: scalefactor: " << *scale 
                           << "  offset: " << *base);
            } else {
              throw IOException("Missing scaling tags");
            }
          }
        }else {
          throw IOException("Missing slot name tag");
        }
      }
    } else {
      throw IOException("Missing default slot name tag");
    }
  } else {
    throw IOException("Missing slot tag");
  }
  char* s;
  TIFFGetField(tfile, jpk_tags::EncoderUnit(newformat.GetSubimage()), &count, &s);
  LOG_INFO("jpk_plugin: units: " << s);


  Real* sx;
  Real* sy;
  TIFFSetDirectory(tfile,0);
  if (TIFFGetField(tfile, jpk_tags::GRID_U_LENGTH, &count, &sx)
    &&TIFFGetField(tfile, jpk_tags::GRID_V_LENGTH, &count, &sy)) {
      image.SetPixelSampling(geom::Vec3(*sx/static_cast<Real>(image.GetSize()[0])*Units::m,
                                        *sy/static_cast<Real>(image.GetSize()[1])*Units::m,
                                        1));
    LOG_INFO("jpk_plugin: sx: " << *sx << "  sy: "<< *sy);
  }
  else {
    throw IOException("Missing size tags");
  }
  // mirror image in case of origin bottom left
  // code duplicated from io_tiff because jpk format ignores standard tiff tag for origin
  uint16_t* reflect;
  uint32_t* dummy;
 if(TIFFGetField(tfile,jpk_tags::GRID_REFLECT,&dummy,&reflect)){
    if(reflect){
      int maxsize=std::max<int>(image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[1]);
      img::ImageHandle tmp=CreateImage(img::Extent(img::Size(maxsize,maxsize,image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
      tmp.Paste(image);

      tmp.ApplyIP(img::alg::Transform(img::alg::Transformation(geom::Mat4(1,0,0,0,
                                                                                   0,-1,0,image.GetExtent().GetSize()[0]-1,
                                                                                   0,0,1,0,
                                                                                   0,0,0,1))));
      image.Reset(img::Extent(img::Size(image.GetExtent().GetSize()[1],image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
      image.Paste(tmp);
    }
  }
  TIFFClose(tfile);
 }

void MapIOJpkHandler::Import(img::MapHandle& image, std::istream& loc,const ImageFormatBase& formatstruct)
{
  JPK form;
  JPK& formatjpk = form;
  if (formatstruct.GetFormatString()==JPK::FORMAT_STRING) {
    formatjpk = formatstruct.As<JPK>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }
  String tmptype;
  JPK newformat=formatjpk;
  if(formatjpk.GetSubimage()==-1){
    newformat.SetSubimage(1);
  }

  TIFF* tfile=open_subimage_stream(loc,newformat);
  load_image_data(tfile,image,newformat);
  int count;
  Real* scale;
  Real* base;
  int* num_slots;
  char* def_slot_name;
  char* slot_name;

  if (TIFFGetField(tfile, jpk_tags::SLOT_COUNT, &count, &num_slots)) {
    if (TIFFGetField(tfile, jpk_tags::DEFAULT_SLOT, &count, &def_slot_name)) {
      for (int i=0;i<(*num_slots);++i) {
        if (TIFFGetField(tfile, jpk_tags::SlotName(i), &count, &slot_name)){
          if(String(def_slot_name).compare(String(slot_name)) == 0){
            if (TIFFGetField(tfile, jpk_tags::ScalingMultiply(i), &count, &scale)
                &&TIFFGetField(tfile, jpk_tags::ScalingOffset(i), &count, &base)) {
              Real units_scale=1.0;
              char *encunit;
              if (TIFFGetField(tfile, jpk_tags::EncoderUnit(i), &count, &encunit)){
                if(String(encunit).compare("m") == 0){
                  units_scale=Units::m;
                }
              }
              image*=(*scale*units_scale);
              image+=(*base*units_scale);
              LOG_INFO("jpk_plugin: scalefactor: " << *scale << "  offset: "
                            << *base);
            }else {
              throw IOException("Missing scaling tags");
            }
          }
        }else {
          throw IOException("Missing slot name tag");
        }
      }
    } else {
      throw IOException("Missing default slot name tag");
    }
  } else {
    throw IOException("Missing slot tag");
  }
  char* s;
  TIFFGetField(tfile, jpk_tags::EncoderUnit(newformat.GetSubimage()), &count, &s);
  LOG_INFO("jpk_plugin: units: " << s);


  Real* sx;
  Real* sy;
  TIFFSetDirectory(tfile,0);
  if (TIFFGetField(tfile, jpk_tags::GRID_U_LENGTH, &count, &sx)
    &&TIFFGetField(tfile, jpk_tags::GRID_V_LENGTH, &count, &sy)) {
      image.SetPixelSampling(geom::Vec3(*sx/static_cast<Real>(image.GetSize()[0])*Units::m,
                                        *sy/static_cast<Real>(image.GetSize()[1])*Units::m,
                                        1));
    LOG_INFO("jpk_plugin: sx: " << *sx << "  sy: "<< *sy);
  }
  else {
    throw IOException("Missing size tags");
  }
  // mirror image in case of origin bottom left
  // code duplicated from io_tiff because jpk format ignores standard tiff tag for origin
  uint16_t* reflect;
  uint32_t* dummy;
 if(TIFFGetField(tfile,jpk_tags::GRID_REFLECT,&dummy,&reflect)){
    if(reflect){
      int maxsize=std::max<int>(image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[1]);
      img::ImageHandle tmp=img::CreateImage(img::Extent(img::Size(maxsize,maxsize,image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
      tmp.Paste(image);

      tmp.ApplyIP(img::alg::Transform(img::alg::Transformation(geom::Mat4(1,0,0,0,
                                                                                   0,-1,0,image.GetExtent().GetSize()[0]-1,
                                                                                   0,0,1,0,
                                                                                   0,0,0,1))));
      image.Reset(img::Extent(img::Size(image.GetExtent().GetSize()[1],image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
      image.Paste(tmp);
    }
  }
  TIFFClose(tfile);
 }

void MapIOJpkHandler::Export(const img::MapHandle& image , const boost::filesystem::path& location,const ImageFormatBase& formatstruct) const
{
  const String& filename = location.string();
  TIFF* tfile=TIFFOpen(filename.c_str(),"w");
  if(!tfile) {
    throw IOException("could not open " + filename + String(" for writing"));
  }

  TIF form;
  TIF& formattif = form;
  if (formatstruct.GetFormatString()==JPK::FORMAT_STRING) {
    formattif = formatstruct.As<TIF>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  do_export(image,tfile,formattif);
}

void MapIOJpkHandler::Export(const img::MapHandle& image, std::ostream& loc,const ImageFormatBase& formatstruct) const
{
  const String& filename = "stream";
  TIFF* tfile=TIFFClientOpen(filename.c_str(),
                             "wm",
                             &loc,
                             &detail::CustomTIFFReadProcOStream,
                             &detail::CustomTIFFWriteProcOStream,
                             &detail::CustomTIFFSeekProcOStream,
                             &detail::CustomTIFFCloseProc,
                             &detail::CustomTIFFSizeProcOStream,
                             &detail::CustomTIFFMapFileProc,
                             &detail::CustomTIFFUnmapFileProc);
  if(!tfile) {
    throw IOException("could not open " + filename + String(" for writing"));
  }

  TIF form;
  TIF& formattif = form;
  if (formatstruct.GetFormatString()==JPK::FORMAT_STRING) {
    formattif = formatstruct.As<TIF>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  do_export(image,tfile,formattif);
}

bool MapIOJpkHandler::MatchContent(unsigned char* header)
{
  // avoid confusion with actual tiff file
  return false;
}

bool MapIOJpkHandler::MatchType(const ImageFormatBase& formatstruct)
{
 if(formatstruct.GetFormatString()=="defined_jpk") {
    return true;
  }
  return false;
}

bool MapIOJpkHandler::MatchSuffix(const String& loc)
{
    if(detail::FilenameEndsWith(loc,".jpk") ) {
      return true;
    }
    return false;
}


}} // ns



