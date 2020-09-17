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

/*
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include <cstdio>
#include <sstream>
#include <algorithm>

#include <stdarg.h>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <ost/img/image_state.hh>
#include <ost/img/progress.hh>
#include <ost/units.hh>
#include <ost/img/normalizer_impl.hh>
#include <ost/img/alg/transform.hh>
#include <ost/log.hh>
#include <ost/io/io_exception.hh>

#include "map_io_tiff_handler.hh"
#include "tiff_util.hh"

namespace ost { namespace io {

String TIF::FORMAT_STRING="defined_tiff";

 TIF::TIF(boost::logic::tribool normalize_on_save, Format bit_depth, bool sign, bool phasecolor, int subimage):
  ImageFormatBase(FORMAT_STRING),
  normalize_on_save_(normalize_on_save),
  bit_depth_(bit_depth),
  signed_(sign),
  phasecolor_(phasecolor),
  subimage_(subimage)
{
}

TIF::TIF(String format_string, boost::logic::tribool  normalize_on_save, Format bit_depth,bool sign, bool phasecolor, int subimage):
  ImageFormatBase(FORMAT_STRING),
  normalize_on_save_(normalize_on_save),
  bit_depth_(bit_depth),
  signed_(sign),
  phasecolor_(phasecolor),
  subimage_(subimage)
{
}

Format TIF::GetBitDepth() const
{
  return bit_depth_;
}

void TIF::SetBitDepth (Format bitdepth)
{
  bit_depth_ = bitdepth;
}

bool TIF::GetSigned () const
{
  return signed_;
}

void TIF::SetSigned (bool sign)
{
  signed_ = sign;
}

bool TIF::GetPhasecolor () const
{
  return phasecolor_;
}

void TIF::SetPhasecolor (bool phasecolor)
{
  phasecolor_ = phasecolor;
}

boost::logic::tribool TIF::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void TIF::SetNormalizeOnSave(boost::logic::tribool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}

Format TIF::GetFormat () const
{
  return bit_depth_;
}

void TIF::SetFormat (Format bit_depth)
{
  bit_depth_ = bit_depth;
}

int TIF::GetSubimage () const
{
  return subimage_;
}

void TIF::SetSubimage (int subimage)
{
  subimage_ = subimage;
}

Real TIF::GetMaximum() const
{
  if(GetSigned()){
    switch(GetBitDepth()){
    case OST_BIT8_FORMAT:
      return 127.0;
      break;
    case OST_BIT16_FORMAT:
    case OST_DEFAULT_FORMAT:
      return 32767.0;
      break;
    case OST_BIT32_FORMAT:
      return 2147483647.0;
      break;
    default:
      return 1.0;
      break;
    }
  }else{
    switch(GetBitDepth()){
    case OST_BIT8_FORMAT:
      return 255.0;
      break;
    case OST_BIT16_FORMAT:
    case OST_DEFAULT_FORMAT:
      return 65535.0;
      break;
    case OST_BIT32_FORMAT:
      return 4294967295.0;
      break;
    default:
      return 1.0;
      break;
    }
  }
  return 1.0;
}

Real TIF::GetMinimum() const
{
  if(GetSigned()){
    switch(GetBitDepth()){
    case OST_BIT8_FORMAT:
      return -128.0;
      break;
    case OST_BIT16_FORMAT:
    case OST_DEFAULT_FORMAT:
      return -32768.0;
      break;
    case OST_BIT32_FORMAT:
      return -2147483648.0;
      break;
    default:
      return 0.0;
      break;
    }
  }
  return 0.0;
}

void MapIOTiffHandler::Import(img::MapHandle& image, const boost::filesystem::path& location,const ImageFormatBase& formatstruct)
{
  TIF form;
  TIF& formattif = form;
  if (formatstruct.GetFormatString()==TIF::FORMAT_STRING) {
    formattif = formatstruct.As<TIF>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  LOG_VERBOSE("I/O Tiff: using library version " << TIFFGetVersion());
  detail::tiff_warning_handler_wrapper twhw;
  TIFF* tfile=open_subimage_file(location,formattif);
  load_image_data(tfile,image,formattif);
  TIFFClose(tfile);
}

void MapIOTiffHandler::Import(img::ImageHandle& image, std::istream& location, const ImageFormatBase& formatstruct)
{
  TIF form;
  TIF& formattif = form;
  if (formatstruct.GetFormatString()==TIF::FORMAT_STRING) {
    formattif = formatstruct.As<TIF>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  LOG_VERBOSE("I/O Tiff: using library version " << TIFFGetVersion());
  detail::tiff_warning_handler_wrapper twhw;
  TIFF* tfile=open_subimage_stream(location,formattif);
  load_image_data(tfile,image,formattif);
  TIFFClose(tfile);
}

void MapIOTiffHandler::Export(const img::MapHandle& image , const boost::filesystem::path& location,const ImageFormatBase& formatstruct) const
{
  const String& filename = location.string();
  TIFF* tfile=TIFFOpen(filename.c_str(),"w");
  if(!tfile) {
    throw IOException("could not open " + filename + String(" for writing"));
  }
  TIF form;
  TIF& formattif = form;
  if (formatstruct.GetFormatString()==TIF::FORMAT_STRING) {
    formattif = formatstruct.As<TIF>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  do_export(image,tfile,formattif);
}

void MapIOTiffHandler::Export(const img::MapHandle& image, std::ostream& loc,const ImageFormatBase& formatstruct) const
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
  if (formatstruct.GetFormatString()==TIF::FORMAT_STRING) {
    formattif = formatstruct.As<TIF>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  do_export(image,tfile,formattif);
}

void MapIOTiffHandler::do_export(const img::MapHandle& image,TIFF* tfile,TIF& formattif) const
{
  if (image.IsFrequency()==true)
  {
    formattif.SetSigned(true);
  }

  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if(boost::logic::indeterminate(formattif.GetNormalizeOnSave())){
    switch(formattif.GetBitDepth()){
    case OST_BIT8_FORMAT:
    case OST_BIT16_FORMAT:
    case OST_BIT32_FORMAT:
    case OST_DEFAULT_FORMAT:
      norm = img::alg::CreateLinearRangeNormalizer(image,formattif.GetMinimum(),formattif.GetMaximum());
      LOG_VERBOSE("Autodetecting normalization for export: normalization on");
      break;
    case OST_FLOAT_FORMAT:
    case OST_DOUBLE_FORMAT:
    LOG_VERBOSE("Autodetecting normalization for export: normalization off");
      break;
    }
  } else if (formattif.GetNormalizeOnSave() == true) {
    LOG_VERBOSE("Normalization used for export.");
    norm = img::alg::CreateLinearRangeNormalizer(image,formattif.GetMinimum(),formattif.GetMaximum());
  } else {
    LOG_VERBOSE("No normalization used for export.");
  }

  detail::tiff_warning_handler_wrapper twhw;

  img::image_state::RealSpatialImageState *isr = NULL ;
  img::image_state::ComplexSpatialImageState *isc = NULL;
  img::image_state::WordSpatialImageState *isw = NULL;

  boost::function<void (TIFF *,img::image_state::RealSpatialImageState*,unsigned int,unsigned int,unsigned int,unsigned int,const  img::NormalizerPtr& )> fsr;
  boost::function<void (TIFF *,img::image_state::ComplexSpatialImageState*,unsigned int,unsigned int,unsigned int,unsigned int,const  img::NormalizerPtr&)> fsc;
  boost::function<void (TIFF *,img::image_state::WordSpatialImageState*,unsigned int,unsigned int,unsigned int,unsigned int,const  img::NormalizerPtr&)> fsw;

  unsigned int width=image.GetSize().GetWidth();
  unsigned int height=image.GetSize().GetHeight();
  uint16 spp=1;
  uint16 fmt =0;
  uint16 bpp =0 ;
  img::Point ori=image.GetSpatialOrigin();
  geom::Vec3 sampling=image.GetPixelSampling();
  float xreso=sampling[0]/Units::cm,yreso=sampling[1]/Units::cm;
  float xpos=std::max<Real>(0.0,xreso*ori[0]),ypos=std::max<Real>(0.0,yreso*ori[1]); //tiff file format only allows positivie origins, negative origins are lost here
  TIFFSetField(tfile,TIFFTAG_IMAGEWIDTH,width);
  TIFFSetField(tfile,TIFFTAG_IMAGELENGTH,height);
  TIFFSetField(tfile,TIFFTAG_SAMPLESPERPIXEL,spp);
  TIFFSetField(tfile,TIFFTAG_RESOLUTIONUNIT,RESUNIT_CENTIMETER);
  TIFFSetField(tfile,TIFFTAG_XPOSITION,xpos);
  TIFFSetField(tfile,TIFFTAG_XRESOLUTION,xreso);
  TIFFSetField(tfile,TIFFTAG_YPOSITION,ypos);
  TIFFSetField(tfile,TIFFTAG_YRESOLUTION,yreso);
  TIFFSetField(tfile,TIFFTAG_PLANARCONFIG,1);
  TIFFSetField(tfile,TIFFTAG_ORIENTATION,ORIENTATION_TOPLEFT);
  TIFFSetField(tfile,TIFFTAG_PHOTOMETRIC,1); // 1=zero is black
  if(image.GetDomain()!=img::SPATIAL){
    throw IOException("Frequency export not implemented.");
  }
  if(image.GetType()==img::REAL) {
    isr= dynamic_cast<img::image_state::RealSpatialImageState*>(image.ImageStatePtr().get());
    if(formattif.GetSigned()){ // signed
      switch(formattif.GetBitDepth()){
      case OST_BIT8_FORMAT:
        fmt=SAMPLEFORMAT_INT;
        bpp=8;
        fsr=detail::do_tiff_write<Real,int8,img::image_state::RealSpatialImageState>;
        break;
      case OST_BIT16_FORMAT:
        fmt=SAMPLEFORMAT_INT;
        bpp=16;
        fsr=detail::do_tiff_write<Real,int16,img::image_state::RealSpatialImageState>;
        break;
      case OST_BIT32_FORMAT:
        fmt=SAMPLEFORMAT_INT;
        bpp=32;
        fsr=detail::do_tiff_write<Real,int32,img::image_state::RealSpatialImageState>;
        break;
      case OST_FLOAT_FORMAT:
        fmt=SAMPLEFORMAT_IEEEFP;
        bpp=32;
        fsr=detail::do_tiff_write<Real,float,img::image_state::RealSpatialImageState>;
        break;
      case OST_DOUBLE_FORMAT:
        fmt=SAMPLEFORMAT_IEEEFP;
        bpp=64;
        fsr=detail::do_tiff_write<Real,Real,img::image_state::RealSpatialImageState>;
        break;
      case OST_DEFAULT_FORMAT:
        fmt=SAMPLEFORMAT_INT;
        bpp=16;
        fsr=detail::do_tiff_write<Real,int16,img::image_state::RealSpatialImageState>;
        break;
      }
    }else{ // unsigned
      switch(formattif.GetBitDepth()){
      case OST_BIT8_FORMAT:
        fmt=SAMPLEFORMAT_UINT;
        bpp=8;
        fsr=detail::do_tiff_write<Real,uint8,img::image_state::RealSpatialImageState>;
        break;
      case OST_BIT16_FORMAT:
        fmt=SAMPLEFORMAT_UINT;
        bpp=16;
        fsr=detail::do_tiff_write<Real,uint16,img::image_state::RealSpatialImageState>;
        break;
      case OST_BIT32_FORMAT:
        fmt=SAMPLEFORMAT_UINT;
        bpp=32;
        fsr=detail::do_tiff_write<Real,uint32,img::image_state::RealSpatialImageState>;
        break;
      case OST_FLOAT_FORMAT:
        fmt=SAMPLEFORMAT_IEEEFP;
        bpp=32;
        fsr=detail::do_tiff_write<Real,float,img::image_state::RealSpatialImageState>;
        break;
      case OST_DOUBLE_FORMAT:
        fmt=SAMPLEFORMAT_IEEEFP;
        bpp=64;
        fsr=detail::do_tiff_write<Real,Real,img::image_state::RealSpatialImageState>;
        break;
      case OST_DEFAULT_FORMAT:
        fmt=SAMPLEFORMAT_UINT;
        bpp=16;
        fsr=detail::do_tiff_write<Real,uint16,img::image_state::RealSpatialImageState>;
        break;
      }
    }
  } else if(image.GetType()==img::COMPLEX) {
    isc= dynamic_cast<img::image_state::ComplexSpatialImageState*>(image.ImageStatePtr().get());
    switch(formattif.GetBitDepth()){
    case OST_BIT8_FORMAT:
      fmt=SAMPLEFORMAT_COMPLEXINT;
      bpp=16;
      fsc=detail::do_tiff_write<Complex,std::complex<int8>,img::image_state::ComplexSpatialImageState>;
      break;
    case OST_BIT16_FORMAT:
      fmt=SAMPLEFORMAT_COMPLEXINT;
      bpp=32;
      fsc=detail::do_tiff_write<Complex,std::complex<int16>,img::image_state::ComplexSpatialImageState>;
      break;
    case OST_BIT32_FORMAT:
      fmt=SAMPLEFORMAT_COMPLEXINT;
      bpp=64;
      fsc=detail::do_tiff_write<Complex,std::complex<int32>,img::image_state::ComplexSpatialImageState>;
      break;
    case OST_FLOAT_FORMAT:
      fmt=SAMPLEFORMAT_COMPLEXIEEEFP;
      bpp=64;
      fsc=detail::do_tiff_write<Complex,std::complex<float>,img::image_state::ComplexSpatialImageState>;
      break;
    case OST_DOUBLE_FORMAT:
      fmt=SAMPLEFORMAT_COMPLEXIEEEFP;
      bpp=128;
      fsc=detail::do_tiff_write<Complex,Complex,img::image_state::ComplexSpatialImageState>;
      break;
    case OST_DEFAULT_FORMAT:
      fmt=SAMPLEFORMAT_COMPLEXIEEEFP;
      bpp=128;
      fsc=detail::do_tiff_write<Complex,Complex,img::image_state::ComplexSpatialImageState>;
      break;
    }
  } else if(image.GetType()==img::WORD) {
    isw= dynamic_cast<img::image_state::WordSpatialImageState*>(image.ImageStatePtr().get());
    fsw=detail::do_tiff_write<uint,uint,img::image_state::WordSpatialImageState>;
    fmt=SAMPLEFORMAT_UINT;
    bpp=16;
  }

  TIFFSetField(tfile,TIFFTAG_BITSPERSAMPLE,bpp);
  unsigned int rowsperstrip=std::max<int>(1,8192/(width*bpp/8)); //ca. 8 kb per strip (adobe tiff specs recomendation)
  TIFFSetField(tfile,TIFFTAG_ROWSPERSTRIP,rowsperstrip);
  TIFFSetField(tfile,TIFFTAG_SAMPLEFORMAT,fmt);
  unsigned int stripcount=static_cast<int>(ceil( height/static_cast<Real>(rowsperstrip)));
  img::Progress::Instance().Register(this,stripcount,100);

  if(image.GetType()==img::REAL) {
    for(unsigned int strip=0;strip<stripcount;++strip)
    {
      fsr(tfile,isr,rowsperstrip,width,height,strip,norm.GetImplementation());
      img::Progress::Instance().AdvanceProgress(this);
    }
  } else if(image.GetType()==img::COMPLEX) {
    for(unsigned int strip=0;strip<stripcount;++strip)
    {
      fsc(tfile,isc,rowsperstrip,width,height,strip,norm.GetImplementation());
      img::Progress::Instance().AdvanceProgress(this);
    }
  } else if(image.GetType()==img::WORD) {
    for(unsigned int strip=0;strip<stripcount;++strip)
    {
      fsw(tfile,isw,rowsperstrip,width,height,strip,norm.GetImplementation());
      img::Progress::Instance().AdvanceProgress(this);
    }
  }

  img::Progress::Instance().DeRegister(this);
  TIFFClose(tfile);
}

bool MapIOTiffHandler::MatchContent(unsigned char* header)
{
  if(
     (header[0]==0x4D &&
      header[1]==0x4D &&
      header[2]==0x00 &&
      header[3]==0x2A)
     ||
     (header[0]==0x49 &&
      header[1]==0x49 &&
      header[2]==0x2A &&
      header[3]==0x00)
     ){
      return true;
  }
  return false;
}
bool MapIOTiffHandler::MatchType(const ImageFormatBase& type)
{
  return (type.GetFormatString()=="defined_tiff");
}

bool MapIOTiffHandler::MatchSuffix(const String& loc)
{
  return (detail::FilenameEndsWith(loc,".tif") || detail::FilenameEndsWith(loc,".tiff") );
}


TIFF* MapIOTiffHandler::open_subimage_file(const boost::filesystem::path& location,const TIF& formattif)
{
  unsigned int directory_number=0;
  if(formattif.GetSubimage()!=-1){
    directory_number=formattif.GetSubimage();
  }
  LOG_VERBOSE("I/O Tiff: using directory nr. " << directory_number);
  const String& filename = location.string();
  LOG_VERBOSE("I/O Tiff: calling TIFFOpen with " << filename.c_str());
  TIFF* tfile=TIFFOpen(filename.c_str(),"r");
  if(!tfile) {
    throw IOException("could not open " + filename + String(" for reading"));
  }
  if(directory_number>0) {
    TIFFSetDirectory(tfile,directory_number);
  }
  if(Logger::Instance().GetVerbosityLevel()>=5) {
    TIFFPrintDirectory(tfile,stderr,directory_number);
  }
  return tfile;
}

TIFF* MapIOTiffHandler::open_subimage_stream(std::istream& location,const TIF& formattif)
{
  unsigned int directory_number=0;
  if(formattif.GetSubimage()!=-1){
    directory_number=formattif.GetSubimage();
  }
  LOG_VERBOSE("I/O Tiff: using directory nr. " << directory_number);
  const String& filename = "stream";
  LOG_VERBOSE("I/O Tiff: calling TIFFClientOpen with " << filename.c_str());
  TIFF* tfile=TIFFClientOpen(filename.c_str(),
                             "rm",
                             &location,
                             &detail::CustomTIFFReadProcIStream,
                             &detail::CustomTIFFWriteProcIStream,
                             &detail::CustomTIFFSeekProcIStream,
                             &detail::CustomTIFFCloseProc,
                             &detail::CustomTIFFSizeProcIStream,
                             &detail::CustomTIFFMapFileProc,
                             &detail::CustomTIFFUnmapFileProc);
  if(!tfile) {
    throw IOException("could not open stream for reading");
  }
  if(directory_number>0) {
    TIFFSetDirectory(tfile,directory_number);
  }
  if(Logger::Instance().GetVerbosityLevel()>=5) {
    TIFFPrintDirectory(tfile,stderr,directory_number);
  }
  return tfile;
}

void MapIOTiffHandler::load_image_data(TIFF* tfile, img::ImageHandle& image,  const TIF& formattif)
{
  unsigned int width,height,rps;
  uint16 bpp,spp,plc,ori,reso=RESUNIT_NONE,fmt=SAMPLEFORMAT_UINT;
  float xpos=0.0,ypos=0.0;
  float xreso=1.0,yreso=1.0;

  TIFFGetField(tfile,TIFFTAG_IMAGEWIDTH,&width);
  TIFFGetField(tfile,TIFFTAG_IMAGELENGTH,&height);
  TIFFGetField(tfile,TIFFTAG_BITSPERSAMPLE,&bpp);
  if(!TIFFGetField(tfile,TIFFTAG_SAMPLESPERPIXEL,&spp)) {
    spp=0;
  }
  TIFFGetField(tfile,TIFFTAG_PLANARCONFIG,&plc);
  if(!TIFFGetField(tfile,TIFFTAG_ORIENTATION,&ori)) {
    ori=0;
  }
  TIFFGetField(tfile,TIFFTAG_ROWSPERSTRIP,&rps);

  LOG_INFO("I/O Tiff: Header: " << width << "x" << height << " " << bpp << "bpp");

  img::Extent image_extent(img::Point(0,0),img::Size(width,height));

  if(!TIFFGetField(tfile,TIFFTAG_RESOLUTIONUNIT,&reso)){
    reso=RESUNIT_NONE;
  }
  if(!(TIFFGetField(tfile,TIFFTAG_XPOSITION,&xpos) &&
       TIFFGetField(tfile,TIFFTAG_YPOSITION,&ypos))){
        xpos=0.0;
        ypos=0.0;
  }
  if(!(TIFFGetField(tfile,TIFFTAG_XRESOLUTION,&xreso) &&
       TIFFGetField(tfile,TIFFTAG_YRESOLUTION,&yreso))){
    xreso=1.0;
    yreso=1.0;
  }
  xreso = std::abs(xreso)<1e-20 ? 1.0 : xreso;
  yreso = std::abs(yreso)<1e-20 ? 1.0 : yreso;

  if(!TIFFGetField(tfile,TIFFTAG_SAMPLEFORMAT,&fmt))
  {
    // todo check subtype
    LOG_INFO("I/O Tiff: Unknown sample format: assuming uint");
    fmt=SAMPLEFORMAT_UINT;
  }

  LOG_INFO("I/O Tiff: Header: spp=" << spp << "; rps=" << rps << "; plc=" << plc << "; ori=" << ori);
  LOG_INFO("I/O Tiff: Header: pos= " << xpos << "," << ypos << "; reso= " << xreso << "," << yreso << ";");

  if(bpp!=8 && bpp!=16 && bpp!=32 && bpp!=64 && bpp!=128) {
    throw IOException("bits per sample must be 8, 16, 32, 64 or 128");
  }

  if(spp==0) {
    LOG_INFO("I/O Tiff: Samples per pixel is zero, assuming one");
    spp=1;
  }

  if(spp!=1) {
    if(spp>4) {
      throw IOException("I/O Tiff: unexpected spp>4 in header");
    }
    LOG_INFO("I/O Tiff: Samples per pixel>1, reading first channel only");
  }

  if(TIFFIsTiled(tfile)) {
    throw IOException("tiled tiff image format not implemented");
  }

  tdata_t buf;
  tstrip_t strip;
  uint stripsize=TIFFStripSize(tfile);
  uint stripcount=TIFFNumberOfStrips(tfile);

  LOG_INFO("I/O Tiff: " << "found " << stripcount << " strips in tiff file");
  LOG_INFO("I/O Tiff: " << "allocating " << stripsize << " bytes for temporary read buffer");
  buf = _TIFFmalloc(stripsize);

  if(image.GetType()==img::WORD) {
    LOG_INFO("I/O Tiff: " << "resetting target image to WORD " << image_extent);
    image.Reset(image_extent, img::WORD, img::SPATIAL);

    img::image_state::WordSpatialImageState *is =
        dynamic_cast<img::image_state::WordSpatialImageState*>(image.ImageStatePtr().get());
    if(!is) {
      throw IOException("unexpected failure of dynamic_cast in tiff io");
    }

    unsigned int current_row = 0;

    LOG_VERBOSE("I/O Tiff: " << "importing data");
    img::Progress::Instance().Register(this,stripcount,100);
    for (strip = 0; strip < stripcount; ++strip) {
      int cread=TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) -1);

      if(bpp==8) {
        int rowcount = cread/(width*spp);
        detail::do_tiff_read<uchar,short,img::image_state::WordSpatialImageState>(buf,rowcount,width,is,current_row,spp);
      } else if(bpp==16) {
        uint rowcount = cread/(width*2*spp);
        detail::do_tiff_read<uint16,short,img::image_state::WordSpatialImageState>(buf,rowcount,width,is,current_row,spp);
      } else if(bpp==32) {
        uint rowcount = cread/(width*4*spp);
        detail::do_tiff_read<uint32,short,img::image_state::WordSpatialImageState>(buf,rowcount,width,is,current_row,spp);
      }
      img::Progress::Instance().AdvanceProgress(this);
    }
    LOG_VERBOSE("I/O Tiff: " << "done");
    img::Progress::Instance().DeRegister(this);
    uint16 photometric;
    if(TIFFGetField(tfile,TIFFTAG_PHOTOMETRIC,&photometric)){
      if(photometric==PHOTOMETRIC_PALETTE ){
        uint16* lookuptable_red;
        uint16* lookuptable_green;
        uint16* lookuptable_blue;
        TIFFGetField(tfile,TIFFTAG_COLORMAP,&lookuptable_red,&lookuptable_green,&lookuptable_blue);
        for(unsigned short* ptr = is->Data().GetData(); ptr<is->Data().GetEnd(); ++ptr) {
          long idx=static_cast<long>(*ptr);
          *ptr = static_cast<unsigned short>((lookuptable_red[idx]+lookuptable_green[idx]+lookuptable_blue[idx])/3.0);
        }
      }
    }
  } else {
    img::image_state::RealSpatialImageState *isr=0;
    img::image_state::ComplexSpatialImageState *isc=0;
    if(fmt==SAMPLEFORMAT_COMPLEXINT || fmt==SAMPLEFORMAT_COMPLEXIEEEFP){
      LOG_INFO("I/O Tiff: " << "resetting target image to complex spatial " << image_extent);
      image.Reset(image_extent,img::COMPLEX,img::SPATIAL);
      isc = dynamic_cast<img::image_state::ComplexSpatialImageState*>(image.ImageStatePtr().get());
      if(!isc) {
        throw IOException("unexpected failure of dynamic_cast in tiff io");
      }
    }else{
      LOG_INFO("I/O Tiff: " << "resetting target image to Real img::SPATIAL" << image_extent);
      image.Reset(image_extent, img::REAL, img::SPATIAL);

      isr= dynamic_cast<img::image_state::RealSpatialImageState*>(image.ImageStatePtr().get());
      if(!isr) {
        throw IOException("unexpected failure of dynamic_cast in tiff io");
      }
    }


    unsigned int current_row = 0;

    LOG_VERBOSE("I/O Tiff: " << "importing data");
    img::Progress::Instance().Register(this,stripcount,100);
    for (strip = 0; strip < stripcount; ++strip) {
      int cread=TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) -1);

      if(bpp==8) {
        int rowcount = cread/(width*spp);
        if(fmt==SAMPLEFORMAT_INT){
          detail::do_tiff_read<int8,Real,img::image_state::RealSpatialImageState>(buf,rowcount,width,isr,current_row,spp);
        }else if(fmt==SAMPLEFORMAT_UINT){
          detail::do_tiff_read<uint8,Real,img::image_state::RealSpatialImageState>(buf,rowcount,width,isr,current_row,spp);
        }else{
          std::ostringstream ostrs;
          ostrs << "TIFF images with "<<bpp << " bits per pixel and sample format "<< fmt<< " are not supported.";
          throw IOException(ostrs.str());
        }
      } else if(bpp==16) {
        uint rowcount = cread/(width*2*spp);
        if(fmt==SAMPLEFORMAT_INT){
          detail::do_tiff_read<int16,Real,img::image_state::RealSpatialImageState>(buf,rowcount,width,isr,current_row,spp);
        }else if(fmt==SAMPLEFORMAT_UINT){
          detail::do_tiff_read<uint16,Real,img::image_state::RealSpatialImageState>(buf,rowcount,width,isr,current_row,spp);
        }else if(fmt==SAMPLEFORMAT_COMPLEXINT){
          detail::do_tiff_read<detail::complexint8,Complex,img::image_state::ComplexSpatialImageState>(buf,rowcount,width,isc,current_row,spp);
        }else{
          std::ostringstream ostrs;
          ostrs << "TIFF images with "<<bpp << " bits per pixel and sample format "<< fmt<< " are not supported.";
          throw IOException(ostrs.str());
        }
      } else if(bpp==32) {
        uint rowcount = cread/(width*4*spp);
        if(fmt==SAMPLEFORMAT_INT){
          detail::do_tiff_read<int32,Real,img::image_state::RealSpatialImageState>(buf,rowcount,width,isr,current_row,spp);
        }else if(fmt==SAMPLEFORMAT_UINT){
          detail::do_tiff_read<uint32,Real,img::image_state::RealSpatialImageState>(buf,rowcount,width,isr,current_row,spp);
        }else if(fmt==SAMPLEFORMAT_IEEEFP){
          detail::do_tiff_read<float,Real,img::image_state::RealSpatialImageState>(buf,rowcount,width,isr,current_row,spp);
        }else if(fmt==SAMPLEFORMAT_COMPLEXINT){
          detail::do_tiff_read<detail::complexint16,Complex,img::image_state::ComplexSpatialImageState>(buf,rowcount,width,isc,current_row,spp);
        }else{
          std::ostringstream ostrs;
          ostrs << "TIFF images with "<<bpp << " bits per pixel and sample format "<< fmt<< " are not supported.";
          throw IOException(ostrs.str());
        }
      } else if(bpp==64) {
        uint rowcount = cread/(width*8*spp);
        if(fmt==SAMPLEFORMAT_IEEEFP){
          detail::do_tiff_read<Real,Real,img::image_state::RealSpatialImageState>(buf,rowcount,width,isr,current_row,spp);
        }else if(fmt==SAMPLEFORMAT_COMPLEXINT){
          detail::do_tiff_read<detail::complexint32,Complex,img::image_state::ComplexSpatialImageState>(buf,rowcount,width,isc,current_row,spp);
        }else if(fmt==SAMPLEFORMAT_COMPLEXIEEEFP){
          detail::do_tiff_read<std::complex<float>,Complex,img::image_state::ComplexSpatialImageState>(buf,rowcount,width,isc,current_row,spp);
        }else{
          std::ostringstream ostrs;
          ostrs << "TIFF images with "<<bpp << " bits per pixel and sample format "<< fmt<< " are not supported.";
          throw IOException(ostrs.str());
        }
      } else if(bpp==128) {
        uint rowcount = cread/(width*16*spp);
        if(fmt==SAMPLEFORMAT_COMPLEXIEEEFP){
          // std::complex<double> refers to the dataformat used in the tiff file and is independend of single/double precision compilation
          detail::do_tiff_read<std::complex<double>,Complex,img::image_state::ComplexSpatialImageState>(buf,rowcount,width,isc,current_row,spp);
        }else{
          std::ostringstream ostrs;
          ostrs << "TIFF images with "<<bpp << " bits per pixel and sample format "<< fmt<< " are not supported.";
          throw IOException(ostrs.str());
        }
      }
      img::Progress::Instance().AdvanceProgress(this);
    }
    LOG_VERBOSE("I/O Tiff: " << "done");

    img::Progress::Instance().DeRegister(this);
    uint16 photometric;
    if(TIFFGetField(tfile,TIFFTAG_PHOTOMETRIC,&photometric)){
      if(photometric==PHOTOMETRIC_PALETTE && image.GetType()==img::REAL){
        uint16* lookuptable_red;
        uint16* lookuptable_green;
        uint16* lookuptable_blue;
        TIFFGetField(tfile,TIFFTAG_COLORMAP,&lookuptable_red,&lookuptable_green,&lookuptable_blue);
        for(Real* ptr = isr->Data().GetData(); ptr<isr->Data().GetEnd(); ++ptr) {
          long idx=static_cast<long>(*ptr);
          *ptr = (lookuptable_red[idx]+lookuptable_green[idx]+lookuptable_blue[idx])/3.0;
        }
      }
    }
  }
  switch(ori){
  case ORIENTATION_TOPRIGHT:
    image.ApplyIP(img::alg::Mirror(img::Plane::YZ));
    break;
  case ORIENTATION_BOTRIGHT:

    image.ApplyIP(img::alg::Transform(img::alg::Transformation(geom::Mat4(-1,0,0,image.GetExtent().GetSize()[0]-1,
                                                                                 0,-1,0,image.GetExtent().GetSize()[1]-1,
                                                                                 0,0,1,0,
                                                                                 0,0,0,1))));
    break;
  case ORIENTATION_BOTLEFT:
    image.ApplyIP(img::alg::Mirror(img::Plane::XZ));
    break;
  case ORIENTATION_LEFTTOP:{
    int maxsize=std::max<int>(image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[1]);
    img::ImageHandle tmp=CreateImage(img::Extent(img::Size(maxsize,maxsize,image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
    tmp.Paste(image);
    tmp.ApplyIP(img::alg::Transform(img::alg::Transformation(geom::Mat4(0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1))));
    image.Reset(img::Extent(img::Size(image.GetExtent().GetSize()[1],image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
    image.Paste(tmp);
  }
    break;
  case ORIENTATION_RIGHTTOP:{
    int maxsize=std::max<int>(image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[1]);
    img::ImageHandle tmp=CreateImage(img::Extent(img::Size(maxsize,maxsize,image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
    tmp.Paste(image);
    tmp.ApplyIP(img::alg::Transform(ost::img::alg::Transformation(geom::Mat4(0,-1,0,image.GetExtent().GetSize()[1]-1,
                                                                                 1,0,0,0,
                                                                                 0,0,1,0,
                                                                                 0,0,0,1))));
    image.Reset(img::Extent(img::Size(image.GetExtent().GetSize()[1],image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
    image.Paste(tmp);
  }
    break;
  case ORIENTATION_RIGHTBOT:{
    int maxsize=std::max<int>(image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[1]);
    img::ImageHandle tmp=CreateImage(img::Extent(img::Size(maxsize,maxsize,image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
    tmp.Paste(image);
    tmp.ApplyIP(img::alg::Transform(img::alg::Transformation(geom::Mat4(0,-1,0,image.GetExtent().GetSize()[1]-1 ,-1,0,0,image.GetExtent().GetSize()[0]-1 ,0,0,1,0,  0,0,0,1))));
    image.Reset(img::Extent(img::Size(image.GetExtent().GetSize()[1],image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
    image.Paste(tmp);
  }
    break;
  case ORIENTATION_LEFTBOT:{
    int maxsize=std::max<int>(image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[1]);
    img::ImageHandle tmp=CreateImage(img::Extent(img::Size(maxsize,maxsize,image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
    tmp.Paste(image);
    tmp.ApplyIP(img::alg::Transform(img::alg::Transformation(geom::Mat4(0,1,0,0,
                                                                                 -1,0,0,image.GetExtent().GetSize()[0]-1,
                                                                                 0,0,1,0,
                                                                                 0,0,0,1))));
    image.Reset(img::Extent(img::Size(image.GetExtent().GetSize()[1],image.GetExtent().GetSize()[0],image.GetExtent().GetSize()[2])),image.GetType(),image.GetDomain());
    image.Paste(tmp);
  }
    break;
  case ORIENTATION_TOPLEFT: //fallthrough intended
  default:
    break;
  }
  switch(reso){
  case RESUNIT_NONE:
    break;
  case RESUNIT_INCH:
    image.SetPixelSampling(geom::Vec3(xreso*Units::in,yreso*Units::in,1));
    image.SetSpatialOrigin(img::Point(static_cast<int>(round(xpos/xreso)),static_cast<int>(round(ypos/yreso)),0));
    break;
  case RESUNIT_CENTIMETER:
    image.SetPixelSampling(geom::Vec3(xreso*Units::cm,yreso*Units::cm,1));
    image.SetSpatialOrigin(img::Point(static_cast<int>(round(xpos/xreso)),static_cast<int>(round(ypos/yreso)),0));
    break;
  }
  _TIFFfree(buf);
}

}} // ns
