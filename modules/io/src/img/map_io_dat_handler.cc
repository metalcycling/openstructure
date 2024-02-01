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
#include <ost/stdint.hh>
#include <boost/shared_array.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <ost/img/image.hh>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/img/progress.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/convert.hh>
#include <ost/io/converting_streams.hh>
#include <cassert>

#include "map_io_dat_handler.hh"

namespace ost { namespace io {

String DAT::FORMAT_STRING = "defined_dat";

DAT::DAT(bool normalize_on_save, Format bit_depth, bool sign, Endianess end):
    ImageFormatBase(FORMAT_STRING),
    normalize_on_save_(normalize_on_save),
    bit_depth_(bit_depth),
    signed_(sign),
    endianness_(end)
{
}

Format DAT::GetBitDepth() const
{
  return bit_depth_;
}

void DAT::SetBitDepth (Format bitdepth)
{
  bit_depth_ = bitdepth;
}

bool DAT::GetSigned () const
{
  return signed_;
}

void DAT::SetSigned (bool sign)
{
  signed_ = sign;
}

Endianess DAT::GetEndianess() const
{
  return endianness_;
}

void DAT::SetEndianess(Endianess end)
{
  endianness_ = end;
}

bool DAT::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void DAT::SetNormalizeOnSave(bool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}

Real DAT::GetMaximum() const
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

Real DAT::GetMinimum() const
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

bool MapIODatHandler::MatchContent(unsigned char* header)
{
  return false;
}

bool MapIODatHandler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()==DAT::FORMAT_STRING) {
    return true;
  }
  return false;
}

bool MapIODatHandler::MatchSuffix(const String& loc)
{
    if(detail::FilenameEndsWith(loc,".dat") || detail::FilenameEndsWith(loc,".img") ) {
      return true;
    }
    return false;
}

void MapIODatHandler::Import(img::MapHandle& sh, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct )
{
  boost::filesystem::ifstream infile(loc, std::ios::binary);
  if(!infile) {
    throw IOException("could not open "+loc.string());
  }
  this->Import(sh,infile,formatstruct);
  infile.close();
}

template <int CONVERSIONTYPE,typename DATATYPE>
void real_filler(img::MapHandle& sh, std::istream& file)
{
  BinaryIStream<CONVERSIONTYPE> file_bin(file);
  img::Size size = sh.GetSize();
  char this_dummy; //create dummy variable to give to img::Progress as this
  img::Progress::Instance().Register(&this_dummy,size[1],100);
  for(unsigned int row=0;row<size[1];row++) {
    for(unsigned int column=0;column<size[0];column++) {
      DATATYPE value;
      file_bin >> value;
      sh.SetReal(img::Point(column,row),static_cast<Real>(value));
    }
    img::Progress::Instance().AdvanceProgress(&this_dummy);
  }
  img::Progress::Instance().DeRegister(&this_dummy);
}

template <int CONVERSIONTYPE,typename DATATYPE>
void real_dumper(const img::MapHandle& sh, std::ostream& file, const DAT& formatdat)
{
  BinaryOStream<CONVERSIONTYPE> file_bin(file);
  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if (formatdat.GetNormalizeOnSave() == true) {
    norm = img::alg::CreateLinearRangeNormalizer(sh,formatdat.GetMinimum(),formatdat.GetMaximum());
  }
  img::Size size = sh.GetSize();
  char this_dummy; //create dummy variable to give to img::Progress as this
  img::Progress::Instance().Register(&this_dummy,size[1],100);
    for(unsigned int row=0;row<size[1];row++) {
      for(unsigned int column=0;column<size[0];column++)
      {
        file_bin << static_cast<DATATYPE>(norm.Convert(sh.GetReal(img::Point(column,row))));
      }
      img::Progress::Instance().AdvanceProgress(&this_dummy);
    }
  img::Progress::Instance().DeRegister(&this_dummy);
}


void MapIODatHandler::Import(img::MapHandle& sh, std::istream& file, const ImageFormatBase& formatstruct)
{

  DAT form;
  DAT& formatdat = form;
  if (formatstruct.GetFormatString()==DAT::FORMAT_STRING) {
    formatdat = formatstruct.As<DAT>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  uint curr_pos = file.tellg();
  file.seekg(0,std::ios::beg);
  uint start = file.tellg();
  file.seekg(0,std::ios::end);
  uint end = file.tellg();
  uint  filesize = end-start;
  file.seekg(curr_pos,std::ios::beg);

  double datasize=1.0;
  switch(formatdat.GetBitDepth()){
  case OST_BIT8_FORMAT:
  case OST_DEFAULT_FORMAT:
    datasize=1.0;
    break;
  case OST_BIT16_FORMAT:
    datasize=2.0;
    break;
  case OST_BIT32_FORMAT:
    datasize=4.0;
    break;
  case OST_FLOAT_FORMAT:
    datasize=sizeof(float);
    break;
  case OST_DOUBLE_FORMAT:
    datasize=sizeof(double);
    break;
  }
  unsigned int size=static_cast<unsigned int>(sqrt(filesize/datasize));
   // default origin is top-left
  sh.Reset(img::Extent(img::Point(0,0),img::Size(size,size)), img::REAL, img::SPATIAL);

  if(formatdat.GetSigned()){
    switch(formatdat.GetBitDepth()){
    case OST_BIT8_FORMAT:
    case OST_DEFAULT_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,int8_t>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,int8_t>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,int8_t>(sh,file);
        break;
      }
      break;
    case OST_BIT16_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,int16_t>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,int16_t>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,int16_t>(sh,file);
        break;
      }
      break;
    case OST_BIT32_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,int32_t>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,int32_t>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,int32_t>(sh,file);
        break;
      }
      break;
    case OST_FLOAT_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,float>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,float>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,float>(sh,file);
        break;
      }
      break;
    case OST_DOUBLE_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,double>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,double>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,double>(sh,file);
        break;
      }
      break;
    }
  }else{
    switch(formatdat.GetBitDepth()){
    case OST_BIT8_FORMAT:
    case OST_DEFAULT_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,uint8_t>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,uint8_t>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,uint8_t>(sh,file);
        break;
      }
      break;
    case OST_BIT16_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,uint16_t>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,uint16_t>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,uint16_t>(sh,file);
        break;
      }
      break;
    case OST_BIT32_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,uint32_t>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,uint32_t>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,uint32_t>(sh,file);
        break;
      }
      break;
    case OST_FLOAT_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,float>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,float>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,float>(sh,file);
        break;
      }
      break;
    case OST_DOUBLE_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_filler<OST_BIG_ENDIAN,double>(sh,file);
        break;
      case OST_LITTLE_ENDIAN:
        real_filler<OST_LITTLE_ENDIAN,double>(sh,file);
        break;
      case OST_VAX_DATA:
        real_filler<OST_VAX_DATA,double>(sh,file);
        break;
      }
      break;
    }
  }
}

void MapIODatHandler::Export(const img::MapHandle& mh2,
                                  const boost::filesystem::path& loc,const ImageFormatBase& formatstruct) const
{
  boost::filesystem::ofstream outfile(loc, std::ios::binary);
  if(!outfile)
  {
    throw IOException("could not open "+loc.string());
  }
  this->Export(mh2,outfile,formatstruct);
  outfile.close();
}

void MapIODatHandler::Export(const img::MapHandle& sh, std::ostream& file,const ImageFormatBase& formatstruct) const
{

  DAT form;
  DAT& formatdat = form;
  if (formatstruct.GetFormatString()==DAT::FORMAT_STRING) {
    formatdat = formatstruct.As<DAT>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  img::Size size = sh.GetSize();

  if (size[2]!=1 || (size[0]!=size[1])) {
    throw IOException("DAT IO: Dat format only supports square 2D images");
  }

  if(formatdat.GetSigned()){
    switch(formatdat.GetBitDepth()){
    case OST_BIT8_FORMAT:
    case OST_DEFAULT_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,int8_t>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,int8_t>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,int8_t>(sh,file,formatdat);
        break;
      }
      break;
    case OST_BIT16_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,int16_t>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,int16_t>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,int16_t>(sh,file,formatdat);
        break;
      }
      break;
    case OST_BIT32_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,int32_t>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,int32_t>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,int32_t>(sh,file,formatdat);
        break;
      }
      break;
    case OST_FLOAT_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,float>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,float>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,float>(sh,file,formatdat);
        break;
      }
      break;
    case OST_DOUBLE_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,double>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,double>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,double>(sh,file,formatdat);
        break;
      }
      break;
    }
  }else{
    switch(formatdat.GetBitDepth()){
    case OST_BIT8_FORMAT:
    case OST_DEFAULT_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,uint8_t>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,uint8_t>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,uint8_t>(sh,file,formatdat);
        break;
      }
      break;
    case OST_BIT16_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,uint16_t>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,uint16_t>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,uint16_t>(sh,file,formatdat);
        break;
      }
      break;
    case OST_BIT32_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,uint32_t>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,uint32_t>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,uint32_t>(sh,file,formatdat);
        break;
      }
      break;
    case OST_FLOAT_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,float>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,float>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,float>(sh,file,formatdat);
        break;
      }
      break;
    case OST_DOUBLE_FORMAT:
      switch(formatdat.GetEndianess()) {
      case OST_BIG_ENDIAN:
        real_dumper<OST_BIG_ENDIAN,double>(sh,file,formatdat);
        break;
      case OST_LITTLE_ENDIAN:
        real_dumper<OST_LITTLE_ENDIAN,double>(sh,file,formatdat);
        break;
      case OST_VAX_DATA:
        real_dumper<OST_VAX_DATA,double>(sh,file,formatdat);
        break;
      }
      break;
    }
  }
}

}} // namespaces


