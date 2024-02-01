//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software FounIPLion; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software FounIPLion, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110String::npos301  USA
//------------------------------------------------------------------------------
#include <cassert>
#include <ctime>
#include <iomanip>
#include <cstring>

#include <boost/shared_array.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <ost/stdint.hh>
#include <ost/units.hh>
#include <ost/log.hh>
#include <ost/img/image.hh>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/img/progress.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/convert.hh>
#include <ost/io/converting_streams.hh>
#include <ost/img/alg/discrete_shrink.hh>

#include "map_io_ipl_handler.hh"

namespace ost { namespace io {

String IPL::FORMAT_STRING = "defined_ipl";

IPL::IPL(bool normalize_on_save, Format bit_depth):
    ImageFormatBase(FORMAT_STRING),
    normalize_on_save_(normalize_on_save),
    bit_depth_(OST_DEFAULT_FORMAT)
{
  this->SetBitDepth(bit_depth);
}

Format IPL::GetBitDepth() const
{
  return bit_depth_;
}

void IPL::SetBitDepth (Format bitdepth)
{
  if( ! (bitdepth==OST_BIT16_FORMAT || bitdepth==OST_BIT32_FORMAT || bitdepth==OST_DEFAULT_FORMAT))
  {
    throw IOException("Unsupported bit depth for IPL file format.");
  }

  bit_depth_ = bitdepth;
}

bool IPL::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void IPL::SetNormalizeOnSave(bool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}

Real IPL::GetMaximum() const
{
  switch(bit_depth_){
  case OST_BIT32_FORMAT:
    return 4294967295.0;
  default:
    return 65535.0;
  }
}

Real IPL::GetMinimum() const
{
  return 0.0;
}

bool MapIOIPLHandler::MatchContent(unsigned char* header)
{
  String magic_token("DITABIS micron Data File");
  if(magic_token.compare(0,magic_token.size(),reinterpret_cast<char*>(header),magic_token.size())==0)
  {
    return true;
  }
  return false;
}

namespace detail{

class IPLHeader{
public:
  IPLHeader():
    date(),
    header_length(2048),
    size_x(),
    size_y(),
    bit_depth(),
    resolution_x(),
    resolution_y(),
    magnification(1),
    thumb_nail_zoom(10),
    channel("PMT LOWSENS"),
    params("dummy.set"),
    format("2 2 2 2  Standard.fmt"),
    laser(30),
    gain(20000),
    offset_correction(true),
    offset(0),
    comment()
  {}
  IPLHeader(const img::ConstImageHandle& im,Format bit_depth):
    date(),
    header_length(2048),
    size_x(im.GetSize()[0]),
    size_y(im.GetSize()[1]),
    bit_depth(bit_depth==OST_BIT32_FORMAT ? 4: 2),
    resolution_x(im.GetSpatialSampling()[0]),
    resolution_y(im.GetSpatialSampling()[1]),
    magnification(1),
    thumb_nail_zoom(10),
    channel("PMT LOWSENS"),
    params("dummy.set"),
    format("2 2 2 2  Standard.fmt"),
    laser(30),
    gain(20000),
    offset_correction(true),
    offset(0),
    comment()
  {}
  String date;
  int header_length;
  int size_x;
  int size_y;
  int bit_depth;
  Real resolution_x;
  Real resolution_y;
  int magnification;
  int thumb_nail_zoom;
  String channel;
  String params;
  String format;
  int laser;
  int gain;
  bool offset_correction;
  int offset;
  String comment;
};

std::ostream& operator<< (std::ostream& out, const IPLHeader& h )
{
  uint start_pos = out.tellp();
  out << "DITABIS micron Data File\r\n";
  time_t rawtime=time(NULL);
  char * timestr = asctime(localtime(&rawtime));
  timestr[strlen(timestr)-1] = '\0';
  out << "CREATED = "<< timestr <<" \r\n";
  out << "HEADER = "<< h.header_length <<" \r\n";
  // x and y get swapped here (follows the behaviour of the original conversion software)
  out << "YPIXEL = "<< h.size_x <<" \r\n";
  out << "XPIXEL = "<< h.size_y <<" \r\n";
  out << "BYTE PER PIXEL = "<< h.bit_depth <<" \r\n";
  // x and y get swapped here (follows the behaviour of the original conversion software)
  out << "XRESOLUTION = "<< std::setprecision(0)<< h.resolution_y/Units::nm <<" (" << std::setprecision(2)<< h.resolution_y/Units::nm<<") \r\n";
  out << "YRESOLUTION = "<< std::setprecision(0)<< h.resolution_x/Units::nm <<" (" << std::setprecision(2)<< h.resolution_x/Units::nm<<") \r\n";
  out << "MAGNIFICATION = "<< h.magnification <<" \r\n";
  out << "THUMB-NAIL-ZOOM = "<< h.thumb_nail_zoom <<" \r\n";
  out << "CHANNEL = "<< h.channel <<" \r\n";
  out << "PARAMS = "<< h.params <<" \r\n";
  out << "FORMAT = "<< h.format <<" \r\n";
  out << "LASER = "<< h.laser <<" \r\n";
  out << "GAIN = "<< h.gain <<" \r\n";
  if(h.offset_correction){
    out << "OFFSET CORRECTION = YES \r\n";
  }else{
    out << "OFFSET CORRECTION = NO \r\n";
  }
  out << "OFFSET = "<< h.offset <<" \r\n";
  out << "COMMENT = Created by OpenStructure \r\n";
  out << " \r\n";
  uint fillsize=h.header_length-out.tellp()+start_pos;
  std::vector<char> empty(fillsize,'\0');
  out.write(&empty[0],fillsize);
  return out;
}

std::istream& operator>> (std::istream& in, IPLHeader& h)
{
  String line;
  uint start_pos = in.tellg();
  do{
    std::getline(in,line);
    if(line.find("DITABIS micron Data File")!=String::npos){
      //ignore
    }else if(line.find("CREATED")!=String::npos){
      h.date=line.substr(line.find("=")+2);
    }else if(line.find("HEADER")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.header_length;
    // x and y get swapped here (follows the behaviour of the original conversion software)
    }else if(line.find("XPIXEL")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.size_y;
    }else if(line.find("YPIXEL")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.size_x;
    }else if(line.find("BYTE PER PIXEL")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.bit_depth;
      // x and y get swapped here (follows the behaviour of the original conversion software)
    }else if(line.find("XRESOLUTION")!=String::npos){
      std::istringstream ( line.substr(line.find("(")+1) ) >> h.resolution_y;
      h.resolution_y*=Units::nm;
    }else if(line.find("YRESOLUTION")!=String::npos){
      std::istringstream ( line.substr(line.find("(")+1) ) >> h.resolution_x;
      h.resolution_x*=Units::nm;
    }else if(line.find("MAGNIFICATION")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.magnification;
    }else if(line.find("THUMB-NAIL-ZOOM")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.thumb_nail_zoom;
    }else if(line.find("CHANNEL")!=String::npos){
      h.channel=line.substr(line.find("=")+2);
    }else if(line.find("PARAMS")!=String::npos){
      h.params=line.substr(line.find("=")+2);
    }else if(line.find("FORMAT")!=String::npos){
      h.format=line.substr(line.find("=")+2);
    }else if(line.find("LASER")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.laser;
    }else if(line.find("GAIN")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.gain;
    }else if(line.find("OFFSET CORRECTION")!=String::npos){
      if(line.substr(line.find("=")+2).find("YES")!=String::npos){
        h.offset_correction=true;
      }else{
        h.offset_correction=false;
      }
    }else if(line.find("OFFSET")!=String::npos){
      std::istringstream ( line.substr(line.find("=")+2) ) >> h.offset;
    }else if(line.find("COMMENT")!=String::npos){
      h.comment=line.substr(line.find("=")+2);
    }else if(line.find(" ")!=String::npos){
      //ignore
    }else{
      LOG_ERROR("IPL import: unknown header line: " << line);
    }
  }while(in.peek()!=0);
  uint fillsize=h.header_length-in.tellg()+start_pos;
  std::vector<char> empty(fillsize,'\0');
  in.read(&empty[0],fillsize);
  return in;
}

}//ns

bool MapIOIPLHandler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()==IPL::FORMAT_STRING) {
    return true;
  }
  return false;
}

bool MapIOIPLHandler::MatchSuffix(const String& loc)
{
    if(detail::FilenameEndsWith(loc,".IPL") || detail::FilenameEndsWith(loc,".ipl") ) {
      return true;
    }
    return false;
}

void MapIOIPLHandler::Import(img::MapHandle& sh, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct )
{
  boost::filesystem::ifstream infile(loc, std::ios::binary);
  if(!infile) {
    throw IOException("could not open "+loc.string());
  }
  this->Import(sh,infile,formatstruct);
  infile.close();
}

template <typename DATATYPE>
void real_filler(img::image_state::RealSpatialImageState&  isi, std::istream& file)
{
  BinaryIStream<OST_LITTLE_ENDIAN> file_bin(file);
  img::Size size = isi.GetSize();
  char this_dummy; //create dummy variable to give to img::Progress as this
  img::Progress::Instance().Register(&this_dummy,size[1],100);
  for(unsigned int row=0;row<size[1];row++) {
    for(unsigned int column=0;column<size[0];column++) {
      DATATYPE value;
      file_bin >> value;
      isi.Value(img::Point(column,row))=static_cast<Real>(value);
    }
    img::Progress::Instance().AdvanceProgress(&this_dummy);
  }
  img::Progress::Instance().DeRegister(&this_dummy);
}

template <typename DATATYPE>
void real_dumper( const img::ConstImageHandle& sh, std::ostream& file, const IPL& formatIPL, int shrinksize)
{
  img::image_state::RealSpatialImageState *isi=dynamic_cast<img::image_state::RealSpatialImageState*>(sh.ImageStatePtr().get());
  if(! isi){
    throw(IOException("IPL export: dynamic cast failed in real dumper."));
  }
  BinaryOStream<OST_LITTLE_ENDIAN> file_bin(file);
  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if (formatIPL.GetNormalizeOnSave() == true) {
    norm = img::alg::CreateLinearRangeNormalizer(sh,formatIPL.GetMinimum(),formatIPL.GetMaximum());
  }
  img::Size size = isi->GetSize();
  img::ImageHandle thumbnail=sh.Apply(img::alg::DiscreteShrink(img::Size(shrinksize,shrinksize,1)));
  img::image_state::RealSpatialImageState *thumb_isi=dynamic_cast<img::image_state::RealSpatialImageState*>(thumbnail.ImageStatePtr().get());
  if(! thumb_isi){
    throw(IOException("IPL export: dynamic cast failed in real dumper."));
  }

  char this_dummy; //create dummy variable to give to img::Progress as this
  img::Progress::Instance().Register(&this_dummy,size[1]+1,100);
  for(unsigned int row=0;row<size[1];row++) {
    for(unsigned int column=0;column<size[0];column++)
    {
      file_bin << static_cast<DATATYPE>(norm.Convert(isi->Value(ost::img::Point(column,row,0))));
    }
    img::Progress::Instance().AdvanceProgress(&this_dummy);
  }
  img::Progress::Instance().AdvanceProgress(&this_dummy);
  img::Size thumb_size = thumb_isi->GetSize();
  for(unsigned int row=0;row<thumb_size[1];row++) {
    for(unsigned int column=0;column<thumb_size[0];column++)
    {
      file_bin << static_cast<DATATYPE>(norm.Convert(thumb_isi->Value(ost::img::Point(column,row,0))));
    }
    img::Progress::Instance().AdvanceProgress(&this_dummy);
  }
  img::Progress::Instance().DeRegister(&this_dummy);
}


void MapIOIPLHandler::Import(img::MapHandle& sh, std::istream& file, const ImageFormatBase& formatstruct)
{

  IPL form;
  IPL& formatIPL = form;
  if (formatstruct.GetFormatString()==IPL::FORMAT_STRING) {
    formatIPL = formatstruct.As<IPL>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  detail::IPLHeader header;
  file >> header;

  sh.Reset(img::Extent(img::Point(0,0),img::Size(header.size_x,header.size_y)), img::REAL, img::SPATIAL);
  sh.SetSpatialSampling(geom::Vec3(header.resolution_x,header.resolution_y,1.0));
  img::image_state::RealSpatialImageState * isi;
  if(! (isi=dynamic_cast<img::image_state::RealSpatialImageState*>(sh.ImageStatePtr().get()))) {
    throw IOException("internal error in IPL io: expected RealSpatialImageState");
  }

  if(header.bit_depth==4){
    real_filler<uint32_t>(*isi,file);
  }else{
    real_filler<uint16_t>(*isi,file);
  }
}

void MapIOIPLHandler::Export(const img::MapHandle& mh2,
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

void MapIOIPLHandler::Export(const img::MapHandle& sh, std::ostream& file,const ImageFormatBase& formatstruct) const
{

  IPL form;
  IPL& formatIPL = form;
  if (formatstruct.GetFormatString()==IPL::FORMAT_STRING) {
    formatIPL = formatstruct.As<IPL>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }
  if (sh.GetSize()[2]!=1 || sh.GetDomain()!=img::SPATIAL || sh.GetType()!=img::REAL) {
    throw IOException("IPL IO: IPL format only supports spatial 2D images.");
  }
  detail::IPLHeader header(sh,formatIPL.GetBitDepth());
  file << header;
  if(header.bit_depth==4){
    real_dumper<uint32_t>(sh,file,formatIPL,header.thumb_nail_zoom);
  }else{
    real_dumper<uint16_t>(sh,file,formatIPL,header.thumb_nail_zoom);
  }
}

}} // namespaces


