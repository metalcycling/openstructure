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
  Author: Ansgar Philippsen
*/

#include <fstream>
#include <sstream>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/filesystem/fstream.hpp>

#include <ost/base.hh>
#include <ost/img/normalizer_impl.hh>
#include <ost/img/image.hh>
#include <ost/img/progress.hh>
#include <ost/units.hh>
#include <ost/log.hh>
#include <ost/io/io_exception.hh>



#include "map_io_nanoscope_handler.hh"

namespace ost { namespace io {

String Nanoscope::FORMAT_STRING = "defined_nanoscope";

Nanoscope::Nanoscope(int subimage):
  ImageFormatBase(FORMAT_STRING),
  subimage_(subimage)
{
}

int Nanoscope::GetSubimage() const
{
  return subimage_;
}

void Nanoscope::SetSubimage (int subimage)
{
  subimage_ = subimage;
}

namespace detail {

struct NHeader
{
  NHeader():
    multih(1.0), multih_f(false),
    multid(1.0), multid_f(false),
    multi(1.0), multi_f(false),
    nmsize(1.0), nmsize_f(false),
    px(0), px_f(false),
    py(0), py_f(false),
    sampling(1.0), sampling_f(false),
    doffset(0), doffset_f(false),
    dlength(0), dlength_f(false)
  {}

  Real multih;
  bool multih_f;
  Real multid;
  bool multid_f;
  Real multi;
  bool multi_f;
  Real nmsize;
  bool nmsize_f;
  uint px;
  bool px_f;
  uint py;
  bool py_f;
  Real sampling;
  bool sampling_f;
  int doffset;
  bool doffset_f;
  int dlength;
  bool dlength_f;
  Real zscale;
  Real zscale_f;
};

void header_filler(NHeader& h, std::istream& in, int inum)
{
  static const boost::regex r_image("Operating mode: ([a-zA-Z]*)");
  static const boost::regex r_itag("\\*Ciao image list");
  static const boost::regex r_multih("@Sens. Zscan: V (.*) nm/V");
  static const boost::regex r_nmsize("Scan size: (.*) (.*) nm");
  static const boost::regex r_multid("@Sens. Deflection: V (.*) nm/V");
  static const boost::regex r_doffset("Data offset: ([[:d:]]*)");
  static const boost::regex r_dlength("Data length: ([[:d:]]*)");
  static const boost::regex r_px("Samps/line: ([[:d:]]*)");
  static const boost::regex r_py("Number of lines: ([[:d:]]*)");
  static const boost::regex r_itype("@2:Image Data: S \\[(.+)\\] \"(.+)\"");
  static const boost::regex r_zscale("@2:Z scale: V \\[Sens. Zscan\\] \\((.*) V/LSB\\) (.*) V");

  String line;
  boost::match_results<std::string::const_iterator> match;
  boost::match_flag_type match_flags = boost::match_default;
  int icount=0;
  while(std::getline(in,line)) {
    if(line[0]=='\\') {
      line[0]=' ';
      //std::cerr << line << std::endl;
      if(regex_search(line, match, r_image, match_flags)) {
  if(match[1]!="Image") {
    throw IOException("nanoscope_io: operating mode must be 'Image'");
  }
      } else if(regex_search(line, match, r_itag, match_flags)) {
  ++icount;
      } else if(regex_search(line, match, r_multih, match_flags)) {
  h.multih = boost::lexical_cast<Real>(match[1]);
  h.multih_f=true;
      } else if(regex_search(line, match, r_nmsize, match_flags)) {
  h.nmsize = boost::lexical_cast<Real>(match[1]);
  h.nmsize_f=true;
      } else if(regex_search(line, match, r_multid, match_flags)) {
  h.multid = boost::lexical_cast<Real>(match[1]);
  h.multid_f=true;
      } else if(regex_search(line, match, r_doffset, match_flags)) {
  if(icount==inum) {
    h.doffset = boost::lexical_cast<int>(match[1]);
    h.doffset_f=true;
  }
      } else if(regex_search(line, match, r_dlength, match_flags)) {
  if(icount==inum) {
    h.dlength = boost::lexical_cast<int>(match[1]);
    h.dlength_f=true;
  }
      } else if(regex_search(line, match, r_px, match_flags)) {
  if(icount==inum) {
    h.px = boost::lexical_cast<int>(match[1]);
    h.px_f=true;
  }
      } else if(regex_search(line, match, r_py, match_flags)) {
  if(icount==inum) {
    h.py = boost::lexical_cast<int>(match[1]);
    h.py_f=true;
  }
      } else if(regex_search(line, match, r_itype, match_flags)) {
  if(icount==inum) {
    if(match[1]=="Height") {
      h.multi=h.multih;
    } else {
      h.multi=h.multid;
    }
    h.multi_f=true;
  }
      } else if(regex_search(line, match, r_zscale, match_flags)) {
  if(icount==inum) {
    h.zscale = boost::lexical_cast<Real>(match[1]);
    h.zscale_f=true;
  }
      }
    }
  }
}

void print_header (const NHeader& header, int inum)
{
  if(Logger::Instance().GetVerbosityLevel()>2) {
    LOG_INFO("io_nanoscope: header dump for image " << inum);
    LOG_INFO(" px      : " << header.px);
    LOG_INFO(" py      : " << header.py);
    LOG_INFO(" nmsize  : " << header.nmsize);
    LOG_INFO(" sampling: " << header.sampling << " nm");
    LOG_INFO(" zscale  : " << header.zscale);
    LOG_INFO(" multih  : " << header.multih);
    LOG_INFO(" multid  : " << header.multid);
    LOG_INFO(" multi   : " << header.multi);
    LOG_INFO(" doffset : " << header.doffset);
    LOG_INFO(" dlength : " << header.dlength);
  }
}


void real_filler(img::ImageHandle& image, std::istream& ff,
                 const NHeader& header)
{
  ff.clear();
  ff.seekg(0);
  ff.seekg(header.doffset, std::ios::beg);

  std::vector<char> bdata(header.px*2); // short format
  for(uint row=0;row<header.py;++row) {

    ff.read(&bdata[0],bdata.size());
    if(ff.fail()) {
      std::stringstream msg;
      msg << boost::format("unexpected failure during stream input at row %1% / %2%") %row %header.py ;
      throw IOException(msg.str());
    }
    short* pdata = reinterpret_cast<short*>(&bdata[0]);
    for(uint col=0;col<header.px;++col) {
      image.SetReal(img::Point(col,header.py-row-1),static_cast<Real>(pdata[col]));
    }
  }

  image*=header.multi;
}

} // ns detail

bool MapIONanoscopeHandler::MatchContent(unsigned char* header)
{
  return false;
}

bool MapIONanoscopeHandler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()=="defined_mrc") {
    return true;
  }
  return false;
}

bool MapIONanoscopeHandler::MatchSuffix(const String& loc)
{
    if(detail::FilenameEndsWith(loc,".mod") ) {
      return true;
    }
    return false;
}

void MapIONanoscopeHandler::Import(img::MapHandle& sh, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct )
{
  boost::filesystem::ifstream infile(loc, std::ios::binary);
  if(!infile)
  {
    throw IOException("could not open "+loc.string());
  }
  this->Import(sh,infile,formatstruct);
  infile.close();
}

void MapIONanoscopeHandler::Import(img::MapHandle& image, std::istream& ff, const ImageFormatBase& formatstruct)
{
  Nanoscope form;
  Nanoscope& formatnanoscope = form;
  if (formatstruct.GetFormatString()==Nanoscope::FORMAT_STRING) {
    formatnanoscope = formatstruct.As<Nanoscope>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  int inum = formatnanoscope.GetSubimage();
  inum = (inum<1) ? 1 : inum;

  detail::NHeader header;

  header_filler(header,ff,inum);

  std::ostringstream msg;

  if(!header.multi_f) {
    msg << "nanoscope_io: image type not found" << std::endl;
  }
  if(!header.zscale_f) {
    msg << "nanoscope_io: Z scaling factor not found" << std::endl;
  } else {
    header.multi*=header.zscale;
  }
  if(!header.multih_f) {
    msg << "nanoscope_io: piezo z-calibration not found" << std::endl;
  }
  if(!header.multid_f) {
    msg << "nanoscope_io: deflection sensitivity not found" << std::endl;
  }
  if(!header.nmsize_f) {
    msg << "nanoscope_io: image size not found" << std::endl;
  }
  if(!header.px_f) {
    msg << "nanoscope_io: points per line not found" << std::endl;
  }
  if(!header.py_f) {
    msg << "nanoscope_io: number of lines not found" << std::endl;
  }
  if(!header.doffset_f) {
    msg << "nanoscope_io: data offset not found" << std::endl;
  }
  if(!header.dlength_f) {
    msg << "nanoscope_io: data length not found" << std::endl;
  }

  if(!msg.str().empty()) {
    throw IOException (msg.str());
  }

  header.sampling = header.nmsize/static_cast<Real>(header.px);

  detail::print_header(header,inum);

  image.Reset(img::Size(header.px,header.py), img::REAL, img::SPATIAL);
  image.SetSpatialSampling(header.sampling*Units::nm);


  detail::real_filler(image, ff,header);

}

void  MapIONanoscopeHandler::Export(const img::MapHandle& sh, const boost::filesystem::path& loc, const ImageFormatBase& formatstruct) const
{
  throw IOException("This format is not supported for Export");
}


void MapIONanoscopeHandler::Export(const img::MapHandle& sh, std::ostream& loc,const ImageFormatBase& formatstruct) const
{
  throw IOException("This format is not supported for Export");
}



}} // namespaces
