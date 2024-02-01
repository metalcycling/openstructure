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
  Author: Marco Biasini
*/

#include <cstdio>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>
#include <fstream>
#include <ost/log.hh>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <ost/io/swap_util.hh>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/converting_streams.hh>
#include <ost/img/progress.hh>

#include "map_io_df3_handler.hh"

namespace bf = boost::filesystem;

namespace ost { namespace io {

String DF3::FORMAT_STRING="defined_df3";

DF3::DF3(bool normalize_on_save):
    ImageFormatBase(FORMAT_STRING),
    normalize_on_save_(normalize_on_save)
{
}

bool DF3::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void DF3::SetNormalizeOnSave(bool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}

Real DF3::GetMaximum() const
{
  return std::numeric_limits<unsigned short>::max();
}

Real DF3::GetMinimum() const
{
  return 0;
}

void DF3MapIOHandler::Import(img::MapHandle& mh, const bf::path& loc,
                             const ImageFormatBase& formatstruct)
{
}

void DF3MapIOHandler::Import(img::MapHandle& mh, std::istream& infile,
                             const ImageFormatBase& formatstruct)
{
}


void DF3MapIOHandler::Export(const img::MapHandle& mh2,
                             const bf::path& loc,
                             const ImageFormatBase& format) const
{
  boost::filesystem::ofstream outfile(loc, std::ios::binary);
  if(!outfile)
  {
    throw IOException("could not open "+loc.string());
  }
  BinaryOStream<OST_BIG_ENDIAN> big_endian_stream(outfile);
  this->Export(mh2, big_endian_stream, format);
  outfile.close();

}  
  
void DF3MapIOHandler::Export(const img::MapHandle& mh2,
                             std::ostream& outfile,
                             const ImageFormatBase& format) const
{  
  DF3 default_df3;
  DF3& fmt=default_df3;
  if (format.GetFormatString()==DF3::FORMAT_STRING) {
    fmt=format.As<DF3>();
  } else {
    assert(format.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }
  if (mh2.GetType()==img::COMPLEX) {
    throw IOException("DF3 doesn't support complex-valued maps");
  }
  img::alg::Normalizer norm=img::alg::CreateNoOpNormalizer();
  if (fmt.GetNormalizeOnSave()) {
    
    norm=img::alg::CreateLinearRangeNormalizer(mh2, fmt.GetMinimum(),
                                               fmt.GetMaximum());
  }
  img::Size size=mh2.GetSize();
  for (size_t i=0; i<3; ++i) {
    unsigned short v=size[i];
    Convert<OST_BIG_ENDIAN,unsigned short>::FromIP(&v);
    outfile.write(reinterpret_cast<const char*>(&v), sizeof(unsigned short));
  }
  for (img::ExtentIterator i(mh2.GetExtent()); !i.AtEnd(); ++i) {
    Real norm_value=norm.Convert(mh2.GetReal(i));
    unsigned short v=static_cast<unsigned short>(norm_value);
    Convert<OST_BIG_ENDIAN,unsigned short>::FromIP(&v);
    outfile.write(reinterpret_cast<const char*>(&v), sizeof(unsigned short));
  }
}

bool DF3MapIOHandler::MatchContent(unsigned char* header)
{
  return false;
}

bool DF3MapIOHandler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()==DF3::FORMAT_STRING) {
    return true;
  }
  return false;
}

bool DF3MapIOHandler::MatchSuffix(const String& loc)
{
  return detail::FilenameEndsWith(loc,".df3");
}



}} //namespaces
