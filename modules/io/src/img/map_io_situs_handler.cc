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
  Author: Valerio Mariani, Ansgar Philippsen
*/

#include <cstdio>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <ost/log.hh>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <ost/io/swap_util.hh>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/io/io_exception.hh>
#include <ost/img/progress.hh>

#include "map_io_situs_handler.hh"

// TODO use pixel sampling from header

namespace bf = boost::filesystem;

namespace ost { namespace io {

String Situs::FORMAT_STRING="defined_situs";

Situs::Situs(bool normalize_on_save):
    ImageFormatBase(FORMAT_STRING),
    normalize_on_save_(normalize_on_save)
{
}

bool Situs::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void Situs::SetNormalizeOnSave(bool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}


namespace detail {

using boost::format;

struct situs_header {
  Real dwidth;
  Real gridx,gridy,gridz;
  unsigned extx,exty,extz;
};

void prep_header(situs_header& header, const img::Size& size, const geom::Vec3& pixel_sampling, const geom::Vec3& start_coord)
{
  header.dwidth=pixel_sampling[0];
  header.extx=size[0];
  header.exty=size[1];
  header.extz=size[2];
  header.gridx=start_coord[0];
  header.gridy=start_coord[1];
  header.gridz=start_coord[2];
}

void print_header(const situs_header& h)
{
  std::ostringstream s;
  s << "situs header: " << std::endl;
  s << format(" voxel width: %1%") % h.dwidth << std::endl;
  s << format(" size x y z: %1% %2% %3%") % h.extx % h.exty % h.extz << std::endl;
  s << format(" origin x y z: %1% %2% %3%") % h.gridx % h.gridy % h.gridz;
  LOG_INFO(s.str());
}


void header_filler(std::istream& f, situs_header& header)
{
  f >> header.dwidth;
  f >> header.gridx;
  f >> header.gridy;
  f >> header.gridz;
  f >> header.extx;
  f >> header.exty;
  f >> header.extz;
}

void header_dumper(std::ostream& f, const situs_header& header)
{
  f << str(format ("%10f") % header.dwidth) << " ";
  f << str(format ("%10f") % header.dwidth) << header.gridx << " ";
  f << str(format ("%10f") % header.dwidth) << header.gridy << " ";
  f << str(format ("%10f") % header.dwidth) << header.gridz << " ";
  f << header.extx << " ";
  f << header.exty << " ";
  f << header.extz << " ";
  f << "\n";
  f << "\n";
}

void real_filler(std::istream& f, const situs_header& header, img::ImageHandle& mh,Real scale)
{
  Real raw;
  img::Point pnt;
  char this_dummy; //create dummy variable to give to Progress as this
  img::Progress::Instance().Register(&this_dummy,header.extx*header.exty*header.extz,100);
  img::Size size = mh.GetSize();
  for(uint z_i=0;z_i<size[2];++z_i) {
    for(uint y_i=0;y_i<size[1];++y_i) {
      for(uint x_i=0;x_i<size[0];++x_i) {
        f >> raw;
        mh.SetReal(img::Point(x_i,y_i,z_i),raw);
      }
      img::Progress::Instance().AdvanceProgress(&this_dummy);
    }
  }
  img::Progress::Instance().DeRegister(&this_dummy);
}

void real_dumper(std::ostream& f,  const situs_header& header, const img::ImageHandle& mh, const  img::alg::Normalizer& norm)
{
  char this_dummy; //create dummy variable to give to Progress as this
  img::Progress::Instance().Register(&this_dummy,header.extx*header.exty*header.extz,100);
  int counter=0;
  uint ext_z=header.extz;
  uint ext_y=header.exty;
  uint ext_x=header.extx;
  for(uint z_i=0;z_i<ext_z;++z_i) {
    for(uint y_i=0;y_i<ext_y;++y_i) {
      for(uint x_i=0;x_i<ext_x;++x_i) {
        Real raw=norm.Convert(mh.GetReal(img::Point(x_i,y_i,z_i)));
        f << " " << format("%10f") % raw << " ";
        ++counter;
        if (counter == 10) {
          f << "\n";
          counter=0;
        }
      }
      img::Progress::Instance().AdvanceProgress(&this_dummy);
    }
  }
  f << "\n";
  img::Progress::Instance().DeRegister(&this_dummy);
}

} // detail



void MapIOSitusHandler::Import(img::MapHandle& mh, const bf::path& loc,const ImageFormatBase& formatstruct)
{
  boost::filesystem::ifstream infile(loc);
  if(!infile)
  {
    throw IOException("could not open "+loc.string());
  }

  this->Import(mh,infile,formatstruct);
  infile.close();
}

void MapIOSitusHandler::Import(img::MapHandle& mh, std::istream& infile,const ImageFormatBase& formatstruct)
{
  Situs form;
  Situs& formatsitus = form;
  if (formatstruct.GetFormatString()==Situs::FORMAT_STRING) {
    formatsitus = formatstruct.As<Situs>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  detail::situs_header header;
  detail::header_filler(infile,header);
  detail::print_header(header);

  img::Size msize;
  msize[0]=header.extx;
  msize[1]=header.exty;
  msize[2]=header.extz;

  geom::Vec3 mstart;
  mstart[0]=header.gridx;
  mstart[1]=header.gridy;
  mstart[2]=header.gridz;

  // MAYBE REFACTOR

  mh.Reset(msize);
  mh.SetSpatialSampling(geom::Vec3(header.dwidth,header.dwidth,header.dwidth));

  LOG_VERBOSE("resulting image extent: " << mh.GetExtent());

  detail::real_filler(infile,header,mh,1.0);

}


void MapIOSitusHandler::Export(const img::MapHandle& mh2,
                                  const bf::path& loc,const ImageFormatBase& formatstruct) const
{
  boost::filesystem::ofstream outfile(loc);
  if(!outfile)
  {
    throw IOException("could not open "+loc.string());
  }

  this->Export(mh2,outfile,formatstruct);
  outfile.close();
  
}  
  
void MapIOSitusHandler::Export(const img::MapHandle& mh2,
                                  std::ostream& outfile,const ImageFormatBase& formatstruct) const
{  
  Situs form;
  Situs& formatsitus = form;
  if (formatstruct.GetFormatString()==Situs::FORMAT_STRING) {
    formatsitus = formatstruct.As<Situs>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if (formatsitus.GetNormalizeOnSave() == true) {
    norm = img::alg::CreateLinearRangeNormalizer(mh2,formatsitus.GetMinimum(),formatsitus.GetMaximum());
  }

  detail::situs_header header;
  geom::Vec3 spacing=mh2.GetSpatialSampling();
  if( (spacing[0] != spacing[1]) ||
      (spacing[1] != spacing[2]) ||
      (spacing[2] != spacing[0]))
  {
    throw IOException("this format does not support anisotropic spatial sampling");
  }

  prep_header(header,mh2.GetSize(),spacing,mh2.GetAbsoluteOrigin());
  detail::header_dumper(outfile,header);
  detail::real_dumper(outfile,header,mh2,norm);
}

bool MapIOSitusHandler::MatchContent(unsigned char* header)
{
  return false;
}

bool MapIOSitusHandler::MatchType(const ImageFormatBase& type)
{
  if(type.GetFormatString()=="defined_situs") {
    return true;
  }
  return false;
}

bool MapIOSitusHandler::MatchSuffix(const String& loc)
{
  if(detail::FilenameEndsWith(loc,".situs") || detail::FilenameEndsWith(loc,".sit") ) {
    return true;
  }
  return false;
}



}} //namespaces
