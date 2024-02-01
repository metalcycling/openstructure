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
#include <cstdio>
#include <utility>
#include <map>
#include <fstream>
#include <sstream>

#include <ost/log.hh>
#include <ost/string_ref.hh>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/format.hpp>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/io/io_exception.hh>
#include <ost/img/map.hh>

#include "map_io_dx_handler.hh"

namespace bf = boost::filesystem;

namespace ost { namespace io {

using boost::format;

namespace {

/*

commented out to silence compiler warning, as it is not needed anywhere

bool IEquals(const StringRef& a, const StringRef& b)
{
  if (a.size()!=b.size()) {
    return false;
  }
  for (size_t i=0; i<a.size(); ++i) {
    if (toupper(a[i])!=b[i]) {
      return false;
    }
  }
  return true;
}
*/

}

String DX::FORMAT_STRING="defined_dx";

DX::DX (bool normalize_on_save):
  ImageFormatBase(FORMAT_STRING),
  normalize_on_save_(normalize_on_save) 
{
}

bool DX::GetNormalizeOnSave() const
{
  return normalize_on_save_;
}

void DX::SetNormalizeOnSave(bool normalize_on_save)
{
  normalize_on_save_ = normalize_on_save;
}

void MapIODxHandler::Import(img::MapHandle& mh, const bf::path& loc,const ImageFormatBase& form)
{
  boost::filesystem::ifstream infile(loc);
  if(!infile)
  {
    throw IOException("could not open "+loc.string());
  }
  boost::iostreams::filtering_stream<boost::iostreams::input> in;
  if (boost::iequals(".gz", boost::filesystem::extension(loc))) {
    in.push(boost::iostreams::gzip_decompressor());
  }
  in.push(infile);
  this->Import(mh,in,form);
  infile.close();
}

void MapIODxHandler::Import(img::MapHandle& mh, std::istream& infile, const ImageFormatBase& formatstruct)
{
  DX form;
  DX& formatdx = form;
  if (formatstruct.GetFormatString()==DX::FORMAT_STRING) {
    formatdx = formatstruct.As<DX>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }

  String line;
  int u_size=0,v_size=0,w_size=0;
  Real x_orig=0.0, y_orig=0.0, z_orig=0.0;
  Real u_spacing=1.0, v_spacing=1.0, w_spacing=1.0;
  size_t delta_count=0;
  Real deltas[][3]={{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
  int num_gridpoints;
  img::MapHandle mh2;
  std::vector<String> tokens;
  while (std::getline(infile,line)) {
    if (line.empty()) {
      continue;
    }
    // read gridpoints line
    if (boost::iequals(line.substr(0,35), "object 1 class gridpositions counts")) {
      LOG_DEBUG("DXImport: reading gridpoints line [" << line << "]");
      boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
      int tokens_size = (int)tokens.size();
      if(tokens_size < 3) {
        String msg="Bad gridpoints line: Can't read number of grid points";
        throw IOException(msg);
      }
      try {
        u_size=boost::lexical_cast<int>(boost::trim_copy(tokens[tokens_size-3]));
        v_size=boost::lexical_cast<int>(boost::trim_copy(tokens[tokens_size-2]));
        w_size=boost::lexical_cast<int>(boost::trim_copy(tokens[tokens_size-1]));
      } catch(boost::bad_lexical_cast&) {
        format fmer  = format("Bad gridpoints line: Can't convert number of grid points '%s' to integral constant.") % line;
        throw IOException(fmer.str());
      }
    }
    // read grid origin line
    else if (boost::iequals(line.substr(0,6), "origin")) {
      LOG_DEBUG("DXImport: reading origin line [" << line << "]");
      boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
      int tokens_size = (int)tokens.size();
      if(tokens_size < 3) {
        String msg="Bad center line: Can't read number of grid points";
        throw IOException(msg);
      }
      try {
        x_orig=boost::lexical_cast<Real>(boost::trim_copy(tokens[tokens_size-3]));
        y_orig=boost::lexical_cast<Real>(boost::trim_copy(tokens[tokens_size-2]));
        z_orig=boost::lexical_cast<Real>(boost::trim_copy(tokens[tokens_size-1]));
      } catch(boost::bad_lexical_cast&) {
        format fmer  = format("Bad origin line: Can't convert String of origin '%s' to Real constant.") % line;
        throw IOException(fmer.str());
      }
    }
    // read grid spacing line
    else if (boost::iequals(line.substr(0,5), "delta")) {
      LOG_DEBUG("DXImport: reading delta line [" << line << "]");
      boost::split(tokens, line, boost::is_any_of("\t "), boost::token_compress_on);
      if(tokens.size() < 4) {
        String msg="DXImport: expected delta line with 3 floats";
        throw IOException(msg);
      }
      try {
        if(delta_count<3) {
          for(size_t i=0;i<3;++i) {
            deltas[delta_count][i] = boost::lexical_cast<Real>(tokens[i+1]);
          }
          ++delta_count;
        }
      } catch(boost::bad_lexical_cast&) {
        format fmer= format("Bad spacing line: Can't convert String of delta"
                   " '%s' to Real constants.") % line;
        throw IOException(fmer.str());
      }
    }
    // read number of data points
    if (boost::iequals(line.substr(0,25), "object 3 class array type")) {
      boost::split(tokens, line, boost::is_any_of(" "), boost::token_compress_on);
      int tokens_size = (int)tokens.size();
      try {
        num_gridpoints=boost::lexical_cast<int>(boost::trim_copy(tokens[tokens_size-3]));
      } catch(boost::bad_lexical_cast&) {
        format fmer = format ("Bad gridpoints line: Can't convert number of grid points '%s' to integral constant.") % line;
        throw IOException(fmer.str());
      }
      if(u_size*v_size*w_size != num_gridpoints) {
        format fmer = format("Bad gridpoints line: Number of gridpoints (%i) does not"
                   " correlate with the grid size (%i, %i, %i)") % num_gridpoints % u_size % v_size % w_size;
        throw IOException(fmer.str());
      }

      // at this point enough info should be available to create map
      if(u_size>1e5 || v_size>1e5 || w_size>1e5) {
        format fmer=format("DXImport: nonsense mapsize read (%d %d %d)") % u_size % v_size % w_size;
        throw IOException(fmer.str());
      }
      LOG_DEBUG("DXImport: creating map of size " << img::Size(u_size,v_size,w_size));
      LOG_DEBUG("          at absolute origin " << geom::Vec3(x_orig, y_orig, z_orig));
      mh2 = CreateMap(img::Size(u_size,v_size,w_size));
      mh2.SetAbsoluteOrigin(geom::Vec3(x_orig, y_orig, z_orig));
      if(delta_count==3) {
        u_spacing=deltas[0][0];
        v_spacing=deltas[1][1];
        w_spacing=deltas[2][2];
        LOG_DEBUG("          and spatial sampling " << geom::Vec3(u_spacing,v_spacing,w_spacing));
        mh2.SetSpatialSampling(geom::Vec3(u_spacing,v_spacing,w_spacing));
      } 

      // and read in the values
      // TODO: this is glacially slow
      Real value=0;
      for(int i=0; i<num_gridpoints; i+=3) {
        std::getline(infile,line);
        StringRef curr_line(line.c_str(), line.size());
        std::vector<StringRef> fields=curr_line.split(' ');
        for (size_t j=0; j<fields.size(); j++) {
          std::pair<bool, float> result=fields[j].trim().to_float();
          if (!result.first) {
            throw IOException((format("Bad value line: Can't convert grid point value '%s' to Real constant.") % line).str());
          }
          value=result.second;
          mh2.SetReal(img::Point(((i+j)/(v_size*w_size))%u_size,((i+j)/w_size)%v_size, (i+j)%w_size), value);
        }
      }
    }
  }
  mh.Swap(mh2);
}


void MapIODxHandler::Export(const img::MapHandle& mh2,
                                  const bf::path& loc,const ImageFormatBase& formatstruct) const
{
  boost::filesystem::ofstream outfile(loc);
  if(!outfile) {
    throw IOException("could not open "+loc.string()+" for writing");
  }

  this->Export(mh2,outfile,formatstruct);
}

void MapIODxHandler::Export(const img::MapHandle& mh2,
                                  std::ostream& outfile,const ImageFormatBase& formatstruct) const
{
  DX form;
  DX& formatdx = form;
  if (formatstruct.GetFormatString()==DX::FORMAT_STRING) {
    formatdx = formatstruct.As<DX>();
  } else {
    assert (formatstruct.GetFormatString()==UndefinedImageFormat::FORMAT_STRING);
  }
  std::ostringstream outstring;
  int u_size = mh2.GetExtent().GetSize()[0];
  int v_size = mh2.GetExtent().GetSize()[1];
  int w_size = mh2.GetExtent().GetSize()[2];
  int num_gridpoints = u_size*v_size*w_size;
  outstring << "# Data created by OpenStructure" << std::endl;
  outstring << "object 1 class gridpositions counts " << u_size << " " << v_size << " " << w_size << std::endl;
  geom::Vec3 ori=mh2.IndexToCoord(img::Point(0,0,0));
  geom::Vec3 spa=mh2.GetSpatialSampling();
  outstring << format("origin %.6e %.6e %.6e") % ori[0] % ori[1] % ori[2] << std::endl;
  outstring << format("delta %.6e 0.000000e+00 0.000000e+00") % spa[0] << std::endl;
  outstring << format("delta 0.000000e+00 %.6e 0.000000e+00") % spa[1] << std::endl; 
  outstring << format("delta 0.000000e+00 0.000000e+00 %.6e") % spa[2] << std::endl; 
  outstring << "object 2 class gridconnections counts " << u_size << " " << v_size << " " << w_size << std::endl;
  outstring << "object 3 class array type Real rank 0 items " << num_gridpoints << " data follows" << std::endl;

  img::alg::Normalizer norm = img::alg::CreateNoOpNormalizer();
  if (formatdx.GetNormalizeOnSave() == true) {
    norm = img::alg::CreateLinearRangeNormalizer(mh2,formatdx.GetMinimum(),formatdx.GetMaximum());
  }

  for(int i=0; i<num_gridpoints; i+=3) {
    for(int j=0; j<3; j++) {  // three values per line
      if((i+j)<num_gridpoints) {
        outstring << format("%.6e ") % norm.Convert(mh2.GetReal(img::Point(((i+j)/(v_size*w_size))%u_size,((i+j)/w_size)%v_size, (i+j)%w_size)));
      }
    }
    outstring << std::endl;
  }

  outfile << outstring.str();
}

bool MapIODxHandler::MatchContent(unsigned char* header)
{
  return false;
}
bool MapIODxHandler::MatchType(const ImageFormatBase& formatstruct)
{
  if(formatstruct.GetFormatString()=="defined_dx") {
    return true;
  }
  return false;
}
bool MapIODxHandler::MatchSuffix(const String& loc)
{
	if(detail::FilenameEndsWith(loc,".dx") || detail::FilenameEndsWith(loc,".dx.gz")) {
    return true;
  }
  return false;
}

}} // ns
