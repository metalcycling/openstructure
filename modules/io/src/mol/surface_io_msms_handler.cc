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
#include <utility>
#include <map>
#include <fstream>
#include <sstream>

#include <ost/log.hh>
#include <boost/filesystem/convenience.hpp>
#include "surface_io_msms_handler.hh"

namespace bf = boost::filesystem;

namespace ost { namespace io {

namespace {

std::pair<bf::path,bf::path> detect_files(const bf::path& loc)
{
  return std::make_pair(bf::change_extension(loc, ".vert"),
                        bf::change_extension(loc, ".face"));
}

}

void SurfaceIOMSMSHandler::Import(mol::SurfaceHandle& sh, const bf::path& loc)
{
  std::pair<bf::path,bf::path> pp = detect_files(loc);
  LOG_DEBUG("loading files " << pp.first.string() << " and " << pp.second.string());

  std::ifstream fvert(pp.first.string().c_str());
  if(!fvert) {
    LOG_ERROR("could not open " << pp.first.string());
    return;
  }
  std::ifstream fface(pp.second.string().c_str());
  if(!fface) {
    LOG_ERROR("could not open " << pp.second.string());
    return;
  }

  String line;
  unsigned int file_vertice_count=0;
  while(std::getline(fvert,line)) {
    if(line[0]!='#') {
      std::istringstream linestream(line);
      int dummy1;
      float dummy2,dummy3;
      linestream >> file_vertice_count >> dummy1 >> dummy2 >> dummy3;
      break;
    }
  }

  unsigned int file_face_count=0;
  while(std::getline(fface,line)) {
    if(line[0]!='#') {
      std::istringstream linestream(line);
      int dummy1;
      float dummy2,dummy3;
      linestream >> file_face_count >> dummy1 >> dummy2 >> dummy3;
      break;
    }
  }

  unsigned int vertice_count=0;
  std::map<unsigned int,mol::SurfaceTriID> vertice_map;
  while(std::getline(fvert,line)) {
    std::istringstream linestream(line);
    geom::Vec3 v,n;
    int x[3];
    linestream >> v[0] >> v[1] >> v[2] >> n[0] >> n[1] >> n[2] >> x[0] >> x[1] >> x[2];
    mol::SurfaceTriID tid = sh.AddVertex(mol::SurfaceVertex(v,n,0));
    ++vertice_count; // msms counts vertices starting at 1
    vertice_map[vertice_count]=tid;
  }

  unsigned int face_count=0;
  while(std::getline(fface,line)) {
    std::istringstream linestream(line);
    int id1,id2,id3;
    int t[2];
    linestream >> id1 >> id2 >> id3 >> t[0] >> t[1];
    ++face_count;
    sh.AddTri(vertice_map[id1],vertice_map[id2],vertice_map[id3]);
  }

  LOG_INFO("loaded " << vertice_count << " vertices and " << face_count << " faces");
}

void SurfaceIOMSMSHandler::Export(const mol::SurfaceHandle& sh,
                                  const bf::path& loc) const
{
  // empty
}

bool SurfaceIOMSMSHandler::ProvidesImport(const bf::path& loc, 
                                          const String& format)
{
  if(format=="auto") {
    std::pair<bf::path,bf::path> pp = detect_files(loc);
    return bf::exists(pp.first) && bf::exists(pp.second);
  } else if(format=="msms") {
    return true;
  }
  return false;
}

bool SurfaceIOMSMSHandler::ProvidesExport(const bf::path& loc, 
                                          const String& format)
{
  return false;
}

}} // ns
