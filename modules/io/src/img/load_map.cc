//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
//
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
#include <ost/log.hh>
#include <ost/img/map.hh>
#include <ost/img/alg/stat_min_max.hh>
#include <ost/io/io_manager.hh>
#include "load_map.hh"
#include "map_io_handler.hh"

namespace ost { namespace io {

img::ImageHandle LoadImage(const boost::filesystem::path& loc)
{
  UndefinedImageFormat undefined;
  return LoadImage(loc,undefined);
}

img::ImageHandle LoadImage(const boost::filesystem::path& loc, const ImageFormatBase& formatstruct)
{
  if(!boost::filesystem::exists(loc)){
    throw IOException("file not found: " + loc.string());
  }

  LOG_DEBUG("creating MapIOHandler for " << loc.string());

  MapIOHandlerPtr map_io = IOManager::Instance().FindMapImportHandlerFile(loc,formatstruct);

  if(!map_io) {
    throw IOUnknownFormatException("could not find io-plugin for " + loc.string());
  }

  img::ImageHandle ih = CreateImage(img::Extent(),img::REAL,img::SPATIAL);

  LOG_DEBUG("calling import on map io handle");
  map_io->Import(ih,loc,formatstruct);

  // TODO: fix this in regard to automatic copy of ConstData info obj
  //info::InfoGroup g = ih.Info().Root().RetrieveGroup("file");
  //g.CreateItem("filename",filename);

  return ih;
}

void SaveImage(const img::ImageHandle& image, const boost::filesystem::path& loc)
{
  UndefinedImageFormat undefined;
  SaveImage(image,loc,undefined);
}

void SaveImage(const img::ImageHandle& image, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct)
{
  MapIOHandlerPtr map_io = IOManager::Instance().FindMapExportHandlerFile(loc,formatstruct);

  if(!map_io) {
    throw IOUnknownFormatException("could not find io-plugin for " + loc.string());
  }
  LOG_DEBUG("calling export on map io handle");
  map_io->Export(image, loc,formatstruct);
}


img::MapHandle LoadMap(const boost::filesystem::path& loc, const ImageFormatBase& formatstruct)
{
  return LoadImage(loc,formatstruct);
}

img::MapHandle LoadMap(const boost::filesystem::path& loc)
{
  return LoadImage(loc);
}

void SaveMap(const img::ImageHandle& image, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct )
{
  SaveImage(image,loc,formatstruct);
}

void SaveMap(const img::ImageHandle& image, const boost::filesystem::path& loc)
{
  SaveImage(image,loc);
}
}} // ns
