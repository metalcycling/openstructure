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
#ifndef OST_IO_MAP_IO_PLUGIN_H
#define OST_IO_MAP_IO_PLUGIN_H

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <ost/io/module_config.hh>
#include <ost/img/map.hh>
#include <ost/img/alg/normalizer.hh>
#include <ost/io/img/image_format.hh>
#include <ost/io/io_utils.hh>

namespace ost { namespace io {

class DLLEXPORT_OST_IO MapIOHandler {
public:
  virtual ~MapIOHandler() {};
  virtual void Import(img::MapHandle& surf, const boost::filesystem::path& loc, const ImageFormatBase& formatstruct) = 0;
  virtual void Import(img::MapHandle& surf, std::istream& stream,const ImageFormatBase& formatstruct) = 0;
  virtual void Export(const img::MapHandle& ent, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct) const = 0;
  virtual void Export(const img::MapHandle& ent, std::ostream& stream,const ImageFormatBase& formatstruct) const = 0;
};

typedef boost::shared_ptr<MapIOHandler> MapIOHandlerPtr;

class DLLEXPORT_OST_IO MapIOHandlerFactoryBase {
public:
  virtual ~MapIOHandlerFactoryBase() {}
  virtual bool MatchContent(unsigned char* header) const = 0;
  virtual bool MatchType(const ImageFormatBase& type) const = 0;
  virtual bool MatchSuffix(const String& loc) const =0 ;
  virtual MapIOHandlerPtr Create() const = 0 ;
  virtual String GetFormatName() const =0;
  virtual String GetFormatDescription() const =0;
  virtual bool ProvidesImport() const = 0;
  virtual bool ProvidesExport() const = 0;
};

typedef boost::shared_ptr<MapIOHandlerFactoryBase> MapIOHandlerFactoryBasePtr;

template <class HANDLER>
class MapIOHandlerFactory: public MapIOHandlerFactoryBase
{

  virtual bool MatchContent(unsigned char* header) const {
    return HANDLER::MatchContent(header);
  }

  virtual bool MatchType(const ImageFormatBase& type) const {
    return HANDLER::MatchType(type);
  }

  virtual bool MatchSuffix(const String& loc) const {
    return HANDLER::MatchSuffix(loc);
  }

  virtual String GetFormatName() const {
     return HANDLER::GetFormatName();
  }

  virtual String GetFormatDescription() const {
     return HANDLER::GetFormatDescription();
  }
  virtual bool ProvidesImport() const
  {
    return HANDLER::ProvidesImport();
  }
  
  virtual bool ProvidesExport() const
  {
    return HANDLER::ProvidesExport();
  }
  
  virtual MapIOHandlerPtr Create() const {
    return MapIOHandlerPtr(new HANDLER);
  }
};


}} // ns

#endif
