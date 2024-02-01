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
#ifndef OST_IO_SURFACE_IO_PLUGIN_H
#define OST_IO_SURFACE_IO_PLUGIN_H

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>
#include <ost/io/module_config.hh>
#include <ost/mol/surface.hh>

namespace ost { namespace io {

class DLLEXPORT_OST_IO SurfaceIOHandler {
public:
  virtual ~SurfaceIOHandler() {}
  virtual void Import(mol::SurfaceHandle& surf, const boost::filesystem::path& loc) = 0;
  virtual void Export(const mol::SurfaceHandle& ent, const boost::filesystem::path& loc) const = 0;
};

typedef boost::shared_ptr<SurfaceIOHandler> SurfaceIOHandlerPtr;

class DLLEXPORT_OST_IO SurfaceIOHandlerFactoryBase {
public:
  virtual ~SurfaceIOHandlerFactoryBase() {}
  virtual bool ProvidesImport(const boost::filesystem::path& loc, const String& type) const = 0;
  virtual bool ProvidesExport(const boost::filesystem::path& loc, const String& type) const = 0;
  virtual SurfaceIOHandlerPtr Create() const = 0;
  virtual String GetFormatName() const =0;
  virtual String GetFormatDescription() const =0;
};

typedef boost::shared_ptr<SurfaceIOHandlerFactoryBase> SurfaceIOHandlerFactoryBasePtr;

template <class HANDLER>
class SurfaceIOHandlerFactory: public SurfaceIOHandlerFactoryBase
{
  virtual bool ProvidesImport(const boost::filesystem::path& loc, const String& type) const {
    return HANDLER::ProvidesImport(loc,type);
  }

  virtual bool ProvidesExport(const boost::filesystem::path& loc, const String& type) const {
    return HANDLER::ProvidesExport(loc,type);
  }

  virtual SurfaceIOHandlerPtr Create() const {
    return SurfaceIOHandlerPtr(new HANDLER);
  }

  virtual String GetFormatName() const {
     return HANDLER::GetFormatName();
  }

  virtual String GetFormatDescription() const {
     return HANDLER::GetFormatDescription();
  }

};


}} // ns

#endif
