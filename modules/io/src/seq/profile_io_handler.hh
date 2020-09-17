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

/*
  Author: Gerardo Tauriello
 */

#ifndef OST_IO_PROFILE_IO_HANDLER_HH
#define OST_IO_PROFILE_IO_HANDLER_HH

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>

#include <ost/io/module_config.hh>
#include <ost/io/io_utils.hh>
#include <ost/seq/profile_handle.hh>

namespace ost { namespace io {

//! pure abstract base class for profile io handlers
class DLLEXPORT_OST_IO ProfileIOHandler {
public:
  virtual ~ProfileIOHandler() {}

  virtual void Import(seq::ProfileHandle& prof,
                      const boost::filesystem::path& loc) = 0;

  virtual void ImportFromString(seq::ProfileHandle& prof,
                                const String& data) = 0;

  virtual void Export(const seq::ProfileHandle& prof,
                      const boost::filesystem::path& loc) const = 0;
};

typedef boost::shared_ptr<ProfileIOHandler> ProfileIOHandlerPtr;

//! pure abstract base class for creation of a specfic io handler
class DLLEXPORT_OST_IO ProfileIOHandlerFactoryBase {
public:
  virtual ~ProfileIOHandlerFactoryBase() {}
  virtual bool ProvidesImport(const boost::filesystem::path& loc, 
                              const String& format="auto") const = 0;
  virtual bool ProvidesExport(const boost::filesystem::path& loc, 
                              const String& format="auto") const = 0;
  virtual ProfileIOHandlerPtr Create() const = 0;
  virtual String GetFormatName() const =0;
  virtual String GetFormatDescription() const =0;
};

typedef boost::shared_ptr<ProfileIOHandlerFactoryBase> ProfileIOHandlerFactoryBasePtr;

template <class HANDLER>
class ProfileIOHandlerFactory: public ProfileIOHandlerFactoryBase {
  virtual bool ProvidesImport(const boost::filesystem::path& loc, 
                              const String& format="auto") const {
    return HANDLER::ProvidesImport(loc, format);
  }

  virtual bool ProvidesExport(const boost::filesystem::path& loc, 
                              const String& format="auto") const {
    return HANDLER::ProvidesExport(loc, format);
  }

  virtual ProfileIOHandlerPtr Create() const {
    return ProfileIOHandlerPtr(new HANDLER);
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
