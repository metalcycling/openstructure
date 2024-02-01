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
#ifndef OST_IO_ENTITY_IO_PLUGIN_H
#define OST_IO_ENTITY_IO_PLUGIN_H

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>


#include <ost/io/module_config.hh>
#include <ost/io/io_utils.hh>

namespace ost { 
  
namespace mol {

class EntityView;
class EntityHandle;

}

namespace io {

//! pure abstract base class for entity io handlers
class DLLEXPORT_OST_IO EntityIOHandler {
public:
  virtual ~EntityIOHandler() {}

  // import data from file into provided entity handle
  /*
    the handler is expected to convert atom names to IUPAC standard
  */
  virtual void Import(mol::EntityHandle& ent, 
                      const boost::filesystem::path& loc)=0;
  // export data from entity view to provided file
  virtual void Export(const mol::EntityView& ent, 
                      const boost::filesystem::path& loc) const = 0;
  // import data from provided stream
  virtual void Import(mol::EntityHandle& ent, 
                      std::istream& stream)=0;
  // export data from entity view to provided stream
  virtual void Export(const mol::EntityView& ent, 
                      std::ostream& stream) const=0;
  virtual bool RequiresProcessor() const=0;
};

typedef boost::shared_ptr<EntityIOHandler> EntityIOHandlerP;


//! pure abstract base class for creation of a specfic entity io handle
class DLLEXPORT_OST_IO EntityIOHandlerFactoryBase {
public:
  virtual ~EntityIOHandlerFactoryBase() {}
  virtual bool ProvidesImport(const boost::filesystem::path& loc, const String& type) const = 0;
  virtual bool ProvidesExport(const boost::filesystem::path& loc, const String& type) const = 0;
  virtual EntityIOHandlerP Create() const = 0;
  virtual String GetFormatName() const =0;
  virtual String GetFormatDescription() const =0;
};

typedef boost::shared_ptr<EntityIOHandlerFactoryBase> EntityIOHandlerFactoryBaseP;

template <class HANDLER>
class EntityIOHandlerFactory: public EntityIOHandlerFactoryBase
{
  virtual bool ProvidesImport(const boost::filesystem::path& loc, const String& type) const {
    return HANDLER::ProvidesImport(loc,type);
  }

  virtual bool ProvidesExport(const boost::filesystem::path& loc, const String& type) const {
    return HANDLER::ProvidesExport(loc,type);
  }

  virtual String GetFormatName() const {
     return HANDLER::GetFormatName();
  }

  virtual String GetFormatDescription() const {
     return HANDLER::GetFormatDescription();
  }

  virtual EntityIOHandlerP Create() const {
    return EntityIOHandlerP(new HANDLER);
  }
};



}} // ns

#endif
