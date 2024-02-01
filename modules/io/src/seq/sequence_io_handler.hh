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
#ifndef OST_IO_SEQUENCE_IO_HANDLER_HH
#define OST_IO_SEQUENCE_IO_HANDLER_HH

#include <boost/shared_ptr.hpp>
#include <boost/filesystem/operations.hpp>


#include <ost/io/module_config.hh>
#include <ost/seq/sequence_list.hh>
#include <ost/io/io_utils.hh>

namespace ost { namespace io {

//! pure abstract base class for entity io handlers
class DLLEXPORT_OST_IO SequenceIOHandler {
public:
  virtual ~SequenceIOHandler() {}


  virtual void Import(seq::SequenceList& sequences,
                      const boost::filesystem::path& loc) = 0;

  virtual void Export(const seq::ConstSequenceList& ent,
                      const boost::filesystem::path& loc) const = 0;
                      
  virtual void Import(seq::SequenceList& ent,
                      std::istream& instream)=0;

  virtual void Export(const seq::ConstSequenceList& ent,
                      std::ostream& ostream) const=0;
};

typedef boost::shared_ptr<SequenceIOHandler> SequenceIOHandlerPtr;


//! pure abstract base class for creation of a specfic alignment io handler
class DLLEXPORT_OST_IO SequenceIOHandlerFactoryBase {
public:
  virtual ~SequenceIOHandlerFactoryBase() {}
  virtual bool ProvidesImport(const boost::filesystem::path& loc, 
                              const String& type) const = 0;
  virtual bool ProvidesExport(const boost::filesystem::path& loc, 
                              const String& type) const = 0;
  virtual SequenceIOHandlerPtr Create() const = 0;
  virtual String GetFormatName() const =0;
  virtual String GetFormatDescription() const =0;
};

typedef boost::shared_ptr<SequenceIOHandlerFactoryBase> SequenceIOHandlerFactoryBasePtr;

template <class HANDLER>
class SequenceIOHandlerFactory: public SequenceIOHandlerFactoryBase {
  virtual bool ProvidesImport(const boost::filesystem::path& loc, 
                              const String& type) const {
    return HANDLER::ProvidesImport(loc,type);
  }

  virtual bool ProvidesExport(const boost::filesystem::path& loc, 
                              const String& type) const {
    return HANDLER::ProvidesExport(loc,type);
  }

  virtual SequenceIOHandlerPtr Create() const {
    return SequenceIOHandlerPtr(new HANDLER);
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
