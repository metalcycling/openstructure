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
#ifndef OST_IO_PROMOD_IO_HANDLER_HH
#define OST_IO_PROMOD_IO_HANDLER_HH

/*
  Author: Marco Biasini
 */

#include <ost/io/seq/sequence_io_handler.hh>

namespace ost { namespace io {

class DLLEXPORT_OST_IO PromodIOHandler : public SequenceIOHandler{
public:
  virtual void Import(seq::SequenceList& ali,
                      const boost::filesystem::path& loc);

  virtual void Export(const seq::ConstSequenceList& ent,
                      const boost::filesystem::path& loc) const;
                      
  virtual void Import(seq::SequenceList& ent,
                      std::istream& instream);

  virtual void Export(const seq::ConstSequenceList& ent,
                      std::ostream& ostream) const;

  static bool ProvidesImport(const boost::filesystem::path& loc, 
                             const String& format="auto");
  static bool ProvidesExport(const boost::filesystem::path& loc, 
                             const String& format="auto");

  static String GetFormatName() { return String("ProMod"); }
  static String GetFormatDescription() { return String("Format used by the SwissModel modeling server"); }
};

typedef SequenceIOHandlerFactory<PromodIOHandler> PromodIOHandlerFactory; 

}}

#endif
