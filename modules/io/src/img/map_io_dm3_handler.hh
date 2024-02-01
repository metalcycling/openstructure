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
#ifndef DX_IO_DM3_IO_PLUGIN_OST_HH
#define DX_IO_DM3_IO_PLUGIN_OST_HH

/*
  Authors: Ansgar Philippsen, Andreas Schenk
*/

#include "map_io_handler.hh"

namespace ost { namespace io {

class DLLEXPORT_OST_IO DM3: public ImageFormatBase
{
 public:

  DM3();
  static String FORMAT_STRING;
};

class DLLEXPORT_OST_IO MapIODm3Handler: public MapIOHandler
{
  public:
    /// \brief Map IO handler to read/write Dm3 map files
    ///
    /// This map IO handler reads and writes Dm3 formatted map files (Digital Micrograph).
    virtual void Import(img::MapHandle& sh, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct);
    virtual void Import(img::MapHandle& sh, std::istream& loc, const ImageFormatBase& formatstruct);
    virtual void Export(const img::MapHandle& sh, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct) const;
    virtual void Export(const img::MapHandle& sh, std::ostream& loc,const ImageFormatBase& formatstruct) const;
    static bool MatchContent(unsigned char* header);
    static bool MatchType(const ImageFormatBase& type);
    static bool MatchSuffix(const String& loc);
    static bool ProvidesImport() { return true; }
    static bool ProvidesExport() { return false; }
    static String GetFormatName() { return String("Dm3"); }
    static String GetFormatDescription() {return String("Format used by Gatan Inc.'s  Digital Micrograph software");}

};

typedef MapIOHandlerFactory<MapIODm3Handler> MapIODm3HandlerFactory;

}} // ns

#endif
