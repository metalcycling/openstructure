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
#ifndef OST_IO_MAP_IO_NANOSCOPE_HANDLER_HH
#define OST_IO_MAP_IO_NANOSCOPE_HANDLER_HH

#include "map_io_handler.hh"

namespace ost { namespace io {

class DLLEXPORT_OST_IO Nanoscope: public ImageFormatBase
{
 public:
  Nanoscope(int subimage = -1);

  static String FORMAT_STRING;
  static String FORMAT_NAME_STRING;
  static String FORMAT_DESCRIPTION_STRING;

  int GetSubimage () const;
  void SetSubimage (int subimage);

 private:

  int subimage_;

};



class MapIONanoscopeHandler: public MapIOHandler {
public:

   public:
  /// \brief Map IO handler to read/write Nanoscope map files
  ///
  /// This map IO handler reads and writes Nanoscope formatted map files
  virtual void Import(img::MapHandle& sh, const boost::filesystem::path& loc,const ImageFormatBase& formatstruct );
  virtual void Import(img::MapHandle& sh, std::istream& loc, const ImageFormatBase& formatstruct);
  virtual void Export(const img::MapHandle& sh, const boost::filesystem::path& loc, const ImageFormatBase& formatstruct) const;
  virtual void Export(const img::MapHandle& sh, std::ostream& loc,const ImageFormatBase& formatstruct) const;
  static bool MatchContent(unsigned char* header);
  static bool MatchType(const ImageFormatBase& type);
  static bool MatchSuffix(const String& loc);
  static bool ProvidesImport() { return true; }
  static bool ProvidesExport() { return true; }
  static String GetFormatName() { return String("Nanoscope"); }
  static String GetFormatDescription() { return String("Format used by software from Veeco"); }
};

typedef MapIOHandlerFactory<MapIONanoscopeHandler> MapIONanoscopeHandlerFactory;

}}

#endif
