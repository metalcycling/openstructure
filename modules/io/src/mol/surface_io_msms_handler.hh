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
#ifndef OST_IO_SURFACE_IO_PLUGIN_MSMS_HH
#define OST_IO_SURFACE_IO_PLUGIN_MSMS_HH

#include "surface_io_handler.hh"

namespace ost { namespace io {

class DLLEXPORT_OST_IO SurfaceIOMSMSHandler: public SurfaceIOHandler
{
 public:
  virtual void Import(mol::SurfaceHandle& sh, const boost::filesystem::path& loc);
  virtual void Export(const mol::SurfaceHandle& sh,
                      const boost::filesystem::path& loc) const;

  static bool ProvidesImport(const boost::filesystem::path& loc, 
                             const String& format="auto");
  static bool ProvidesExport(const boost::filesystem::path& loc, 
                             const String& format="auto");
  static String GetFormatName() { return String("Msms"); }
  static String GetFormatDescription() { return String("Format used by the MSMS program"); }
};

typedef SurfaceIOHandlerFactory<SurfaceIOMSMSHandler> SurfaceIOMSMSHandlerFactory;

}} // ns

#endif
