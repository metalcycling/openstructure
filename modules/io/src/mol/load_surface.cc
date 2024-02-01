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
#include <ost/log.hh>
#include <ost/mol/surface_handle.hh>

#include "load_surface.hh"
#include <ost/io/io_manager.hh>
#include "surface_io_handler.hh"

namespace ost { namespace io {

namespace {

void ImportSurface(mol::SurfaceHandle& sh, const String& fname, const String& type)
{
  LOG_DEBUG("creating EntityIOHandle for " << fname);
  SurfaceIOHandlerPtr surf_io = IOManager::Instance().FindSurfaceImportHandler(fname,type);
  LOG_DEBUG("calling import on surface io handle");
  surf_io->Import(sh,fname);
}

}

mol::SurfaceHandle LoadSurface(const String& filename, const String& type)
{
  mol::SurfaceHandle sh = mol::CreateSurface();
  ImportSurface(sh,filename,type);
  return sh;
}

mol::SurfaceHandle LoadManagedSurface(const String& name,
                                 const String& filename, const String& type)
{
  mol::SurfaceHandle sh = mol::CreateSurface();
  ImportSurface(sh,filename,type);
  return sh;
}

}} // ns
