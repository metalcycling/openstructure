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
#ifndef OST_GUI_FILE_LOADER_HH
#define OST_GUI_FILE_LOADER_HH
#include <vector>


#include <ost/gfx/gfx_object.hh>

#include <ost/gui/module_config.hh>

#include <ost/io/io_exception.hh>
#include <ost/io/mol/entity_io_handler.hh>
#include <ost/io/seq/sequence_io_handler.hh>
#include <ost/io/mol/surface_io_handler.hh>
#include <ost/io/img/map_io_handler.hh>

#include <QString>
#include <QList>
#include <QMap>
namespace ost { namespace gui {

class DLLEXPORT_OST_GUI FileLoader {
private:

  enum ErrorType { DEFAULT = 0,
    IO_LOADING,
    GFX_ADD,
    GFX_MULTIPLE_ADD,
    INFO
  };

  FileLoader();
  static gfx::GfxObjP TryLoadEntity(const QString& filename, io::EntityIOHandlerP handler=io::EntityIOHandlerP(), const QString& selection=QString());
  static gfx::GfxObjP TryLoadSurface(const QString& filename, io::SurfaceIOHandlerPtr handler=io::SurfaceIOHandlerPtr());
  static gfx::GfxObjP TryLoadAlignment(const QString& filename, io::SequenceIOHandlerPtr handler=io::SequenceIOHandlerPtr());
  static gfx::GfxObjP TryLoadMap(const QString& filename, io::MapIOHandlerPtr handler=io::MapIOHandlerPtr());
  static void RunScript(const QString& filename);
  static void LoadPDB(const QString& filename, const QString& selection=QString());
  static void AddToScene(const QString& filename, gfx::GfxObjP obj);
  static void HandleError(const Error& e, ErrorType type, const QString& filename, gfx::GfxObjP obj=gfx::GfxObjP());
  static gfx::GfxObjP NoHandlerFound(const QString& filename);
  virtual ~FileLoader();
  static QList<img::ImageHandle> loaded_images_;

public:
  static void LoadObject(const QString& filename, const QString& selection=QString());
  static std::vector<String> GetSiteLoaderIdents();
};

} }

#endif /* OST_GUI_FILE_LOADER_HH */
