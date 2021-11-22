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
#ifndef OST_IO_MANAGER_H
#define OST_IO_MANAGER_H

#include <vector>

#include <ost/config.hh>
#include <ost/io/module_config.hh>

#include <ost/io/mol/entity_io_handler.hh>
#include <ost/io/seq/sequence_io_handler.hh>
#include <ost/io/seq/profile_io_handler.hh>
#include <ost/io/mol/surface_io_handler.hh>
#include <ost/io/img/map_io_handler.hh>

#include <ost/io/io_exception.hh>

namespace ost { namespace io {

typedef std::vector<EntityIOHandlerFactoryBaseP> EntityIOHFList;
typedef std::vector<SequenceIOHandlerFactoryBasePtr> AlignmentIOFList;
typedef std::vector<ProfileIOHandlerFactoryBasePtr> ProfileIOFList;
typedef std::vector<SurfaceIOHandlerFactoryBasePtr> SurfaceIOFList;
typedef std::vector<MapIOHandlerFactoryBasePtr> MapIOFList;  

/// \brief Central registry for input/output handlers
class DLLEXPORT_OST_IO IOManager {
public:
  /// \brief get entity io handler that is able to import the given file.
  ///
  /// Find EntityIOHandler that provides import for the given file. The method
  /// will return the first EntityIOHandler in chain that is able to import the
  /// file either bases on the \var filename, or, if provided, by \var format.
  ///
  /// \param filename  The filename of the file to export. Most IOHandlers 
  ///    will use the file extension to find out if the can import the file but 
  ///    may also implement more elaborate methods such as magic number matching.
  /// \param format If set to any other value than auto, the IO handlers will 
  ///    treat \c filename as if \c format were the file extension (without 
  ///    the dot).
  /// \return The first EntityIOHandler in chain that provides import or a 
  ///    null pointer if no importer is able to load the  file.
  /// \sa FindEntityExportHandler(const String&,const String&)
  EntityIOHandlerP FindEntityImportHandler(const String& filename, 
                                           const String& format="auto");
                                           
  /// \brief get entity IO handler that is able to export the file.
  ///
  /// \note The same rules as for 
  ///  FindEntityImportHandler(const String&, const String&) apply here.
  EntityIOHandlerP FindEntityExportHandler(const String& filename,
                                           const String& format="auto");
      

  /// \name Alignment/Sequence IO
  //@{
  /// \brief get sequence io handler that is able to import the given file.
  /// \sa FindEntityImportHandler(filename)       
  SequenceIOHandlerPtr FindAlignmentImportHandler(const String& filename,
                                                  const String& format="auto");
                                                   
  SequenceIOHandlerPtr FindAlignmentExportHandler(const String& filename,
                                                  const String& format="auto");

  ProfileIOHandlerPtr FindProfileImportHandler(const String& filename,
                                               const String& format="auto");
                                                   
  ProfileIOHandlerPtr FindProfileExportHandler(const String& filename,
                                               const String& format="auto");

  SurfaceIOHandlerPtr FindSurfaceImportHandler(const String& filename,
                                               const String& format="auto");
  //@}
  /// \brief register entity io handler factory
  /// \sa adding_io_handler
  void RegisterFactory(const EntityIOHandlerFactoryBaseP&);
  /// \brief register aligment io handler factory
  /// \sa adding_io_handler
  void RegisterFactory(const SequenceIOHandlerFactoryBasePtr&);
  /// \brief register profile io handler factory
  /// \sa adding_io_handler
  void RegisterFactory(const ProfileIOHandlerFactoryBasePtr&);
  /// \brief register surface io handler factory
  /// \sa adding_io_handler
  void RegisterFactory(const SurfaceIOHandlerFactoryBasePtr&);  

  /// \brief Get a list with all available EntityHandler
  const EntityIOHFList& GetAvailableEntityHandler() const;

  /// \brief Get a list with all available AlignmentHandler
  const AlignmentIOFList& GetAvailableAlignmentHandler() const;

  /// \brief Get a list with all available ProfileHandler
  const ProfileIOFList& GetAvailableProfileHandler() const;

  /// \brief Get a list with all available SurfaceHandler
  const SurfaceIOFList& GetAvailableSurfaceHandler() const;

  /// \name Image/Map IO
  //@{
  /// \brief register map io handle factory
  /// \sa adding_io_handler
  void RegisterFactory(const MapIOHandlerFactoryBasePtr&);  

  MapIOHandlerPtr FindMapImportHandlerFile(const boost::filesystem::path& loc,
                                           const ImageFormatBase& format);

  MapIOHandlerPtr FindMapImportHandlerStream(std::istream& stream,
                                             const ImageFormatBase& format);

  MapIOHandlerPtr FindMapExportHandlerFile(const boost::filesystem::path& loc,
                                           const ImageFormatBase& format);

  MapIOHandlerPtr FindMapExportHandlerStream(std::istream& stream,
                                             const ImageFormatBase& format);

  const MapIOFList& GetAvailableMapHandler() const;
  //@}

  // singleton interface
  static IOManager& Instance();

private:
  IOManager();
  IOManager(const IOManager&) {}
  IOManager& operator=(const IOManager&) {return *this;}

  EntityIOHFList entity_iohf_list_;
  AlignmentIOFList alignment_io_list_;
  ProfileIOFList profile_io_list_;
  SurfaceIOFList surface_io_list_;
  MapIOFList map_io_list_;
};



}}

#endif
