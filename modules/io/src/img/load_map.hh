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
#ifndef OST_IO_LOAD_MAP_HH
#define OST_IO_LOAD_MAP_HH

#include <ost/io/module_config.hh>
#include <ost/img/map.hh>
#include <ost/img/alg/normalizer_factory.hh>
#include <ost/io/img/image_format.hh>

namespace ost { 

namespace io {
/// \brief Function that loads an image from a file
///
/// This function reads an image file from disk and transfers the data into memory
/// The function automatically determines the file format by looking at the filename
/// extension. In case the extension is not available, the content of the file will
/// be examined. This analysis will be sufficient to determine the format type for
/// many files. In case the file format cannot be determined, the function will throw an error.
DLLEXPORT_OST_IO img::ImageHandle LoadImage(const boost::filesystem::path& loc);

/// \brief Function that loads an image from a file with a specific format
///
/// This function reads an image file from disk and transfers the data into memory
/// The function assumes that the file has the specified format
DLLEXPORT_OST_IO img::ImageHandle LoadImage(const boost::filesystem::path& loc, const ImageFormatBase& formatstruct);

/// \brief Function that saves an image to a file
///
/// This function saves an image on disk. The function automatically determines the file format by looking at the filename
/// extension. In case the extension is not availabel, the function will throw an error
/// \sa \ref spoke_pattern.py "Spoke Pattern Example"
DLLEXPORT_OST_IO void SaveImage(const img::ImageHandle& image, const boost::filesystem::path& loc);
/// \brief Function that saves an image to a file using a specified format
///
/// This function saves an image on disk using the format specified in the function call
/// \sa \ref spoke_pattern.py "Spoke Pattern Example"
DLLEXPORT_OST_IO void SaveImage(const img::ImageHandle& image, const boost::filesystem::path& loc, const ImageFormatBase& formatstruct);

/// \brief Alias for the equivalent LoadImage function
DLLEXPORT_OST_IO img::MapHandle LoadMap(const boost::filesystem::path& loc);
/// \brief Alias for the equivalent LoadImage function
DLLEXPORT_OST_IO img::MapHandle LoadMap(const boost::filesystem::path& loc, const ImageFormatBase& formatstruct);
/// \brief Alias for the equivalent SaveImage function
DLLEXPORT_OST_IO void SaveMap(const img::ImageHandle& image, const boost::filesystem::path& loc);
/// \brief Alias for the equivalent SaveImage function
DLLEXPORT_OST_IO void SaveMap(const img::ImageHandle& image, const boost::filesystem::path& loc, const ImageFormatBase& formatstruct);
}} // ns

#endif
