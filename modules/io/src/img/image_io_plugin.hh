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

/*
  IO plugin interface


  Authors: Ansgar Philippsen, Andreas Schenk
*/

#ifndef IMG_IMAGE_IO_PLUGIN
#define IMG_IMAGE_IO_PLUGIN

#include <boost/filesystem/operations.hpp>
#include <boost/shared_ptr.hpp>

#include <img/base.hh>
#include <img/image.hh>
#include <img/normalizer_impl.hh>
#include <img/message.hh>
#include "image_format.hh"

namespace ost { namespace img { namespace io {

using boost::filesystem::path;

//!  abstract base class to derive image IO plugins from
class ImageIOPlugin {
public:
  // no explicit constructor
  // but virtual destructor
  virtual ~ImageIOPlugin() {}

  /*!
    attempt to detect wether the given file and/or type matches
    the format implemented by this plugin, return true if so, 
    false otherwise
  */
  virtual bool MatchContent(unsigned  char* header) = 0;
  virtual bool MatchType(const string& type) = 0;
  virtual bool MatchSuffix(const string& suffix) = 0;

  //! import from the given file into the given Image
  virtual bool Import(const path& location, ImageHandle& i, const ImageFormat& formatstruct) = 0;

  //! export from the given Image into the given file
  virtual bool Export(const path& location, const ConstImageHandle& i, const ImageFormat& formatstruct,const  NormalizerPtr& nptr) = 0;

  virtual Real GetMinimum(const ImageFormat& formatstruct) const = 0;
  virtual Real GetMaximum(const ImageFormat& formatstruct) const = 0;
};

// types of factories
typedef ImageIOPlugin* image_io_plugin_create();
typedef void image_io_plugin_destroy(ImageIOPlugin*);

typedef boost::shared_ptr<ImageIOPlugin> ImageIOPluginPtr;

}}} // namespaces

#endif
