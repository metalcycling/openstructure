//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2023 by the OpenStructure authors
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

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem.hpp>

#include <ost/io/mol/star_writer.hh>

namespace ost{ namespace io{

StarWriter::StarWriter(std::ostream& stream): filename_("<stream>") {
  if(!stream) {
    std::stringstream ss;
    ss << "Cannot open stream: [Errno " << errno << "] "
       << strerror(errno) << std::endl;
    throw IOException(ss.str());
  }
  stream_.push(stream);
}


StarWriter::StarWriter(const String& filename): filename_(filename),
                                                fstream_(filename.c_str()) {
  if (!fstream_) {
    std::stringstream ss;
    ss << "Cannot open " << filename_ << ": [Errno " << errno << "] "
       << strerror(errno) << std::endl;
    throw IOException(ss.str());
  }
  if (boost::iequals(".gz", boost::filesystem::extension(filename))) {
    stream_.push(boost::iostreams::gzip_compressor());
  }
  stream_.push(fstream_);
}

void StarWriter::Write(const String& data_name) {
  // write data header
  stream_ << "data_" << data_name << std::endl;

  for(auto star_obj : categories_to_write_) {
    star_obj->ToStream(stream_);
    stream_ << String("#") << std::endl;
  }
}

}} // ns
