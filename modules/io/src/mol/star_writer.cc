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

void StarWriter::Write(const String& data_name, std::ostream& stream) {
  if(!stream) {
      throw IOException("[Errno " + std::to_string(errno) + "] " +
                        std::string(strerror(errno)) +
                        ": <stream>");
  }
  // write data header
  stream << "data_" << data_name << std::endl;
  // write StarWriterObjects
  for(auto star_obj : categories_to_write_) {
    star_obj->ToStream(stream);
    stream << String("#") << std::endl;
  }
}


void StarWriter::Write(const String& data_name, const String& filename) {
  std::ofstream fstream(filename.c_str());
  if (!fstream) {
      throw IOException("[Errno " + std::to_string(errno) + "] " +
                        std::string(strerror(errno)) +
                        ": '" + filename + "'");
  }
  boost::iostreams::filtering_stream<boost::iostreams::output> stream;
  if (boost::iequals(".gz", boost::filesystem::extension(filename))) {
    stream.push(boost::iostreams::gzip_compressor());
  }
  stream.push(fstream);
  this->Write(data_name, stream);
}

}} // ns
