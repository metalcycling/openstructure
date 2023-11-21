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

#include <ost/io/mol/star_writer.hh>

namespace ost{ namespace io{

StarWriter::StarWriter(std::ostream& stream): filename_("<stream>"),
                                              file_open_(true) {
  if(!stream) {
    file_open_ = false;
  }
  stream_.push(stream);
}


StarWriter::StarWriter(const String& filename): filename_(filename),
                                                file_open_(true),
                                                fstream_(filename.c_str()) {
  if (!fstream_) {
    file_open_ = false;
  }
  stream_.push(fstream_);
}

void StarWriter::Write(const String& data_name) {
  if (!file_open_) {
    throw IOException("yolo");
  }

  // write data header
  stream_ << "data_" << data_name << std::endl;

  for(auto star_obj : categories_to_write_) {
    star_obj->ToStream(stream_);
    stream_ << String("#") << std::endl;
  }
}

}} // ns
