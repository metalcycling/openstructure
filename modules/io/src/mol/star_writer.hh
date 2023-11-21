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
#ifndef OST_IO_STAR_WRITER_HH
#define OST_IO_STAR_WRITER_HH

#include <map>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <ost/string_ref.hh>
#include <ost/io/mol/star_base.hh>


namespace ost { namespace io {

class DLLEXPORT_OST_IO StarWriter {
public:
  StarWriter(std::ostream& stream);
  StarWriter(const String& filename);
  virtual ~StarWriter() { }

  void Push(StarObject* obj) { categories_to_write_.push_back(obj); }
  void Write(const String& data_name);
private:
  String filename_;
  bool file_open_;
  std::ofstream fstream_;
  boost::iostreams::filtering_stream<boost::iostreams::output> stream_;
  std::vector<StarObject*> categories_to_write_;
};

}} // ns

#endif
