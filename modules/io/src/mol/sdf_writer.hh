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
/*
  Author: Tobias Schmidt
 */
#ifndef OST_IO_SDF_WRITER_HH
#define OST_IO_SDF_WRITER_HH

#include <iostream>
#include <sstream>
#include <iomanip>

#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <ost/log.hh>
#include <ost/conop/conop.hh>
#include <ost/io/io_exception.hh>
#include <ost/mol/entity_visitor.hh>

namespace ost { namespace io {

class DLLEXPORT_OST_IO SDFWriter : public mol::EntityViewVisitor {
public:
  SDFWriter(std::ostream& ostream);
  SDFWriter(const String& filename);
  SDFWriter(const boost::filesystem::path& filename);

  void Write(const mol::EntityView& ent);
  void Write(const mol::EntityHandle& ent);

private:
  virtual bool VisitChain(const mol::ChainView& chain);

  std::ofstream      outfile_;
  std::ostream&      ostr_;
  int                counter_;
  std::map<long,int> atom_indices_;
};

}}

#endif
