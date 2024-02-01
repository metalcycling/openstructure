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
#ifndef OST_IO_ENTITY_IO_PLUGIN_MAE_H
#define OST_IO_ENTITY_IO_PLUGIN_MAE_H

#include <ost/mol/entity_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/io/mol/entity_io_handler.hh>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem/fstream.hpp>

namespace ost { namespace io {

class DLLEXPORT_OST_IO MAEReader {
public:
  MAEReader(const boost::filesystem::path& loc);
  
  void Import(mol::EntityHandle& ent);

private:

  void parse_and_add_atom(mol::EntityHandle ent,
                          mol::XCSEditor& editor,
                          char* line,
                          size_t line_len,
                          int i_atom_name,
                          int i_atom_xpos,
                          int i_atom_ypos,
                          int i_atom_zpos,
                          int i_res_name,
                          int i_res_num,
                          int i_chain_name);

  mol::ChainHandle curr_chain_;
  mol::ResidueHandle curr_residue_;
  int chain_count_;
  int residue_count_;
  int atom_count_;
  boost::filesystem::ifstream infile_;
  boost::iostreams::filtering_stream<boost::iostreams::input>  in_;
};

class DLLEXPORT_OST_IO EntityIOMAEHandler: public EntityIOHandler {
public:
  virtual void Import(mol::EntityHandle& ent, const boost::filesystem::path& loc);
  
  virtual void Export(const mol::EntityView& ent, 
                      const boost::filesystem::path& loc) const;
                      
  virtual void Import(mol::EntityHandle& ent, std::istream& stream);

  virtual void Export(const mol::EntityView& ent, std::ostream& stream) const;
  
  static bool ProvidesImport(const boost::filesystem::path& loc, 
                             const String& format="auto");
  static bool ProvidesExport(const boost::filesystem::path& loc, 
                             const String& format="auto");
  virtual bool RequiresProcessor() const;

  static String GetFormatName() { return String("Mae"); }
  static String GetFormatDescription() { return String("MAEstro coordinate file format"); }
};


typedef EntityIOHandlerFactory<EntityIOMAEHandler> EntityIOMAEHandlerFactory;

mol::EntityHandle DLLEXPORT_OST_IO LoadMAE(const String& file_name);

}} // ns

#endif
