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
  Author: Ansgar Philippsen
*/

#ifndef OST_IO_ENTITY_IO_PLUGIN_PQR_H
#define OST_IO_ENTITY_IO_PLUGIN_PQR_H

#include <ost/mol/entity_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/entity_visitor.hh>
#include <ost/io/mol/entity_io_handler.hh>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem/fstream.hpp>

namespace ost { namespace io {

/// \brief CHARMM coordinate file import
class DLLEXPORT_OST_IO PQRReader {
public:
  PQRReader(const boost::filesystem::path& loc);
  
  void Import(mol::EntityHandle& ent);

  std::vector<mol::AtomHandle> GetSequentialAtoms() const;

private:

  void ParseAndAddAtom(const String& line, mol::EntityHandle& h);
  void ParseAndAddAtomExpanded(const String& line, mol::EntityHandle& h);

  std::vector<mol::AtomHandle> sequential_atom_list_;
  mol::ChainHandle curr_chain_;
  mol::ResidueHandle curr_residue_;
  int chain_count_;
  int residue_count_;
  int atom_count_;
  boost::filesystem::ifstream infile_;
  boost::iostreams::filtering_stream<boost::iostreams::input>  in_;
};

/// \brief CHARMM coordinate file export
class DLLEXPORT_OST_IO PQRWriter : public mol::EntityVisitor {
public:
  PQRWriter(const String& filename, bool ext=false);
  PQRWriter(const boost::filesystem::path& filename, bool ext=false);
  PQRWriter(std::ostream& outstream, bool ext=false);

  void Write(const mol::EntityView& ent);
  void Write(const mol::EntityHandle& ent);

  virtual bool VisitAtom(const mol::AtomHandle& atom);
  virtual bool VisitResidue(const mol::ResidueHandle& r);

  void WriteHeader(const mol::EntityView& ent);

private:
  void Init();

  std::ofstream   outfile_;
  std::ostream&   outstream_;
  bool ext_;
  int atom_count_;
  int atom_total_;
  int res_count_;
};

class DLLEXPORT_OST_IO EntityIOPQRHandler: public EntityIOHandler {
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
  virtual bool RequiresBuilder() const;

  static String GetFormatName() { return String("Crd"); }
  static String GetFormatDescription() { return String("CARD format file used by the Charmm software package"); }
};


typedef EntityIOHandlerFactory<EntityIOPQRHandler> EntityIOPQRHandlerFactory;

mol::EntityHandle DLLEXPORT_OST_IO LoadPQR(const String& file_name);

}} // ns

#endif
