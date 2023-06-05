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
#ifndef OST_IO_SDF_READER_HH
#define OST_IO_SDF_READER_HH

#include <tuple>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/filesystem/fstream.hpp>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/io/module_config.hh>

namespace ost { namespace io {




class DLLEXPORT_OST_IO SDFReader {
public:
  SDFReader(const String& filename);
  SDFReader(const boost::filesystem::path& loc);
  SDFReader(std::istream& instream);

  bool HasNext();

  void Import(mol::EntityHandle& ent);

private:
  typedef std::tuple<int, String, String, String, String, String> atom_data;
  typedef std::tuple<String, String, String> bond_data;
  typedef std::tuple<std::vector<String>, std::map<String, String>> v3000_line_tokens;

  void ClearState(const boost::filesystem::path& loc);
  void NextMolecule();

  void ParseHeader(const String& line, int line_num, mol::EntityHandle& ent,
                         mol::XCSEditor& editor);
  void SetCounts(const String& anum, const String bnum, int line_num);

  atom_data ParseAtom(const String& line, int line_num);
  void AddAtom(const atom_data& atom_tuple, int line_num, mol::EntityHandle& ent,
               bool hetatm, mol::XCSEditor& editor);

  bond_data ParseBond(const String& line, int line_num);
  void AddBond(const bond_data& bond_tuple, int line_num, mol::EntityHandle& ent,
                       mol::XCSEditor& editor);

  // V3000 methods
  v3000_line_tokens TokenizeV3000Line(const String& line, int line_num,
                                      int num_posval);
  String CleanupV3000Line(const String& line);
  void ProcessV3000Line(const String& line, mol::EntityHandle& ent,
                       mol::XCSEditor& editor);
  atom_data ParseV3000Atom(const String& line, int line_num);
  bond_data ParseV3000Bond(const String& line, int line_num);
  std::tuple<String, String> ParseV3000Counts(const String& line, int line_num);
  void VerifyV3000Counts();

  String curr_chain_name_;
  mol::ResidueKey curr_res_key_;
  mol::ChainHandle curr_chain_;
  mol::ResidueHandle curr_residue_;
  int chain_count_;
  int residue_count_;
  int atom_count_;
  int bond_count_;
  int line_num;
  boost::filesystem::ifstream infile_;
  std::istream& instream_;
  boost::iostreams::filtering_stream<boost::iostreams::input>  in_;
  String version_;
  bool v3000_atom_block_;
  bool v3000_bond_block_;
};

}}

#endif
