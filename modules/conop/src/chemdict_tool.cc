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
  Author: Marco Biasini
 */

#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#include <ost/io/mol/chemdict_parser.hh>

using namespace ost;

void PrintUsage()
{
  std::cout << "usage: chemdict_tool ACTION <compound-dict> <db> [DIALECT] [OPTIONS]" << std::endl;
  std::cout << "supported actions are:" << std::endl;
  std::cout << "  create  - creates a new db " << std::endl;
  std::cout << "  update  - update existing db" << std::endl;
  std::cout << "supported dialects are: pdb, charmm, amber, opls" << std::endl;
  std::cout << "supported options are:" << std::endl;
  std::cout << "  -i  - ignore compounds reserved by the PDB (01-99, DRG, INH, LIG)" << std::endl;
}

int main(int argc, char const *argv[])
{
  if (argc!=4 && argc!=5 && argc!=6) {
    PrintUsage();
    return 0;
  }
  conop::Compound::Dialect dialect=conop::Compound::PDB;
  bool ignore_reserved=false;
  for (int i = 4; i < argc; i++) {
    String param=argv[i];

    if (param=="charmm") {
      dialect=conop::Compound::CHARMM;
    } else if (param=="pdb") {
      dialect=conop::Compound::PDB;
    } else if (param=="opls") {
      dialect=conop::Compound::OPLS;
    } else if (param=="amber") {
      dialect=conop::Compound::AMBER;
    } else if (param=="-i") {
      ignore_reserved=true;
    } else {
      PrintUsage();
      return 0;
    }
  }
  boost::iostreams::filtering_stream<boost::iostreams::input>  filtered_istream;  
  std::ifstream istream(argv[2]);
  if (! istream.is_open()) {
      std::cout << "Cannot open " << argv[2] << ": [Errno " << errno << "] "
                << strerror(errno) << std::endl;
      return 1;
  }
  if (boost::iequals(".gz", boost::filesystem::extension(argv[2]))) {
    filtered_istream.push(boost::iostreams::gzip_decompressor());
  }
  filtered_istream.push(istream);  
  io::ChemdictParser cdp(filtered_istream, dialect, ignore_reserved);
  conop::CompoundLibPtr compound_lib;
  bool in_mem=false;
  if (!strncmp(argv[1], "create", 6)) {
    compound_lib=conop::CompoundLib::Create(":memory:");
    in_mem=true;
  } else if (!strncmp(argv[1], "update", 6)) {
    compound_lib=conop::CompoundLib::Load(argv[3]);
  } else {
    PrintUsage();
    return 0;
  }
  if (!compound_lib) {
    return 0;
  }
  assert(compound_lib);
  conop::CompoundLibPtr in_mem_lib=in_mem ? compound_lib :
                                   compound_lib->Copy(":memory:");
  compound_lib.reset();
  cdp.SetCompoundLib(in_mem_lib);
  cdp.Parse();
  in_mem_lib->SetChemLibInfo();
  conop::CompoundLibPtr copy = in_mem_lib->Copy(argv[3]);
  if (! copy) {
      std::cout << "Cannot save " << argv[3] << ": [Errno " << errno << "] "
                << strerror(errno) << std::endl;
      return 1;
  }
  return 0;
}
