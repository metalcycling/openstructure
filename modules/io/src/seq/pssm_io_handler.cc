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
  Author: Gerardo Tauriello
 */

#include "pssm_io_handler.hh"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <ost/string_ref.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/io_utils.hh>

namespace ost { namespace io {

void PssmIOHandler::Import(seq::ProfileHandle& prof,
                           const boost::filesystem::path& loc) {
  // open it up
  boost::iostreams::filtering_stream<boost::iostreams::input> in;
  boost::filesystem::ifstream stream(loc);
  if (!stream) {
    throw IOException("Could not open " + loc.string());
  }
  // add unzip if necessary
  if (boost::iequals(".gz", boost::filesystem::extension(loc))) {
    in.push(boost::iostreams::gzip_decompressor());
  }
  in.push(stream);

  this->Import(prof, in);
}

void PssmIOHandler::ImportFromString(seq::ProfileHandle& prof,
                                     const String& data) {
  std::stringstream ss(data);
  this->Import(prof, ss);
}

void PssmIOHandler::Export(const seq::ProfileHandle& prof,
                          const boost::filesystem::path& loc) const {
  throw IOException("Cannot write pssm files.");
}

bool PssmIOHandler::ProvidesImport(const boost::filesystem::path& loc, 
                                  const String& format) {
  if (format=="auto") {
    String match_suf_string = loc.string();
    std::transform(match_suf_string.begin(), match_suf_string.end(),
                   match_suf_string.begin(), tolower);
    if (   detail::FilenameEndsWith(match_suf_string,".pssm")
        || detail::FilenameEndsWith(match_suf_string,".pssm.gz")) {
      return true;
    }
  } else if (format == "pssm") {
    return true;
  }
  return false;
}

bool PssmIOHandler::ProvidesExport(const boost::filesystem::path& loc, 
                                  const String& format) {
  // no writers here
  return false;
}

void PssmIOHandler::Import(seq::ProfileHandle& prof,
                           std::istream& in) {

  // reset profile
  prof.clear();

  // set null model
  prof.SetNullModel(seq::ProfileColumn::BLOSUMNullModel());

  // tmp. storage
  std::string line;
  ost::StringRef sline;
  char olcs[20];
  std::vector<ost::StringRef> chunks;

  // read line-by-line looking for table
  // -> start = line with 40 1-char tokens
  bool table_found = false;
  while (std::getline(in, line)) {
    sline = ost::StringRef(line.c_str(), line.length());
    chunks = sline.split();
    if (chunks.size() == 40) {
      bool got_it = true;
      for (uint i = 0; i < 40; ++i) {
        if (chunks[i].length() != 1) {
          got_it = false;
          break;
        }
        if (i >= 20) {
          olcs[i-20] = chunks[i][0];
        }
      }
      if (got_it) {
        table_found = true;
        break;
      }
    }
  }
  if (!table_found) {
    throw IOException("No ASCII table found in input");
  }

  // parse table (assume: index olc 20xscore 20xfreq)
  while (std::getline(in, line)) {
    sline = ost::StringRef(line.c_str(), line.length());
    chunks = sline.split();
    // stop if not enough entries
    if (chunks.size() < 42) break;
    // parse row (freq in % as integers)
    const char olc = chunks[1][0];
    Real sum_freq = 0;
    Real freqs[20];
    for (uint i = 22; i < 42; ++i) {
      std::pair<bool, int> pbi = chunks[i].to_int();
      if (!pbi.first) {
        throw IOException("Badly formatted line\n" + line);
      }
      sum_freq += pbi.second;
      freqs[i-22] = pbi.second;
    }
    // fill normalized frequencies
    seq::ProfileColumn pc;
    for (uint i = 0; i < 20; ++i) {
      pc.SetFreq(olcs[i], freqs[i]/sum_freq);
    }
    prof.AddColumn(pc, olc);
  }
}

}}
