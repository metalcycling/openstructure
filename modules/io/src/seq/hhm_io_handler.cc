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
  Author: Gerardo Tauriello, Gabriel Studer
 */
  
#include "hhm_io_handler.hh"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <ost/string_ref.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/io_utils.hh>

namespace ost { namespace io {

namespace {
// reset pc with freqs in chunks[i+offset] (i in [0,20[)
// checks for consistency of chunks and values read
seq::ProfileColumn GetColumn(const std::vector<ost::StringRef>& chunks,
                             const uint offset, const char* olcs,
                             const String& exception) {
  // check chunks size
  if (chunks.size() < offset+20) {
    throw IOException(exception);
  }
  // reset pc
  seq::ProfileColumn pc; // all 0
  for (uint i = 0; i < 20; ++i) {
    if (chunks[i+offset] != ost::StringRef("*", 1)) {
      std::pair<bool, int> pbi = chunks[i+offset].to_int();
      if (!pbi.first) {
        throw IOException(exception);
      }
      pc.SetFreq(olcs[i], pow(2.0, -0.001 * pbi.second));
    }
  }
  return pc;
}

seq::HMMDataPtr GetHMMData(const std::vector<ost::StringRef>& chunks, 
                           const std::vector<uint>& idx, 
                           const std::vector<seq::HMMTransition>& transitions,
                           int neff_idx,
                           int neff_i_idx,
                           int neff_d_idx, 
                           const String& exception) {

  // check chunks size
  if(chunks.size() != 10) {
    throw IOException(exception);
  }

  seq::HMMDataPtr data(new seq::HMMData);
  for (uint i = 0; i < idx.size(); ++i) {
    if (chunks[idx[i]] != ost::StringRef("*", 1)) {
      std::pair<bool, int> p = chunks[idx[i]].to_int();
      if (!p.first) {
        throw IOException(exception);
      }
      data->SetProb(transitions[i], pow(2.0, -0.001 * p.second));
    } else {
      data->SetProb(transitions[i], 0.0);
    }
  }

  if (chunks[neff_idx] != ost::StringRef("*", 1)) {
    std::pair<bool, int> p = chunks[neff_idx].to_int();
    if (!p.first) {
      throw IOException(exception);
    }
    data->SetNeff(0.001 * p.second);
  } else {
    data->SetNeff(0.0);
  }

  if (chunks[neff_i_idx] != ost::StringRef("*", 1)) {
    std::pair<bool, int> p = chunks[neff_i_idx].to_int();
    if (!p.first) {
      throw IOException(exception);
    }
    data->SetNeff_I(0.001 * p.second);
  } else {
    data->SetNeff_I(0.0);
  }

  if (chunks[neff_d_idx] != ost::StringRef("*", 1)) {
    std::pair<bool, int> p = chunks[neff_d_idx].to_int();
    if (!p.first) {
      throw IOException(exception);
    }
    data->SetNeff_D(0.001 * p.second);
  } else {
    data->SetNeff_D(0.0);
  }

  return data;
}


} // anon ns

void HhmIOHandler::Import(seq::ProfileHandle& prof,
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

void HhmIOHandler::ImportFromString(seq::ProfileHandle& prof, 
                                    const String& data) {
  std::stringstream ss(data);
  this->Import(prof, ss);
}

void HhmIOHandler::Export(const seq::ProfileHandle& prof,
                          const boost::filesystem::path& loc) const {
  throw IOException("Cannot write hhm files.");
}

bool HhmIOHandler::ProvidesImport(const boost::filesystem::path& loc, 
                                  const String& format) {
  if (format=="auto") {
    String match_suf_string = loc.string();
    std::transform(match_suf_string.begin(), match_suf_string.end(),
                   match_suf_string.begin(), tolower);
    if (   detail::FilenameEndsWith(match_suf_string,".hhm")
        || detail::FilenameEndsWith(match_suf_string,".hhm.gz")) {
      return true;
    }
  } else if (format == "hhm") {
    return true;
  }
  return false;
}

bool HhmIOHandler::ProvidesExport(const boost::filesystem::path& loc, 
                                  const String& format) {
  // no writers here
  return false;
}

void HhmIOHandler::Import(seq::ProfileHandle& prof, std::istream& in) {
  // reset profile
  prof.clear();
  
  // tmp. storage
  std::string line;
  ost::StringRef sline;
  std::string null_line;
  std::string neff_line;
  char olcs[20];
  std::vector<ost::seq::HMMTransition> transitions;
  std::vector<uint> transition_idx; 
  int neff_idx = -1;
  int neff_i_idx = -1;
  int neff_d_idx = -1;
  std::vector<ost::StringRef> chunks;

  // read line-by-line looking for NEFF and NULL
  bool null_found = false;
  bool neff_found = false;
  while (std::getline(in, line)) {
    sline = ost::StringRef(line.c_str(), line.length());
    if (sline.length()>4 && 
        sline.substr(0, 5) == ost::StringRef("NULL ", 5)) {
      null_line = line;
      null_found = true;
      break;
    }
    if (sline.length()>4 && 
        sline.substr(0, 5) == ost::StringRef("NEFF ", 5)) {
      neff_line = line;
      neff_found = true;
    }
  }
  if (!null_found) {
    throw IOException("No NULL found in file");
  }

  // read until we hit HMM, extract olcs, then skip 2 more lines
  bool hmm_found = false;
  while (std::getline(in, line)) {
    sline = ost::StringRef(line.c_str(), line.length());
    chunks = sline.split();
    if (chunks[0] == ost::StringRef("HMM", 3)) {
      // extract olcs
      if (chunks.size() != 21) {
        throw IOException("Badly formatted HMM line: " + line);
      }
      for (uint i = 1; i < 21; ++i) {
        if (chunks[i].length() != 1) {
          throw IOException("Badly formatted HMM line: " + line);
        }
        olcs[i-1] = chunks[i][0];
      }
      // extract hmm data
      std::getline(in, line);
      sline = ost::StringRef(line.c_str(), line.length());
      chunks = sline.split();
      if (chunks.size() != 10) {
        throw IOException("Badly formatted HMM line: " + line);
      }

      for (uint i = 0; i < 10; ++i) {
        if(chunks[i] == ost::StringRef("M->M", 4)) {
          transition_idx.push_back(i);
          transitions.push_back(ost::seq::HMM_M2M);  
        } else if(chunks[i] == ost::StringRef("M->I", 4)) {
          transition_idx.push_back(i);
          transitions.push_back(ost::seq::HMM_M2I); 
        } else if(chunks[i] == ost::StringRef("M->D", 4)) {
          transition_idx.push_back(i);
          transitions.push_back(ost::seq::HMM_M2D); 
        } else if(chunks[i] == ost::StringRef("I->M", 4)) {
          transition_idx.push_back(i);
          transitions.push_back(ost::seq::HMM_I2M); 
        } else if(chunks[i] == ost::StringRef("I->I", 4)) {
          transition_idx.push_back(i);
          transitions.push_back(ost::seq::HMM_I2I); 
        } else if(chunks[i] == ost::StringRef("D->M", 4)) {
          transition_idx.push_back(i);
          transitions.push_back(ost::seq::HMM_D2M); 
        } else if(chunks[i] == ost::StringRef("D->D", 4)) {
          transition_idx.push_back(i);
          transitions.push_back(ost::seq::HMM_D2D); 
        } else if(chunks[i] == ost::StringRef("Neff", 4)) {
          neff_idx = i;
        } else if(chunks[i] == ost::StringRef("Neff_I", 6)) {
          neff_i_idx = i;
        } else if(chunks[i] == ost::StringRef("Neff_D", 6)) {
          neff_d_idx = i;
        } else {
          throw IOException("Badly formatted HMM line: " + line); 
        }
      }

      // check whether we really got everything
      if(neff_idx == -1 || neff_i_idx == -1 || neff_d_idx == -1) {
        throw IOException("Badly formatted HMM line: " + line);
      }

      if(transition_idx.size() != 7) {
        throw IOException("Badly formatted HMM line: " + line);      
      }

      // skip one line 
      std::getline(in, line);
      hmm_found = true;
      break;
    }
  }
  if (!hmm_found) {
    throw IOException("No HMM found in given input");
  }

  // set null model
  chunks = ost::StringRef(null_line.c_str(), null_line.length()).split();
  prof.SetNullModel(GetColumn(chunks, 1, olcs, "Badly formatted NULL line\n"
                              + null_line));

  // set columns
  while (std::getline(in, line)) {
    sline = ost::StringRef(line.c_str(), line.length());
    if (sline.trim().empty()) continue;
    if (sline.trim() == ost::StringRef("//", 2)) break;
    chunks = sline.split();
    // one letter code
    const char olc = chunks[0][0];
    // frequencies
    prof.AddColumn(GetColumn(chunks, 2, olcs, "Badly formatted line\n"
                             + line), olc);
    // get transition probabilities
    std::getline(in, line);
    sline = ost::StringRef(line.c_str(), line.length());
    if (sline.trim().empty()) continue;
    if (sline.trim() == ost::StringRef("//", 2)) break;
    chunks = sline.split();
    prof.back().SetHMMData(GetHMMData(chunks, transition_idx, 
                           transitions, neff_idx, neff_i_idx, neff_d_idx,
                           "Badly formatted line\n" + line));
  }

  // parse neff if it's there, calculate the average of every column otherwise
  if(neff_found) {
    sline = ost::StringRef(neff_line.c_str(), neff_line.length());
    chunks = sline.split();
    if(chunks.size() != 2) {
      throw IOException("Badly formatted line\n" + neff_line);
    }
    std::pair<bool, float> p = chunks[1].to_float();
    if (!p.first) {
      throw IOException("Badly formatted line\n" + neff_line);
    }
    prof.SetNeff(p.second);
  } else {
    Real summed_neff = 0.0;
    for(size_t i = 0; i < prof.size(); ++i) {
      summed_neff += prof[i].GetHMMData()->GetNeff();
    }
    if(prof.size() > 0) {
      prof.SetNeff(summed_neff/prof.size());
    } else{
      prof.SetNeff(1.0);
    }
  }
}

}}
