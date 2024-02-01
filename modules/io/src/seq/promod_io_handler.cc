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
#include <ctype.h>
#include <limits>

#include "promod_io_handler.hh"
#include <ost/io/io_exception.hh>

#include <ost/message.hh>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/erase.hpp>
/*
  Author: Marco Biasini
 */

using boost::format;

namespace ost { namespace io {

void PromodIOHandler::Import(seq::SequenceList& ali,
                              const boost::filesystem::path& loc) {                        
  boost::filesystem::ifstream infile(loc);
  this->Import(ali, infile);
}

void PromodIOHandler::Export(const seq::ConstSequenceList& msa,
                             const boost::filesystem::path& loc) const 
{
  throw IOException("Export of alignments to PROMOD Format is not supported");
}


bool PromodIOHandler::ProvidesImport(const boost::filesystem::path& loc, 
                                    const String& format) {
  if (format=="auto") {
   String match_suf_string=loc.string();
   std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
   if (detail::FilenameEndsWith(match_suf_string,".ali")) {
     return true;
   }
  } else if(format=="promod") {
    return true;
  }
  return false;                
}

bool PromodIOHandler::ProvidesExport(const boost::filesystem::path& loc, 
                                      const String& format) {
  return false;
}

void PromodIOHandler::Import(seq::SequenceList& ali,
                             std::istream& instream)
{
  static const char* whitespace="\t ";
  typedef std::map<String, seq::SequenceHandle> SeqMap;
  SeqMap seq_map;  
  String line;  
  std::vector<int> starts;
  while (std::getline(instream, line)) {
    size_t pos=line.find_first_not_of(whitespace);
    if (pos==String::npos) {
      continue;
    }
    if (pos>0) {
      // we have no use for sequence conservation lines at the moment
      continue;
    }
    size_t seq_id_end=line.find_first_of(whitespace);
    if (seq_id_end==String::npos) {
      throw IOException("Bad PROMOD file: Whitespace expected after"
                        " sequence name");
    }
    String seq_id=line.substr(0, seq_id_end);

    size_t num_start=line.find_first_not_of(whitespace, seq_id_end);
    size_t num_end=line.find_first_of(whitespace, num_start);
    size_t seq_start=line.find_first_not_of(whitespace, num_end);
    size_t seq_end=line.length();
    if (seq_start==String::npos || seq_end==String::npos ||
        !isdigit(line[num_start])) {
      continue;
    }
    String seq_string=line.substr(seq_start, seq_end-seq_start);
    boost::erase_all(seq_string, " ");
    boost::to_upper(seq_string);
    SeqMap::iterator entry=seq_map.find(seq_id);
    if (entry!=seq_map.end()) {
      entry->second.SetString(entry->second.GetString()+seq_string);
    } else {
      seq::SequenceHandle sequence=seq::CreateSequence(seq_id, seq_string);
      seq_map.insert(std::make_pair(seq_id, sequence));
      starts.push_back(seq_start);
    }
  }
  int min_start=std::numeric_limits<int>::max();
  for (std::vector<int>::iterator i=starts.begin(), e=starts.end();
       i!=e; ++i) {
    min_start=std::min(min_start, *i);
  }
  
  if (seq_map.size()==0) {
    throw IOException("Bad PROMOD file: File is empty");
  }
    
  for (SeqMap::iterator i=seq_map.begin(), e=seq_map.end(); i!=e; ++i) {
    int start=starts[std::distance(i, seq_map.begin())];
    if (min_start<start) {
      i->second.SetString(String(start-min_start, '-')+i->second.GetString());
    }    
    ali.AddSequence(i->second);
  }

}

void PromodIOHandler::Export(const seq::ConstSequenceList& ali,
                             std::ostream& ostream) const
{
  throw IOException("Export of alignments to PROMOD Format is not supported");
}

}}
