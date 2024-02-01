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

#include <map>
#include <ost/message.hh>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <ost/seq/sequence_handle.hh>

#include "clustal_io_handler.hh"
#include <ost/io/io_exception.hh>

/*
  Author: Marco Biasini
 */

using boost::format;

namespace ost { namespace io {

void ClustalIOHandler::Import(seq::SequenceList& aln,
                              const boost::filesystem::path& loc) {                        
    boost::filesystem::ifstream infile(loc);    
    this->Import(aln, infile);
}

void ClustalIOHandler::Export(const seq::ConstSequenceList& msa,
                            const boost::filesystem::path& loc) const {
}


bool ClustalIOHandler::ProvidesImport(const boost::filesystem::path& loc, 
                                    const String& format) {
  if (format=="auto") {
   String match_suf_string=loc.string();
   std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
   if (detail::FilenameEndsWith(match_suf_string,".aln")) {
     return true;
   }
  } else if(format=="clustal") {
    return true;
  }
  return false;                
}

bool ClustalIOHandler::ProvidesExport(const boost::filesystem::path& loc, 
                                      const String& format) {
  return false;
}

void ClustalIOHandler::Import(seq::SequenceList& aln,
                              std::istream& instream)
{
  static const char* whitespace="\t\r\n ";  
  int line_num=0;
  String line;
  typedef std::map<String, seq::SequenceHandle> SeqMap;
  std::vector<seq::SequenceHandle> order;
  SeqMap seq_map;
  if (!instream) {
    throw IOException("Can't import CLUSTAL alignment. Inexisting file "
                      "or invalid stream.");
  }
  
  // skip empty lines
  while (std::getline(instream, line)) {
    line_num+=1;
    size_t pos=line.find_first_not_of(whitespace);
    if (pos!=String::npos) {
      break;
    }    
  }
  if (line.find("CLUSTAL")!=0) {
    throw IOException("bad CLUSTAL file. First line must contain CLUSTAL");
  }
  while (std::getline(instream, line)) {
    line_num+=1;
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
      std::stringstream ss;
      ss << "Bad CLUSTAL file on line " << line_num 
         << ": Whitespace expected after sequence name";
      throw IOException(ss.str());
    }
    String seq_id=line.substr(0, seq_id_end);

    size_t start=line.find_first_not_of(whitespace, seq_id_end);
    size_t end=line.find_first_of(whitespace, start);
    end=end==String::npos ? line.length() : end;
    String seq_string=line.substr(start, end-start);      
    SeqMap::iterator entry=seq_map.find(seq_id);
    if (entry!=seq_map.end()) {
      entry->second.SetString(entry->second.GetString()+seq_string);
    } else {
      seq::SequenceHandle sequence=seq::CreateSequence(seq_id, seq_string);
      seq_map.insert(std::make_pair(seq_id, sequence));
      order.push_back(sequence);
    }
  }    
  for (std::vector<seq::SequenceHandle>::iterator i=order.begin(), 
       e=order.end(); i!=e; ++i) {
    aln.AddSequence(*i);
  }
  if (seq_map.size()==0) {
    throw IOException("Bad CLUSTAL file: File is empty");
  }
}

void ClustalIOHandler::Export(const seq::ConstSequenceList& ent,
                              std::ostream& ostream) const
{
  throw IOException("Export of ClustalW alignments to stream is not supported");
}

}}
