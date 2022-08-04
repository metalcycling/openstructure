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
#include "fasta_io_handler.hh"
#include <ost/message.hh>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <ost/io/io_exception.hh>
#include <ost/seq/invalid_sequence.hh>
/*
  Author: Marco Biasini
 */

using boost::format;

namespace ost { namespace io {

void FastaIOHandler::Import(seq::SequenceList& aln,
                            const boost::filesystem::path& loc) 
{
    if (!boost::filesystem::exists(loc)) {
      throw IOException("File does not exist");
    }
    
    boost::filesystem::ifstream infile(loc);
    this->Import(aln, infile);
}

void FastaIOHandler::Export(const seq::ConstSequenceList& msa,
                            const boost::filesystem::path& loc) const {
  boost::filesystem::ofstream outfile(loc);
  this->Export(msa, outfile);
}




bool FastaIOHandler::ProvidesImport(const boost::filesystem::path& loc, 
                                    const String& format) {
  if (format=="auto") {
   String match_suf_string=loc.string();
   std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
   if (detail::FilenameEndsWith(match_suf_string,".fasta") || detail::FilenameEndsWith(match_suf_string,".fa") ||
       detail::FilenameEndsWith(match_suf_string,".fna") || detail::FilenameEndsWith(match_suf_string,".fsa") ||
       detail::FilenameEndsWith(match_suf_string,".fas") )  {
	return true;
   }
  } else if(format=="fasta") {
    return true;
  }
  return false;                
}

bool FastaIOHandler::ProvidesExport(const boost::filesystem::path& loc, 
                           const String& format) {
  return FastaIOHandler::ProvidesImport(loc, format);
}

void FastaIOHandler::Import(seq::SequenceList& aln,
                            std::istream& instream)
{
  const char* error_msg="Bad FASTA file: Expected '>', but '%1%' found.";  
  
  String line;
  std::getline(instream, line);    
  int seq_count=0;
  while (!instream.eof()) {
    // parse header information. cut after first "|"
    if (line.find_first_not_of("\n\t ")==String::npos) {
      std::getline(instream, line);
      continue;
    }
    if (line.length()==0 || line[0]!='>') {
      String error=str(format(error_msg) % line);
      throw IOException(error);
    }
    String name=line.substr(1);
    std::stringstream seq_string;
    while (std::getline(instream, line)) {
      if (!line.empty() && line[0]=='>') {
        break;
      }
      for (String::iterator i=line.begin(), e=line.end(); i!=e; ++i) {
        if (isspace(*i)) {
          continue;
        }
        seq_string << *i;
      }
    }    
    if (!seq_string.str().empty()) {
      try {
        seq::SequenceHandle seq=seq::CreateSequence(name, seq_string.str());
        aln.AddSequence(seq);          
        seq_count+=1;
      } catch (seq::InvalidSequence& e) {
        throw seq::InvalidSequence("Failed for sequence with name " + name + ": " + e.what());
      }
    } else {
      throw IOException("Bad FASTA file: Sequence is empty.");
    }
  }    
  if (seq_count==0) {
    throw IOException("Bad FASTA file: File is empty");
  }                        

}

void FastaIOHandler::Export(const seq::ConstSequenceList& seqs,
                            std::ostream& ostream) const
{
  for (int i=0; i<seqs.GetCount(); ++i) {
    ostream << ">" << seqs[i].GetName() << std::endl;
    ostream << seqs[i].GetString() << std::endl;
  }
}

}}
