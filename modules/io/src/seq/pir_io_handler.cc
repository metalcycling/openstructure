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
#include "pir_io_handler.hh"
#include <ost/message.hh>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

#include <ost/io/io_exception.hh>
#include <ost/seq/invalid_sequence.hh>
/*
  Author: Juergen Haas
 */

using boost::format;

namespace ost { namespace io {

void PirIOHandler::Import(seq::SequenceList& aln,
                            const boost::filesystem::path& loc) 
{
    if (!boost::filesystem::exists(loc)) {
      throw IOException("File does not exist");
    }
    
    boost::filesystem::ifstream infile(loc);
    this->Import(aln, infile);
}

void PirIOHandler::Export(const seq::ConstSequenceList& msa,
                            const boost::filesystem::path& loc) const {
  boost::filesystem::ofstream outfile(loc);
  this->Export(msa, outfile);
}




bool PirIOHandler::ProvidesImport(const boost::filesystem::path& loc, 
                                    const String& format) {
  if (format=="auto") {
   String match_suf_string=loc.string();
   std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
   if (detail::FilenameEndsWith(match_suf_string,".pir"))  {
    return true;
   }
  } else if(format=="pir") {
    return true;
  }
  return false;
}

bool PirIOHandler::ProvidesExport(const boost::filesystem::path& loc, 
                           const String& format) {
  return PirIOHandler::ProvidesImport(loc, format);
}

void PirIOHandler::Import(seq::SequenceList& aln, std::istream& instream)
{
  const char* error_msg="Expected %1%, but '%2%' found.";  
  const char* premature_end="Premature end in PIR file on line %1%";
  const char* empty_seq="PIR file contains empty sequence on line %1%";
  String line;
  int line_num=1, seq_count=0;
  // we can't use a normal while(std::getline(...)) here, because the PIR format 
  // requires us to perform a one-line look-ahead. 
  std::getline(instream, line);  
  while (!instream.eof()) {
    // skip empty lines
    if (line.find_first_not_of("\n\t ")==String::npos) {
      continue;
    }
    if (line[0]!='>') {
      String error=str(format(error_msg) % "'>'" % line);
      throw IOException(error);
    }
    if (line.length()<4) {
      String error=str(format(error_msg) % "seq-type;identifier" % line);
      throw IOException(error);
    }
    String seq_type=line.substr(1, 2);
    if (!(seq_type=="P1" || seq_type=="F1" || seq_type=="DL" || 
          seq_type=="DC" || seq_type=="RL" || seq_type=="RC" || 
          seq_type=="XX")) {
      throw IOException("Bad PIR file: Unknown sequence type"/*" '%1%' on line %2%" % seq_type % line_nu*m*/);
    }
    String name=line.substr(4);
    int subject_line_num=line_num;
    if (!std::getline(instream, line)) {
      throw IOException(str(format(premature_end) % line_num));
    }
    line_num+=1;
    std::stringstream seq_string;
    bool end_seq=false;
    while (!end_seq && std::getline(instream, line)) {
      line_num+=1;
      for (String::iterator i=line.begin(), e=line.end(); i!=e; ++i) {
        if (isspace(*i)) {
          continue;
        }
        if (*i=='*') {
          end_seq=true;
          break;
        }
        seq_string << *i;
      }
    }
    while (std::getline(instream, line)) {
      line_num+=1;
      if (!line.empty() && line[0]=='>') {
        break;
      }
    }
    if (!end_seq) {
      throw IOException(str(format(premature_end) % line_num));
    }
    if (seq_string.str().empty()) {
      throw IOException(str(format(empty_seq) % subject_line_num));
    }
    try {
      seq::SequenceHandle ss=seq::CreateSequence(name, seq_string.str());
      aln.AddSequence(ss);          
      seq_count+=1;
    } catch (seq::InvalidSequence& e) {
      throw e;
    }
  }
  if (seq_count==0) {
    throw IOException("Bad Pir file: File is empty");
  }
}

void PirIOHandler::Export(const seq::ConstSequenceList& seqs,
                            std::ostream& ostream) const
{
  for (int i=0; i<seqs.GetCount(); ++i) {
    ostream << ">P1;" << seqs[i].GetName() << std::endl;
    ostream << "title" << std::endl;
    ostream << seqs[i].GetString() <<"*" << std::endl;
  }
}

}}
