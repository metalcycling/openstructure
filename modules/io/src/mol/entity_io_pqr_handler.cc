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

#include <iostream>
#include <sstream>
#include <iomanip>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/tokenizer.hpp>

#include <ost/log.hh>
#include <ost/conop/conop.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/profile.hh>
#include <ost/boost_filesystem_helper.hh>

#include <ost/io/io_exception.hh>
#include <ost/io/swap_util.hh>
#include <ost/io/mol/io_profile.hh>

#include "entity_io_pqr_handler.hh"

namespace ost { namespace io {

using boost::format;

/// \brief Reader for PQR file format, used in APBS
PQRReader::PQRReader(const boost::filesystem::path& loc):
  sequential_atom_list_(),
  curr_chain_(),
  curr_residue_(),
  chain_count_(0),
  residue_count_(0),
  atom_count_(0),
  infile_(loc),
  in_()
{
  if (boost::iequals(".gz", boost::filesystem::extension(loc))) {
    in_.push(boost::iostreams::gzip_decompressor());    
  }
  in_.push(infile_);  
  if(!infile_) throw IOException("could not open "+loc.string());
}

/// \brief Returns an vector containing all atom handles as observed in the file
std::vector<mol::AtomHandle> PQRReader::GetSequentialAtoms() const
{
  return sequential_atom_list_;
}

/// \brief Performes file import
void PQRReader::Import(mol::EntityHandle& ent)
{
  Profile profile_import("PQRReader::Import");

  String line;
  while(std::getline(in_,line)) {
    if(line.substr(0,4)=="ATOM") {
      ParseAndAddAtom(line,ent);
    }
  }

  LOG_INFO("imported " << chain_count_ << " chains, " << residue_count_
                << " residues, " << atom_count_ << " atoms");  
}

/// \brief Parsing for standard format
void PQRReader::ParseAndAddAtom(const String& line, mol::EntityHandle& ent)
{
  static bool auto_chain=false;
  static boost::char_separator<char> sep(" ");
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);

  //int anum = boost::lexical_cast<int>(boost::trim_copy(line.substr(0,5)));
  String aname = boost::trim_copy(line.substr(12,4));
  String ele = aname.substr(0,1);
  String rname = boost::trim_copy(line.substr(17,3));
  int irnum = boost::lexical_cast<int>(boost::trim_copy(line.substr(22,4)));
  String s_chain = boost::trim_copy(line.substr(21,1));

  // space separated x y z charge rad at column 30+
  Real ax=0.0,ay=0.0,az=0.0,ac=0.0,ar=0.0;
  try {
    std::string tmp(line.substr(30,std::string::npos));
    boost::tokenizer<boost::char_separator<char> > tok(tmp,sep);
    boost::tokenizer<boost::char_separator<char> >::iterator tit=tok.begin();
    ax=boost::lexical_cast<Real>(*tit); ++tit;
    ay=boost::lexical_cast<Real>(*tit); ++tit;
    az=boost::lexical_cast<Real>(*tit); ++tit;
    ac=boost::lexical_cast<Real>(*tit); ++tit;
    ar=boost::lexical_cast<Real>(*tit); ++tit;
  } catch (boost::bad_lexical_cast& e) {
    return;
  }
  
  geom::Vec3 apos(ax,ay,az);

  mol::ResidueKey rkey(rname);
  
  mol::ResNum rnum(irnum);
  
  // determine chain and residue update
  bool update_chain=false;
  bool update_residue=false;
  if(!curr_chain_) {
    update_chain=true;
    update_residue=true;
    if(s_chain.empty()) {
      s_chain="A";
      auto_chain=true;
    }
  } else {
    if(auto_chain) {
      if(curr_residue_) {
        if(rnum<curr_residue_.GetNumber()) {
          update_chain=true;
          update_residue=true;
          std::ostringstream tmp;
          tmp << static_cast<char>(65+chain_count_);
          s_chain=tmp.str();
        }
      }
    } else if(curr_chain_.GetName()!=s_chain) {
      update_chain=true;
      update_residue=true;
    }  
  }

  if(!curr_residue_) {
    update_residue=true;
  } else if(curr_residue_.GetNumber()!=rnum) {
    update_residue=true;
  }

  if(update_chain) {  
    LOG_DEBUG("new chain " << s_chain);      
    curr_chain_=editor.InsertChain(s_chain);
    ++chain_count_;      
  }

  if(update_residue) {
    LOG_DEBUG("new residue " << rkey << " " << rnum);
    curr_residue_=editor.AppendResidue(curr_chain_, rkey, rnum);
    assert(curr_residue_.IsValid());
    ++residue_count_;
  }

  // finally add atom
  LOG_DEBUG("adding atom " << aname << " (" << ele << ") @" << apos);
  mol::AtomHandle ah = editor.InsertAtom(curr_residue_, aname, apos, ele);
  ah.SetRadius(ar);
  ah.SetCharge(ac);
  sequential_atom_list_.push_back(ah);
  ++atom_count_;
}

PQRWriter::PQRWriter(std::ostream& ostream, bool ext) :
     outfile_(), outstream_(ostream), ext_(ext), atom_count_(0), res_count_(0)
{}

PQRWriter::PQRWriter(const boost::filesystem::path& filename, bool ext) :
  outfile_(BFPathToString(filename).c_str()),
  outstream_(outfile_), ext_(ext), atom_count_(0), res_count_(0)
{}

PQRWriter::PQRWriter(const String& filename, bool ext) :
  outfile_(filename.c_str()), outstream_(outfile_), ext_(ext), atom_count_(0),
  res_count_(0)
{}

/// \brief Write header containing standard title and atom count/format row
void PQRWriter::WriteHeader(const mol::EntityView& ent)
{
}

void PQRWriter::Init()
{
  res_count_ = 0;
}

void PQRWriter::Write(const mol::EntityView& ent)
{
  if (!outstream_) {
    throw IOException("Can't write PQR file. Bad output stream");
  }
  this->Init();
  mol::EntityView non_const_view = ent;
  this->WriteHeader(non_const_view);
  non_const_view.Apply(*this);
}

void PQRWriter::Write(const mol::EntityHandle& ent)
{
  if (!outstream_) {
    throw IOException("Can't write PQR file. Bad output stream");
  }
  this->Init();
  mol::EntityView non_const_view = ent.CreateFullView();
  this->WriteHeader(non_const_view);
  non_const_view.Apply(*this);
}


bool PQRWriter::VisitResidue(const mol::ResidueHandle& r)
{
  res_count_++;
  return true;
}

bool PQRWriter::VisitAtom(const mol::AtomHandle& atom)
{
  atom_count_++;
  String e_name=atom.GetEntity().GetName();
  if (e_name=="") {
    e_name="MOL";
  }
	
  mol::ResidueHandle res=atom.GetResidue();
/*
  outstream_  << format("%5i") % atom_count_
              << format("%5i") % res_count_ << " "
	            << format("%-4s") % res.GetKey() << " "
              << format("%-4s") % atom.GetName()
              << format("%10.5f") % atom.GetPos().x
              << format("%10.5f") % atom.GetPos().y
              << format("%10.5f") % atom.GetPos().z << " "
              << format("%-4s") % res.GetChain().GetName() << " "
              << format("%-5i") % res.GetNumber() << " "
              << format("%8.5f") % atom.GetBFactor()
              << std::endl;
              */
 
  return true;
}

bool EntityIOPQRHandler::RequiresBuilder() const
{
  return true;
}

void EntityIOPQRHandler::Import(mol::EntityHandle& ent,
                                const boost::filesystem::path& loc)
{
  PQRReader reader(loc);
  reader.Import(ent);
}

void EntityIOPQRHandler::Export(const mol::EntityView& ent,
                                std::ostream& stream) const
{
  PQRWriter writer(stream);
  writer.Write(ent);
}

void EntityIOPQRHandler::Export(const mol::EntityView& ent,
                                const boost::filesystem::path& loc) const
{
  PQRWriter writer(loc);
  writer.Write(ent);
}

namespace {

bool pqr_handler_is_responsible_for(const boost::filesystem::path& loc,
                                    const String& type) {
  if(type=="auto") {
    String match_suf_string=loc.string();
	std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
    if( detail::FilenameEndsWith(match_suf_string,".pqr") || detail::FilenameEndsWith(match_suf_string,".pqr.gz") ) {
      return true;
    }
  } else if(type=="pqr") {
    return true;
  }

  return false;
}

}

bool EntityIOPQRHandler::ProvidesImport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return pqr_handler_is_responsible_for(loc, type);
}

bool EntityIOPQRHandler::ProvidesExport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return pqr_handler_is_responsible_for(loc, type);
}

mol::EntityHandle LoadPQR(const String& file_name) 
{
  Profile profile_load("LoadPQR");
  IOProfile default_profile=IOProfileRegistry::Instance().GetDefault();
  PQRReader reader(file_name);
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  reader.Import(ent);
  if (default_profile.processor) 
    default_profile.processor->Process(ent);
  return ent;
}


void EntityIOPQRHandler::Import(mol::EntityHandle& ent, 
                                std::istream& stream)
{
  throw IOException("PQR format does not support import from stream");
}




}} // ns


