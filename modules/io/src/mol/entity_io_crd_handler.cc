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

#include <ost/log.hh>
#include <ost/conop/conop.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/profile.hh>
#include <ost/boost_filesystem_helper.hh>

#include <ost/io/io_exception.hh>
#include <ost/io/mol/io_profile.hh>
#include <ost/io/swap_util.hh>

#include "entity_io_crd_handler.hh"

namespace ost { namespace io {

using boost::format;


/// \brief Reader for CHARMM crd file format
///
/// Standard and extended CHARMM format is supported, and file format is
/// automatically detected based on file header or atom number
CRDReader::CRDReader(const boost::filesystem::path& loc):
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
std::vector<mol::AtomHandle> CRDReader::GetSequentialAtoms() const
{
  return sequential_atom_list_;
}


/// \brief Performes file import
void CRDReader::Import(mol::EntityHandle& ent)
{
  Profile profile_import("CRDReader::Import");

  String line;
  // skip over comments at the start
  while(std::getline(in_,line)) {
    if(line[0]!='*') break;
  }

  // line contains atom count
  std::vector<String> line_content;
  boost::trim(line);
  boost::split(line_content,line,boost::is_any_of(" "));

  // expanded charmm CARD format check
  if (line_content.size() > 1 || boost::lexical_cast<int>(line_content[0]) > 99999)
    while(std::getline(in_,line)) {
      if (!boost::trim_copy(line).empty())
        ParseAndAddAtomExpanded(line,ent);
    }
  else
    while(std::getline(in_,line)) {
      if (!boost::trim_copy(line).empty())
        ParseAndAddAtom(line,ent);
    }

  LOG_INFO("imported " << chain_count_ << " chains, " << residue_count_
                << " residues, " << atom_count_ << " atoms");  
}

/// \brief Parsing for standard format
void CRDReader::ParseAndAddAtom(const String& line, mol::EntityHandle& ent)
{
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);

  LOG_TRACE( "line: [" << line << "]" );

  //int anum = boost::lexical_cast<int>(boost::trim_copy(line.substr(0,5)));
  String aname = boost::trim_copy(line.substr(16,4));
  String ele = aname.substr(0,1);
  String rname = boost::trim_copy(line.substr(11,4));
  int irnum = boost::lexical_cast<int>(boost::trim_copy(line.substr(55,4)));
  String s_chain = boost::trim_copy(line.substr(51,4));
  geom::Vec3 apos(boost::lexical_cast<Real>(boost::trim_copy(line.substr(21,9))),
                  boost::lexical_cast<Real>(boost::trim_copy(line.substr(31,9))),
                  boost::lexical_cast<Real>(boost::trim_copy(line.substr(41,9))));

  mol::ResidueKey rkey(rname);
  
  // some postprocessing
  LOG_TRACE( "s_chain: [" << s_chain << "]" );

  mol::ResNum rnum(irnum);
  
  // determine chain and residue update
  bool update_chain=false;
  bool update_residue=false;
  if(!curr_chain_) {
    update_chain=true;
    update_residue=true;
  } else if(curr_chain_.GetName()!=s_chain) {
    update_chain=true;
    update_residue=true;
  }

  if(!curr_residue_) {
    update_residue=true;
  } else if(curr_residue_.GetNumber()!=rnum) {
    update_residue=true;
  }

  if(update_chain) {  
    if (!(curr_chain_=ent.FindChain(s_chain))) {
      LOG_DEBUG("new chain " << s_chain);      
      curr_chain_=editor.InsertChain(s_chain);
      ++chain_count_;      
    }
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
  sequential_atom_list_.push_back(ah);
  ++atom_count_;
}

/// \brief Parsing for extended format
void CRDReader::ParseAndAddAtomExpanded(const String& line, mol::EntityHandle& ent)
{
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  
  LOG_TRACE( "line: [" << line << "]" );
    
  //int anum = boost::lexical_cast<int>(boost::trim_copy(line.substr(0,5)));
  String aname = boost::trim_copy(line.substr(32,8));
  String ele = aname.substr(0,1);
  String rname = boost::trim_copy(line.substr(22,8));
  int irnum = boost::lexical_cast<int>(boost::trim_copy(line.substr(112,8)));
  String s_chain = boost::trim_copy(line.substr(102,8));
  geom::Vec3 apos(boost::lexical_cast<Real>(boost::trim_copy(line.substr(40,20))),
                  boost::lexical_cast<Real>(boost::trim_copy(line.substr(60,20))),
                  boost::lexical_cast<Real>(boost::trim_copy(line.substr(80,20))));
  Real b_factor = boost::lexical_cast<Real>(boost::trim_copy(line.substr(128,12)));
    
  mol::ResidueKey rkey(rname);
    
  // some postprocessing
  LOG_TRACE( "s_chain: [" << s_chain << "]" );
    
  mol::ResNum rnum(irnum);
  
  // determine chain and residue update
  bool update_chain=false;
  bool update_residue=false;
  if(!curr_chain_) {
    update_chain=true;
    update_residue=true;
  } else if(curr_chain_.GetName()!=s_chain) {
    update_chain=true;
    update_residue=true;
  }
    
  if(!curr_residue_) {
    update_residue=true;
  } else if(curr_residue_.GetNumber()!=rnum) {
    update_residue=true;
  }
    
  if(update_chain) {  
    if (!(curr_chain_=ent.FindChain(s_chain))) {
      LOG_DEBUG("new chain " << s_chain);      
      curr_chain_=editor.InsertChain(s_chain);
      ++chain_count_;      
    }
  }
    
  if(update_residue) {
    LOG_DEBUG("new residue " << rkey << " " << rnum);
    curr_residue_=editor.AppendResidue(curr_chain_, rkey, rnum);
    assert(curr_residue_.IsValid());
    ++residue_count_;
  }
    
  // finally add atom
  LOG_DEBUG("adding atom " << aname << " (" << ele << ") @" << apos);
  mol::AtomHandle ah = editor.InsertAtom(curr_residue_, aname, apos, ele,
                                         1.0, b_factor);
  sequential_atom_list_.push_back(ah);
  ++atom_count_;
}

/// \brief CHARMM format writer
///
/// The charmm format writer supports both the standard and the extended
/// format. Standard format supports at maximum 99999 atoms. If less atoms than
/// that are present, the default format is standard format, otherwise the
/// extended format is used. Extended format can be requested by setting the
/// parameter ext to true.
CRDWriter::CRDWriter(std::ostream& ostream, bool ext) :
     outfile_(), outstream_(ostream), ext_(ext), atom_count_(0), res_count_(0)
{}

CRDWriter::CRDWriter(const boost::filesystem::path& filename, bool ext) :
  outfile_(BFPathToString(filename).c_str()),
  outstream_(outfile_), ext_(ext), atom_count_(0), res_count_(0)
{}

CRDWriter::CRDWriter(const String& filename, bool ext) :
  outfile_(filename.c_str()), outstream_(outfile_), ext_(ext), atom_count_(0),
  res_count_(0)
{}

/// \brief Write header containing standard title and atom count/format row
void CRDWriter::WriteHeader(const mol::EntityView& ent)
{
  outstream_  << "* COOR FILE CREATED BY OPENSTRUCTURE" << std::endl;
  outstream_  << "*" << std::endl;
  
  atom_total_ = ent.GetAtomCount();
	
  if (atom_total_ > 99999 || ext_) {
    outstream_  << format("%10i") % ent.GetAtomCount() << "  EXT" << std::endl;
  }
  else {
    outstream_  << format("%5i") % ent.GetAtomCount() << std::endl;
  }
}

/// \brief Performs re-initialisation for writer
///
/// Currently, this only sets the residue count to zero
void CRDWriter::Init()
{
  res_count_ = 0;
}

void CRDWriter::Write(const mol::EntityView& ent)
{
  if (!outstream_) {
    throw IOException("Can't write CRD file. Bad output stream");
  }
  this->Init();
  mol::EntityView non_const_view = ent;
  this->WriteHeader(non_const_view);
  non_const_view.Apply(*this);
}

void CRDWriter::Write(const mol::EntityHandle& ent)
{
  if (!outstream_) {
    throw IOException("Can't write CRD file. Bad output stream");
  }
  this->Init();
  mol::EntityView non_const_view = ent.CreateFullView();
  this->WriteHeader(non_const_view);
  non_const_view.Apply(*this);
}


bool CRDWriter::VisitResidue(const mol::ResidueHandle& r)
{
  res_count_++;
  return true;
}

bool CRDWriter::VisitAtom(const mol::AtomHandle& atom)
{
  atom_count_++;
  String e_name=atom.GetEntity().GetName();
  if (e_name=="") {
    e_name="MOL";
  }
	
  mol::ResidueHandle res=atom.GetResidue();

  if (atom_total_ > 99999 || ext_) {
    outstream_  << format("%10i") % atom_count_
                << format("%10i") % res_count_ << "  "
          	    << format("%-8s") % res.GetKey() << "  "
                << format("%-8s") % atom.GetName()
                << format("%20.10f") % atom.GetPos().x
                << format("%20.10f") % atom.GetPos().y
                << format("%20.10f") % atom.GetPos().z << "  "
                << format("%-8s") % res.GetChain().GetName() << "  "
                << format("%-8i") % res.GetNumber()
                << format("%20.10f") % atom.GetBFactor()
                << std::endl;
  } else {
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
  }
  return true;
}

/// \brief CHARMM file format requires processor
bool EntityIOCRDHandler::RequiresProcessor() const
{
  return true;
}

void EntityIOCRDHandler::Import(mol::EntityHandle& ent,
                                const boost::filesystem::path& loc)
{
  CRDReader reader(loc);
  reader.Import(ent);
}

void EntityIOCRDHandler::Export(const mol::EntityView& ent,
                                std::ostream& stream) const
{
  CRDWriter writer(stream);
  writer.Write(ent);
}

void EntityIOCRDHandler::Export(const mol::EntityView& ent,
                                const boost::filesystem::path& loc) const
{
  CRDWriter writer(loc);
  writer.Write(ent);
}

namespace {

bool crd_handler_is_responsible_for(const boost::filesystem::path& loc,
                                    const String& type) {
  if(type=="auto") {
    String match_suf_string=loc.string();
	std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
    if( detail::FilenameEndsWith(match_suf_string,".crd") || detail::FilenameEndsWith(match_suf_string,".crd.gz") ) {
      return true;
    }
  } else if(type=="crd") {
    return true;
  }

  return false;
}

}

bool EntityIOCRDHandler::ProvidesImport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return crd_handler_is_responsible_for(loc, type);
}

bool EntityIOCRDHandler::ProvidesExport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return crd_handler_is_responsible_for(loc, type);
}

mol::EntityHandle LoadCRD(const String& file_name) 
{
  Profile profile_load("LoadCRD");
  CRDReader reader(file_name);
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  reader.Import(ent);
  IOProfile& prof = IOProfileRegistry::Instance().GetDefault();
  if (prof.processor) {
    prof.processor->Process(ent);
  }
  return ent;
}


void EntityIOCRDHandler::Import(mol::EntityHandle& ent, 
                                std::istream& stream)
{
  throw IOException("CRD format does not support import from stream");
}




}} // ns


