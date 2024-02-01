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
#include <vector>
#include <string>

#include <ctype.h>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/shared_array.hpp>

#include <ost/log.hh>
#include <ost/conop/conop.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/profile.hh>
#include <ost/io/mol/io_profile.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/swap_util.hh>

#include "entity_io_mae_handler.hh"

namespace ost { namespace io {

using boost::format;

MAEReader::MAEReader(const boost::filesystem::path& loc):
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

void MAEReader::Import(mol::EntityHandle& ent)
{
  Profile profile_import("MAEReader::Import");

  static const boost::regex r_ct_block("f_m_ct +\\{");
  static const boost::regex r_atom_block("m_atom\\[.*\\] +\\{");
  static const boost::regex r_delim(":::");
  static const boost::regex r_end(" *\\} *");

  String line;
  bool in_ct_block=false;
  bool in_atom_block=false;
  bool parsing_atoms=false;
  std::vector<std::string> prop_list;
  int i_atom_name=-1;
  int i_atom_xpos=-1;
  int i_atom_ypos=-1;
  int i_atom_zpos=-1;
  int i_res_name=-1;
  int i_res_num=-1;
  int i_chain_name=-1;

  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);

  boost::shared_array<char> xline(new char[512]);
  
  while(std::getline(in_,line)) {
    //line = boost::trim_copy(line);
    if(in_ct_block) {
      if(in_atom_block) {
        if(parsing_atoms) {
          if(boost::regex_search(line,r_delim)) {
            LOG_TRACE( "stopping atom parsing" );
            parsing_atoms=false;
            prop_list.clear();
            i_atom_name=-1;
            i_atom_xpos=-1;
            i_atom_ypos=-1;
            i_atom_zpos=-1;
            i_res_name=-1;
            i_res_num=-1;
            i_chain_name=-1;
          } else {
            // parsing atom line
            memset(xline.get(),0,512);
            strncpy(xline.get(),line.c_str(),511);
            parse_and_add_atom(ent,editor,xline.get(),line.size(),
                               i_atom_name, i_atom_xpos, i_atom_ypos, i_atom_zpos, 
                               i_res_name, i_res_num, i_chain_name);
          }
        } else { // not parsing atoms
          if(boost::regex_search(line,r_delim)) {
            if(i_atom_name==-1 ||
               i_atom_xpos==-1 ||
               i_atom_ypos==-1 ||
               i_atom_zpos==-1 ||
               i_res_name==-1 ||
               i_res_num==-1) {
              throw IOException("missing atom prop");
            }
            LOG_TRACE( "starting atom parsing" );
            parsing_atoms=true;
          } else if(boost::regex_search(line,r_end)) {
            LOG_TRACE( "exiting atom block" );
            in_atom_block=false;
          } else {
            std::string line2=boost::trim_copy(line);
            // parsing property line
            if(line2[0]!='#') {
              int pid=prop_list.size()+1;
              prop_list.push_back(line2);
              LOG_TRACE( "found property '" << prop_list.back() << "' id=" << pid );
              if(line2=="s_m_pdb_atom_name") i_atom_name=pid;
              else if(line2=="r_m_x_coord") i_atom_xpos=pid;
              else if(line2=="r_m_y_coord") i_atom_ypos=pid;
              else if(line2=="r_m_z_coord") i_atom_zpos=pid;
              else if(line2=="s_m_pdb_residue_name") i_res_name=pid;
              else if(line2=="i_m_residue_number") i_res_num=pid;
              else if(line2=="s_m_chain_name") i_chain_name=pid;
              else if(line2=="s_m_pdb_segment_name") {
                // only use this one if s_m_chain_name is not present
                if(i_chain_name<0) {
                  i_chain_name=pid;
                }
              }
            }
          }
        }
      } else { // not in atom block
        if(boost::regex_search(line,r_atom_block)) {
          LOG_TRACE( "entering atom block" );
          in_atom_block=true;
        } else if(boost::regex_search(line,r_end)) {
          LOG_TRACE( "exiting ct block" );
          in_ct_block=false;
        }
      }
    } else { // not in ct block
      if(boost::regex_search(line,r_ct_block)) {
        LOG_TRACE( "entering ct block" );
        in_ct_block=true;
      }
    }
  }

  LOG_INFO("imported " << chain_count_ << " chains, " << residue_count_
                << " residues, " << atom_count_ << " atoms");  
}


void MAEReader::parse_and_add_atom(mol::EntityHandle ent,
                                   mol::XCSEditor& editor,
                                   char* line,
                                   size_t line_len,
                                   int i_atom_name,
                                   int i_atom_xpos,
                                   int i_atom_ypos,
                                   int i_atom_zpos,
                                   int i_res_name,
                                   int i_res_num,
                                   int i_chain_name)  
{
  std::string aname("X");
  std::string rname("UNK");
  std::string cname("X");
  geom::Vec3 apos;
  int irnum=0;

  LOG_TRACE("[" << std::string(line) << "]");

  bool in_string=false;
  char* p0=line;
  char* p2=line+line_len;
  while(isspace(*p0) && p0<=p2) ++p0;
  int tcount=0;
  for(char* p1=p0+1;p1<p2;++p1) {
    if(in_string) {
      if((*p1)=='"') {
        in_string=false;
        char* p3=p1;
        (*p3)='\0';
        while(isspace(*p0) && p0<=p3) ++p0;
        while(isspace(*(p3-1)) && p3>p0) --p3;
        (*p3)='\0';
        if(tcount==i_atom_name) {
          aname=std::string(p0);
        } else if(tcount==i_atom_xpos) {
          apos[0]=atof(p0);
        } else if(tcount==i_atom_ypos) {
          apos[1]=atof(p0);
        } else if(tcount==i_atom_zpos) {
          apos[2]=atof(p0);
        } else if(tcount==i_res_name) {
          rname=std::string(p0);
        } else if(tcount==i_res_num) {
          irnum=atoi(p0);
        } else if(tcount==i_chain_name) {
          cname=std::string(p0);
        }
        ++tcount;
        ++p1;
        p0=p1;
        while(isspace((*p1)) && p1<p2) {++p0;++p1;}
        --p0; --p1;
      }
    } else { // not instring
      if(isspace((*p1))) {
        (*p1)='\0';
        ++p0;
        if(tcount==i_atom_name) {
          aname=std::string(p0);
        } else if(tcount==i_atom_xpos) {
          apos[0]=atof(p0);
        } else if(tcount==i_atom_ypos) {
          apos[1]=atof(p0);
        } else if(tcount==i_atom_zpos) {
          apos[2]=atof(p0);
        } else if(tcount==i_res_name) {
          rname=std::string(p0);
        } else if(tcount==i_res_num) {
          irnum=atoi(p0);
        } else if(tcount==i_chain_name) {
          cname=std::string(p0);
        }
        ++tcount;
        ++p1;
        p0=p1;
        while(isspace((*p1)) && p1<p2) {++p0;++p1;}
        --p0; --p1;
      } else if((*p1)=='"') {
        in_string=true;
        p0=p1+1;
      }
    }
  }

  std::string ele=aname.empty() ? "X" : aname.substr(0,1);
  if(isdigit(ele[0])) ele="H";

  mol::ResidueKey rkey(rname);
  
  mol::ResNum rnum(irnum);
  
  // determine chain and residue update
  bool update_chain=false;
  bool update_residue=false;
  if(!curr_chain_) {
    update_chain=true;
    update_residue=true;
  } else if(curr_chain_.GetName()!=cname) {
    update_chain=true;
    update_residue=true;
  }

  if(!curr_residue_) {
    update_residue=true;
  } else if(curr_residue_.GetNumber()!=rnum) {
    update_residue=true;
  }

  if(update_chain) {  
    if (!(curr_chain_=ent.FindChain(cname))) {
      curr_chain_=editor.InsertChain(cname);
      LOG_TRACE("new chain " << curr_chain_);      
      ++chain_count_;      
    } else {
      LOG_TRACE("old chain " << curr_chain_);      
    }
  }

  if(update_residue) {
    if (!(curr_residue_=curr_chain_.FindResidue(rnum))) {
      if(curr_chain_.GetResidueCount()>0) {
        int loc=curr_chain_.GetResidueCount()-1;
        for(;loc>=0;--loc) {
          if(curr_chain_.GetResidueByIndex(loc).GetNumber()<rnum) break;
        }
        if(loc<0) {
          curr_residue_=editor.InsertResidueBefore(curr_chain_,0,rnum,rkey);
        } else {
          curr_residue_=editor.InsertResidueAfter(curr_chain_,loc,rnum,rkey);
        }
        if(!curr_residue_) {
          // this should not happen...
          curr_residue_=editor.AppendResidue(curr_chain_, rkey, rnum);
        }
      } else {
        curr_residue_=editor.AppendResidue(curr_chain_, rkey, rnum);
      }
      assert(curr_residue_.IsValid());
      LOG_TRACE(" new residue " << curr_residue_);
      ++residue_count_;
    } else {
      LOG_TRACE(" old residue " << curr_residue_);
    }
  }

  // finally add atom
  LOG_TRACE("  atom " << aname << " (" << ele << ") @" << apos);
  ++atom_count_;
  mol::AtomHandle ah = editor.InsertAtom(curr_residue_, aname, apos, ele);
}
  



bool EntityIOMAEHandler::RequiresProcessor() const
{
  return true;
}

void EntityIOMAEHandler::Import(mol::EntityHandle& ent,
                                const boost::filesystem::path& loc)
{
  MAEReader reader(loc);
  reader.Import(ent);
}

void EntityIOMAEHandler::Export(const mol::EntityView& ent,
                                const boost::filesystem::path& loc) const
{
}

namespace {

bool mae_handler_is_responsible_for(const boost::filesystem::path& loc,
                                    const String& type) {
  if(type=="auto") {
    String match_suf_string=loc.string();
    std::transform(match_suf_string.begin(),match_suf_string.end(),match_suf_string.begin(),tolower);
    if( detail::FilenameEndsWith(match_suf_string,".mae") || detail::FilenameEndsWith(match_suf_string,".mae.gz") ) {
      return true;
    }
  } else if(type=="mae") {
    return true;
  }

  return false;
}

}

bool EntityIOMAEHandler::ProvidesImport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return mae_handler_is_responsible_for(loc, type);
}

bool EntityIOMAEHandler::ProvidesExport(const boost::filesystem::path& loc,
                                        const String& type)
{
  return false;
}

mol::EntityHandle LoadMAE(const String& file_name) 
{
  MAEReader reader(file_name);
  mol::EntityHandle ent=mol::CreateEntity();
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  reader.Import(ent);
  IOProfile& prof = IOProfileRegistry::Instance().GetDefault();
  if (prof.processor) {
    prof.processor->Process(ent);
  }
  return ent;
}


void EntityIOMAEHandler::Import(mol::EntityHandle& ent, 
                                std::istream& stream)
{
  throw IOException("MAE format does not support import from stream");
}


void EntityIOMAEHandler::Export(const mol::EntityView& ent, 
                                std::ostream& stream) const
{
  throw IOException("MAE format does not support export to stream");
}

}} // ns
