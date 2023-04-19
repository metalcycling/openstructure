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
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <ost/dyn_cast.hh>
#include <ost/profile.hh>
#include <ost/log.hh>
#include <ost/message.hh>
#include <ost/mol/bond_handle.hh>

#include <ost/conop/conop.hh>
#include <ost/geom/mat3.hh>
#include <ost/io/io_exception.hh>
#include "pdb_reader.hh"



namespace ost { namespace io {

using boost::format;

namespace {

  bool IEquals(const StringRef& a, const StringRef& b)
  {
    if (a.size()!=b.size()) {
      return false;
    }
    for (size_t i=0; i<a.size(); ++i) {
      if (toupper(a[i])!=b[i]) {
        return false;
      }
    }
    return true;
  }

  mol::ResNum to_res_num(int num, char ins_code)
  {
    return mol::ResNum(num, ins_code==' ' ? '\0' : ins_code);
  }

}

PDBReader::PDBReader(std::istream& instream, const IOProfile& profile):
  infile_(), instream_(instream), compnds_(), profile_(profile)
{
  this->Init(boost::filesystem::path(""));
}

PDBReader::PDBReader(const String& filename, const IOProfile& profile)
  : infile_(filename), instream_(infile_), compnds_(), profile_(profile)
{
  this->Init(boost::filesystem::path(filename));
}

PDBReader::PDBReader(const boost::filesystem::path& loc, 
                     const IOProfile& profile):
  infile_(loc), instream_(infile_), compnds_(), profile_(profile)
{
  this->Init(loc);
}

void PDBReader::Init(const boost::filesystem::path& loc)
{
  warned_name_mismatch_=false;
  read_seqres_=false;
  warned_rule_based_=false;
  charmm_style_=profile_.dialect=="CHARMM";
  num_model_records_=0;
  if (boost::iequals(".gz", boost::filesystem::extension(loc))) {
    in_.push(boost::iostreams::gzip_decompressor());
  }
  in_.push(instream_);
  if(!infile_) throw IOException("could not open "+loc.string());
  line_num_=0;
  if(boost::iequals(boost::filesystem::extension(loc), ".pqr")) {
    is_pqr_=true;
  } else {
    is_pqr_=false;
  }
  hard_end_=false;
  skip_next_=false;
  data_continues_=false;
  old_key_="";
  mol_id_=std::make_pair(false, 0);
  
}

void PDBReader::ThrowFaultTolerant(const String& msg) {
  if (profile_.fault_tolerant) {
    LOG_WARNING(msg);
    return;
  }
  throw IOException(msg);
}

void PDBReader::ParseCompndEntry (const StringRef& line, int line_num)
{
  if (line.size()<12) {
    if (profile_.fault_tolerant) {
      LOG_WARNING("invalid COMPND record on line " << line_num 
                  << ": record is too short");
      return;
    }
    std::stringstream ss;
    ss << "invalid COMPND record on line ";
    ss << line_num <<": record is too short";
    throw IOException(ss.str());
  }
  if (line.rtrim().size()>80) {
    if (profile_.fault_tolerant) {
      LOG_WARNING("invalid COMPND record on line " << line_num 
                  << ": record is too long");
      return;
    }
    std::stringstream ss;
    ss << "invalid COMPND record on line ";
    ss << line_num <<": whole record is too long";
    throw IOException(ss.str());
  }
  StringRef entry=line.substr(10,line.size()-10).trim();

  StringRef data, key;
  if (entry.size()<2){
    ThrowFaultTolerant(str(format("invalid COMPND record on line %d, compnd record too small")%line_num));
  }
  char last_letter=entry[entry.size()-1];
  if((last_letter==',') ||
      (last_letter==';')) {
      entry=entry.substr(0, entry.size()-1); //trim last char
  }
  if((entry.find(':')!=entry.end())){
    std::vector<StringRef> fields=entry.split(':');
      key=fields[0].trim();
      old_key_=key.str();
    if (fields.size()>1) {
      data=fields[1].trim();
    }


    if(data.size()<1){
      if (!(key.str()=="MOL_ID")&&!(key.str()=="CHAIN")){
        LOG_WARNING("skipping unsupported COMPND record on line " << line_num<< ": record value"<< key.str()<<" too small");
        if (data_continues_) {
          skip_next_=true;
        } else {
          return;
        }
      }
      ThrowFaultTolerant(str(format("invalid COMPND record on line %d, record after ':' too small")%line_num));
    }
    data_continues_=true;
    if (last_letter==';') {
      data_continues_=false;
    }
  } else if(skip_next_) {
    if (last_letter==';'){
      skip_next_=false;
      data_continues_=false;
      old_key_="";
    } else if (last_letter==','){
      data_continues_=true;
    }
    return;
  } else if (data_continues_){
    data=entry.trim();
    
    if(data.size()<1){
      ThrowFaultTolerant(str(format("invalid COMPND record on line %d, record after\
                                    ':' was empty")%line_num));
    }
    key=StringRef(old_key_.data(), old_key_.size());

  } 
      //currently only these are parsed
  if (!(key.str()=="MOL_ID")&&!(key.str()=="CHAIN")){
    LOG_INFO("reading COMPND record on line " << line_num<< "is not supported");
    if (data_continues_) {
      skip_next_=true;
    } else {
      return;
    }
  }
  
  std::vector<StringRef> chain_list;
  std::vector<String> chains;
  if ((IEquals(key, StringRef("MOL_ID", 6)))) {
    mol_id_=data.trim().to_int();
    if (mol_id_.first) {
      LOG_TRACE("COMPND record on line " << line_num<< " MOL_ID: "<<mol_id_.second);
    }
    if (!mol_id_.first) {
     ThrowFaultTolerant(str(format("invalid COMPND record on line %d")%line_num));
    }
  }
  if (IEquals(key, StringRef("CHAIN", 5))) {
    if (!mol_id_.first) {
    ThrowFaultTolerant(str(format("invalid COMPND record on line %d, CHAIN must be succeeding MOL_ID ")%line_num));
    }
    if (data.find(',')!=data.end()) {
      chain_list=data.split(',');
    } else {
      if(data.size()==1){
        chain_list.push_back(data);
      } else {
       ThrowFaultTolerant(str(format("invalid COMPND record on line %d, CHAIN must be succeeding MOL_ID ")%line_num));
      }
    }
    for (std::vector<StringRef>::const_iterator it = chain_list.begin(); it != chain_list.end(); ++it) {
      chains.push_back(it->trim().str());
    }
    compnds_.push_back(CompndEntry(chains, mol_id_.second));
  }
}


void PDBReader::ParseSeqRes(const StringRef& line, int line_num)
{
  conop::CompoundLibPtr comp_lib;
  comp_lib = conop::Conopology::Instance().GetDefaultLib();
   
  if (!comp_lib) {
    if (!warned_rule_based_) {
      LOG_WARNING("SEQRES import requires a compound library. Ignoring" 
                  "SEQRES records");      
    }
    warned_rule_based_=true;
    return;
  }
  if (!seqres_.IsValid()) {
    seqres_=seq::CreateSequenceList();
  }
  if (line.size()<17) {
    if (profile_.fault_tolerant) {
      LOG_WARNING("invalid SEQRES record on line " << line_num 
                  << ": record is too short");
      return;
    }
    std::stringstream ss;
    ss << "invalid SEQRES record on line ";
    ss << line_num <<": record is too short";
    throw IOException(ss.str());
  }
  String chain(1, line[11]);
  seq::SequenceHandle curr_seq;
  if (seqres_.GetCount()==0 || 
      seqres_[seqres_.GetCount()-1].GetName()!=chain) {
    curr_seq=seq::CreateSequence(chain, "");
    seqres_.AddSequence(curr_seq);

  } else {
    curr_seq=seqres_[seqres_.GetCount()-1];
  }
  for (int i=0; i<14; ++i) {
    size_t start=19+i*4;
    if (line.size()<start+3) {
      return;
    }
    StringRef rname=line.substr(start, 3);
    StringRef trimmed=rname.trim();
    if (trimmed.empty()) {
      return;
    }
    conop::CompoundPtr compound=comp_lib->FindCompound(trimmed.str(), 
                                                       conop::Compound::PDB);
    if (!compound) {
      if (rname!=StringRef("UNK", 3)) {
     
        LOG_WARNING("unknown residue '" << trimmed << "' in SEQRES record. "
                    "Setting one-letter-code to '?'");
      }
      curr_seq.Append('?');
      continue;
    }
    curr_seq.Append(compound->GetOneLetterCode());
  }
}

bool PDBReader::HasNext()
{
  if (hard_end_) {
    return false;
  }
  // since helix and sheet have to appear before any atom/hetatm records
  // a HELIX/SHEET entry implies a next model.
  while (std::getline(in_, curr_line_) && ++line_num_) {
     StringRef curr_line(curr_line_.c_str(), curr_line_.length());
     if (curr_line.size()>5 && 
         (IEquals(curr_line.substr(0, 6), StringRef("ATOM  ", 6)) ||
         (!profile_.no_hetatms &&
          IEquals(curr_line.substr(0, 6),StringRef("HETATM ", 6))) ||
          IEquals(curr_line.substr(0, 6),StringRef("ANISOU ", 6)) ||
         IEquals(curr_line.substr(0, 6), StringRef("SHEET ", 6)) ||
         IEquals(curr_line.substr(0, 6), StringRef("HELIX ", 6)) ||
         IEquals(curr_line.substr(0, 6), StringRef("COMPND", 6)) ||
         IEquals(curr_line.substr(0, 6), StringRef("MODEL ", 6)) ||
         IEquals(curr_line.substr(0, 6), StringRef("SEQRES", 6)) ||
         IEquals(curr_line.substr(0, 6), StringRef("HET   ", 6)))) {
       return true;
     } else if (IEquals(curr_line.rtrim(), StringRef("END", 3))) {
       hard_end_=true;
       return false;
     }
  }
  return false;
}

void PDBReader::Import(mol::EntityHandle& ent,
                       const String& restrict_chains)
{
  LOG_DEBUG("PDBReader: " << profile_);

  Profile profile_import("PDBReader::Import");
  this->ClearState();
  // first read curr_line and then read next...
  restrict_chains_=restrict_chains;
  bool go_on=true;
  do {
    if (curr_line_.empty()) {
      continue;
    }
    StringRef curr_line(curr_line_.c_str(), curr_line_.size());
    switch(curr_line[0]) {
      case 'A':
      case 'a':
        if (curr_line.size()<6) {
          continue;
        }
        if (IEquals(curr_line.substr(0, 6), StringRef("ATOM  ", 6))) {
          LOG_TRACE("processing ATOM entry");
          this->ParseAndAddAtom(curr_line, line_num_, ent, 
                                StringRef("ATOM", 4));
        } else if (IEquals(curr_line.substr(0, 6), StringRef("ANISOU", 6))) {
          if (!charmm_style_) {
            LOG_TRACE("processing ANISOU entry");
            this->ParseAnisou(curr_line, line_num_, ent);            
          }
        }
        break;
         case 'C':
         case 'c':
           if (curr_line.size()<6) {
             LOG_TRACE("skipping entry");
             continue;
           }
           else if (IEquals(curr_line.substr(0, 6), StringRef("CONECT", 6)) &&
                    profile_.read_conect) {
             LOG_TRACE("processing CONECT entry");
             this->ParseConectEntry(curr_line, line_num_, ent);
           }
           else if (IEquals(curr_line.substr(0, 6), StringRef("COMPND", 6))) {
             LOG_TRACE("processing COMPND entry");
             this->ParseCompndEntry(curr_line, line_num_);
           }
          break;
      case 'E':
      case 'e':
        if (curr_line.size()<3) {
          continue;
        }
        if (IEquals(curr_line.rtrim(), StringRef("END", 3))) {
          hard_end_=true;
          go_on=false;
          break;
        }
        if (curr_line.size()<6) {
          continue;
        }
        if (IEquals(curr_line.substr(0, 6), StringRef("ENDMDL", 6))) {
          go_on=false;
          num_model_records_=0;
          break;
        }
      case 'H':
      case 'h':
        if (curr_line.size()<6) {
          continue;
        }
        if (IEquals(curr_line.substr(0, 6), StringRef("HETATM", 6))) {
          if (profile_.no_hetatms)
            continue;
          LOG_TRACE("processing HETATM entry");
          this->ParseAndAddAtom(curr_line, line_num_, ent, 
                                StringRef("HETATM", 6));
        } else if (IEquals(curr_line.substr(0, 6), StringRef("HELIX ", 6))) {
          if (!charmm_style_) {
            this->ParseHelixEntry(curr_line);            
          }
        } else if (IEquals(curr_line.substr(0, 6), StringRef("HET   ", 6))) {
          // remember het entry to mark the residues as ligand during ATOM import
          char chain=curr_line[12];
          std::pair<bool, int> num=curr_line.substr(13, 4).ltrim().to_int();
          if (!num.first) {
            if (profile_.fault_tolerant) {
              LOG_WARNING("Invalid HET entry on line " << line_num_);
              continue;
            } else {
              String msg=str(format("Invalid HET entry on line %d")%line_num_);
              throw IOException(msg);
            }
          }
          hets_.push_back(HetEntry(chain, to_res_num(num.second,
                                   curr_line[17])));
        }
        break;
      case 'M':
      case 'm':
        if (curr_line.size()<6) {
          continue;
        }
        if (IEquals(curr_line.substr(0, 6), StringRef("MODEL ", 6))) {
          ++num_model_records_;          
          if (num_model_records_<2) {
            continue;
          }
          if (profile_.fault_tolerant) {
            go_on=false;
            num_model_records_=1;
            break;
          }
          String msg=str(format("MODEL record without matching ENDMDL on line %d")%line_num_);
          throw IOException(msg);
        }
        break;
      case 'S':
      case 's':
        if (curr_line.size()<6) {
          continue;
        }
        if (IEquals(curr_line.substr(0, 6), StringRef("SHEET ", 6))) {
          if (!charmm_style_) {
            this->ParseStrandEntry(curr_line);
          }
        }
        if (IEquals(curr_line.substr(0, 6), StringRef("SEQRES", 6))) {
          if (read_seqres_) {
            this->ParseSeqRes(curr_line, line_num_);
          }
        }
        break;
      default:
        break;
    }
  }  while (std::getline(in_, curr_line_) && ++line_num_ && go_on);
  LOG_INFO("imported "
               << chain_count_ << " chains, "
               << residue_count_ << " residues, "
               << atom_count_ << " atoms; with "
               << helix_list_.size() << " helices and "
               << strand_list_.size() << " strands");
  this->AssignSecStructure(ent);
  this->AssignMolIds(ent);
  for (HetList::const_iterator i=hets_.begin(), e=hets_.end(); i!=e; ++i) {
    mol::ResidueHandle res=ent.FindResidue(String(1, i->chain), i->num);
    if (res.IsValid()) {
      res.SetIsLigand(true);
    }      
  }
}

void PDBReader::AssignMolIds(mol::EntityHandle ent) {
  LOG_VERBOSE("Assigning MOL_IDs");
  for (CompndList::const_iterator compnd_iterator=compnds_.begin(), e=compnds_.end();
       compnd_iterator!=e; ++compnd_iterator) {
    for (std::vector<String>::const_iterator chain_iterator = compnd_iterator->chains.begin();
                                             chain_iterator!= compnd_iterator->chains.end();
                                             ++chain_iterator) {
      if (restrict_chains_.size()==0 ||
        (restrict_chains_.find(*chain_iterator)!=String::npos)) {
        mol::ChainHandle chain=ent.FindChain(*chain_iterator);
        if (chain) {
          chain.SetIntProp("mol_id", compnd_iterator->mol_id);
        }else if (!profile_.no_hetatms){
          // only throw if no_hetatms flag is False. Some chains might be
          // missing if no_hetatms is true and they're solely composed of
          // hetatms.
          LOG_WARNING("failed to assign MOL_ID to chain: "<<*chain_iterator <<std::endl);
          std::stringstream ss;
          ss << "could not map COMPND record MOL_ID onto chain";
          ss <<*chain_iterator;
          ThrowFaultTolerant(ss.str());
        }
      }
    }
  }
  if (! compnds_.empty()){
    mol::ChainHandleList ch_list=ent.GetChainList();
    for (mol::ChainHandleList::const_iterator chain=ch_list.begin();
         chain!=ch_list.end(); ++chain) {
           //~ skip HETATM only chains!!
      if(chain->IsValid()){
        if (!chain->HasProp("mol_id")) {
          std::stringstream ss;
          ss << "found chain without MOL_ID: ";
          ss << chain->GetName();
          LOG_WARNING(ss.str());
        }
      }
    }
  }
}

void PDBReader::AssignSecStructure(mol::EntityHandle ent)
{
  // and now for the fun part, assigning the strands/helices to residues
  for (HSList::const_iterator i=helix_list_.begin(), e=helix_list_.end();
       i!=e; ++i) {
    mol::ChainHandle chain=ent.FindChain(i->chain);
    if (!chain.IsValid()) {
      LOG_INFO("ignoring helix record for unknown chain "+i->chain);
      continue;
    }
    mol::SecStructure alpha(mol::SecStructure::ALPHA_HELIX);
    chain.AssignSecondaryStructure(alpha, i->start, i->end);
  }

  for (HSList::const_iterator i=strand_list_.begin(), e=strand_list_.end();
       i!=e; ++i) {
    mol::ChainHandle chain=ent.FindChain(i->chain);
    if (!chain.IsValid()) {
      LOG_INFO("ignoring strand record for unknown chain "+i->chain);
      continue;
    }
    mol::SecStructure extended(mol::SecStructure::EXTENDED);
    chain.AssignSecondaryStructure(extended, i->start, i->end);
  }
}

void PDBReader::ClearState()
{
  curr_chain_=mol::ChainHandle();
  curr_residue_=mol::ResidueHandle();
  seqres_=seq::SequenceList();
  chain_count_=0;
  residue_count_=0;
  atom_count_=0;
  hard_end_=false;
  helix_list_.clear();
  strand_list_.clear();
  hets_.clear();
}

bool PDBReader::EnsureLineLength(const StringRef& line, size_t size)
{
  if (line.length()<size) {
    if (profile_.fault_tolerant) {
      return false;
    }    
    throw IOException(str(format("premature end of line %d") %line_num_));
  }  
  return true;
}

bool PDBReader::ParseAtomIdent(const StringRef& line, int line_num, 
                               String& chain_name,  StringRef& res_name,
                               mol::ResNum& resnum, StringRef& atom_name, 
                               char& alt_loc, const StringRef& record_type,
                               int& serial)
{
  if (!this->EnsureLineLength(line, 27)) {
    return false;
  }
  atom_name=line.substr(12, 4).trim();
  if (profile_.calpha_only) {
    if (atom_name!=StringRef("CA", 2)) {
      return false;
    }
  }
  if (charmm_style_) {
    if (line.size()>73) {
      size_t width=std::min(line.size()-72, size_t(4));
      chain_name=line.substr(72, width).trim().str();      
    }
  } else {
    chain_name=String(1, line[21]);    
    if (! restrict_chains_.empty() &&
      restrict_chains_.find(chain_name)==String::npos) {
      return false;
    }    
  }

  alt_loc=line[16];
  res_name=line.substr(17, charmm_style_ ? 4 : 3).trim();
  std::pair<bool, int> res_num=line.substr(22, 4).ltrim().to_int();;
  if (!res_num.first) {
    if (profile_.fault_tolerant) {
      return false;
    }
    throw IOException(str(format("invalid res number on line %d") % line_num));
  }

  char  ins_c=line[26];  
  resnum=to_res_num(res_num.second, ins_c);

  std::pair<bool, int> tmp = line.substr(6, 5).trim().to_int();
  if(!tmp.first) {
    if (profile_.fault_tolerant) {
      return false;
    }
    throw IOException(str(format("invalid atom serial on line %d") % line_num));
  }
  serial = tmp.second;

  return true;
}

void PDBReader::ParseAnisou(const StringRef& line, int line_num, 
                            mol::EntityHandle& ent)
{
  if (!this->EnsureLineLength(line, 77)) {
    return;
  }  
  String chain_name;
  char alt_loc=0;
  StringRef res_name, atom_name;
  mol::ResNum res_num(0);
  int serial = -1;
  if (!this->ParseAtomIdent(line, line_num, chain_name, res_name, res_num, 
                            atom_name, alt_loc, StringRef("ANISOU", 6),
                            serial)) {
    return;                            
  }
  double anisou[6]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  
  for (int i=0;i<6; ++i) {
    std::pair<bool, int> result=line.substr(29+i*7, 6).ltrim().to_int();
    if (!result.first) {
      if (profile_.fault_tolerant) {
        return;
      }
      throw IOException(str(format("invalid ANISOU record on line %d")%line_num));
    }
    anisou[i]=result.second;
  }
  String aname(atom_name.str());
  if (!curr_residue_.IsValid()) {
    if (profile_.fault_tolerant || 
        profile_.calpha_only || 
        profile_.no_hetatms) {
      return;
    } 
    const char* fmt_str="invalid ANISOU record for inexistent atom on line %d";
    throw IOException(str(format(fmt_str) % line_num));      
  }
  mol::AtomHandle atom=curr_residue_.FindAtom(aname);
  if (!atom.IsValid()) {
    if (profile_.fault_tolerant ||
        profile_.calpha_only  ||
        profile_.no_hetatms ||
        warned_name_mismatch_) {
      return;
    }
    const char* fmt_str="invalid ANISOU record for inexistent atom on line %d";
    throw IOException(str(format(fmt_str) % line_num));      
  }
  //get properties which are already set and extend them by adding the ANISOU info
  

  geom::Mat3 mat(anisou[0], anisou[3], anisou[4],
                 anisou[3], anisou[1], anisou[5],
                 anisou[4], anisou[5], anisou[2]);
  mat/=10000;
  atom.SetAnisou(mat);
}

void PDBReader::ParseAndAddAtom(const StringRef& line, int line_num,
                                mol::EntityHandle& ent, 
                                const StringRef& record_type)
{
  if (!this->EnsureLineLength(line, 54)) {
    return;
  }
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  char alt_loc=0;
  String chain_name;
  StringRef res_name, atom_name;
  mol::ResNum res_num(0);
  int serial = -1;
  if (!this->ParseAtomIdent(line, line_num, chain_name, res_name, res_num, 
                            atom_name, alt_loc, record_type, serial)) {
    return;                            
  }
  std::pair<bool, Real> charge, radius;
  std::pair<bool, Real> occ, temp;
  geom::Vec3 apos;

  for (int i=0;i<3; ++i) {
    std::pair<bool, float> result=line.substr(30+i*8, 8).ltrim().to_float();
    if (!result.first) {
      if (profile_.fault_tolerant) {
        return;
      }
      throw IOException(str(format("invalid coordinate on line %d")%line_num));
    }
    apos[i]=result.second;
  }
  // b-factors and occ are replaced by radius and charge if PQR file format
  if(is_pqr_) {
    occ=std::make_pair(true, Real(1.0));
    temp=std::make_pair(true, Real(0.0));
    if (line.length()>=60) {
      charge=line.substr(55,7).ltrim().to_float();
    }
    if (line.length()>=68) {
      radius=line.substr(63,6).ltrim().to_float();
    }
  } else {
    if (line.length()>=60) {
      occ=line.substr(54,6).ltrim().to_float();      
    }
    if (line.length()>=66) {
      temp=line.substr(60, 6).ltrim().to_float();      
    }
  }
  LOG_TRACE( "line: [" << line << "]" );
  String s_ele;
  // determine element from element column (77-78, right justified) if 
  // present otherwise set to empty String. It is up to the builder to 
  // determine the correct element in that case.
  if (line.length()>=78) {                       // element column present
    if(line[76]==' ' && line[77]==' ') {         // both characters are empty
      s_ele="";                                  
    } else if(line[76]!=' ' || line[77]!=' ') {  // at least one character not
                                                 // empty
      if(line[76]==' ' && line[77]!=' ') {       // single character element,
                                                 // right justified
        s_ele=line.substr(77,1).str();           
      } else if(line[76]!=' ' && line[77]==' ') {// single character element,
                                                 // left justified
        s_ele=line.substr(76,1).str();           
      } else {                                   // Real character element
        s_ele=line.substr(76,2).str();            
      }
    }
  }
  String aname(atom_name.str());  
  // some postprocessing
  LOG_TRACE( "s_chain: [" << chain_name << "]" );

  // determine chain and residue update
  bool update_chain=false;
  bool update_residue=false;
  if(!curr_chain_) {
    update_chain=true;
    update_residue=true;
  } else if(curr_chain_.GetName()!=chain_name) {
    update_chain=true;
    update_residue=true;
  }

  if(!curr_residue_) {
    update_residue=true;
  } else if(curr_residue_.GetNumber()!=res_num) {
    update_residue=true;
  }

  if(update_chain) {
    curr_chain_=ent.FindChain(chain_name);
    if(!curr_chain_.IsValid()) {
      LOG_DEBUG("new chain " << chain_name);
      curr_chain_=editor.InsertChain(chain_name);
      ++chain_count_;
    }
    assert(curr_chain_.IsValid());
  }
  if(update_residue) {
    curr_residue_=mol::ResidueHandle();
    if (profile_.join_spread_atom_records) {
      curr_residue_=curr_chain_.FindResidue(res_num);
    }
    if (!curr_residue_.IsValid()) {
      if(curr_chain_.GetResidueCount()>0 && profile_.join_spread_atom_records) {
        int loc=curr_chain_.GetResidueCount()-1;
        for(;loc>=0;--loc) {
          if(curr_chain_.GetResidueByIndex(loc).GetNumber()<res_num) break;
        }
        if(loc<0) {
          curr_residue_=editor.InsertResidueBefore(curr_chain_,0,res_num,res_name.str());
        } else {
          curr_residue_=editor.InsertResidueAfter(curr_chain_,loc,res_num,res_name.str());
        }
        if(!curr_residue_) {
          // this should not happen...
          curr_residue_=editor.AppendResidue(curr_chain_, res_name.str(), res_num);
        }
        LOG_DEBUG("inserted new residue " << res_name << " " << res_num << " after " << loc);
      } else {
        curr_residue_=editor.AppendResidue(curr_chain_, res_name.str(), res_num);
        LOG_DEBUG("appended new residue " << res_name << " " << res_num);
      }
      warned_name_mismatch_=false;
      ++residue_count_; 
    }
    assert(curr_residue_.IsValid());
  }
  // finally add atom
  LOG_DEBUG("adding atom " << aname << " (" << s_ele << " '" << alt_loc << "'" << ") @" << apos);
  mol::AtomHandle ah;
  if (curr_residue_.GetName()!=res_name.str()) {
    if (!profile_.fault_tolerant && alt_loc==' ') {
      std::stringstream ss;
      ss << "error on line " << line_num << ": "
         << "residue with number " << res_num << " has more than one name.";
      throw IOException(ss.str());
    }
    if(!profile_.quack_mode) {
      if (!warned_name_mismatch_) {
        if (alt_loc==' ') {
          LOG_WARNING("Residue with number " << res_num << " has more than one name. "
                      "Ignoring atoms for everything but the first");        
        } else {
          LOG_WARNING("Residue with number " << res_num 
                      << " contains a microheterogeneity. Everything but atoms for "
                      << "the residue '" << curr_residue_.GetName() 
                      << "' will be ignored");
        }
      }
      warned_name_mismatch_=true;
      return;
    }
  }
  Real b=temp.first ? temp.second : 0.0;
  Real o=occ.first ? occ.second : 1.0;
  if (!profile_.quack_mode && alt_loc!=' ') {
    // Check if there is already a atom with the same name.
    mol::AtomHandle me=curr_residue_.FindAtom(aname);
    if (me.IsValid()) {
      try {
        editor.AddAltAtomPos(String(1, alt_loc), me, apos, o, b);
      } catch (Error&) {
        LOG_INFO("Ignoring atom alt location since there is already an atom "
                     "with name " << aname << ", but without an alt loc");
        return;
      }
      return;
    } else {
      ah=editor.InsertAltAtom(curr_residue_, aname,
                              String(1, alt_loc), apos, s_ele, o, b);
      ++atom_count_;
    }
  } else {
    mol::AtomHandle atom=curr_residue_.FindAtom(aname);
    if (atom.IsValid() && !profile_.quack_mode) {
      if (profile_.fault_tolerant) {
        LOG_WARNING("duplicate atom '" << aname << "' in residue " 
                    << curr_residue_);
        return;
      }
      throw IOException("duplicate atom '"+aname+"' in residue "+
                        curr_residue_.GetQualifiedName());
    }
    ah=editor.InsertAtom(curr_residue_, aname, apos, s_ele);
    ++atom_count_;
  }
  if(is_pqr_) {
    if (radius.first) {
      ah.SetRadius(radius.second);
    }
  }
  ah.SetBFactor(b);
  ah.SetOccupancy(o);
  if (charge.first) {
    ah.SetCharge(charge.second);
  }
  ah.SetHetAtom(record_type[0]=='H');
  if(profile_.read_conect && serial != -1) {
    this->amap_[serial] = ah;
  }
}

void PDBReader::ParseHelixEntry(const StringRef& line)
{
  if (!this->EnsureLineLength(line, 38)) {
    return;
  }
  std::pair<bool, int>  start_num=line.substr(21, 4).ltrim().to_int();
  std::pair<bool, int>  end_num=line.substr(33, 4).ltrim().to_int();
  if (!start_num.first || !end_num.first) {
    if (profile_.fault_tolerant) {
      return;
    }    
    throw IOException(str(format("invalid helix entry on line %d") % line_num_));
  }
  LOG_DEBUG("making helix entry: " << start_num.second << ", " 
             << line[25] << " " << end_num.second << " " << line[37]);
  HSEntry hse = {to_res_num(start_num.second, line[25]),
                 to_res_num(end_num.second, line[37]),
                 line.substr(19,1).str()};
  if (restrict_chains_.size()==0 ||
      restrict_chains_.find(hse.chain)!=String::npos) {
    helix_list_.push_back(hse);
  }

}

void PDBReader::ParseStrandEntry(const StringRef& line)
{
  if (!this->EnsureLineLength(line, 38)) {
    return;
  }  
  std::pair<bool, int>  start_num=line.substr(22, 4).ltrim().to_int();
  std::pair<bool, int>  end_num=line.substr(33, 4).ltrim().to_int();
  if (!start_num.first || !end_num.first) {
    if (profile_.fault_tolerant) {
      return;
    }
    throw IOException(str(format("invalid strand entry on line %d")%line_num_));
  }
  LOG_DEBUG("making strand entry: " << start_num.second << ", " << line[26] 
            << " " << end_num.second << " " << line[37] << " chain=" << line[21]);
  HSEntry hse = {to_res_num(start_num.second, line[26]),
                 to_res_num(end_num.second, line[37]),
                 line.substr(21,1).str()};
  if (restrict_chains_.size()==0 ||
     restrict_chains_.find(hse.chain)!=String::npos) {
   strand_list_.push_back(hse);
  }
}

void PDBReader::ParseConectEntry (const StringRef& line, int line_num, mol::EntityHandle& ent)
{
  if (line.size()<16) {
    if (profile_.fault_tolerant) {
      LOG_WARNING("invalid CONECT record on line " << line_num 
                  << ": record is too short");
      return;
    }
    std::stringstream ss;
    ss << "invalid CONECT record on line ";
    ss << line_num <<": record is too short";
    throw IOException(ss.str());
  }
  if (line.rtrim().size()>80) {
    if (profile_.fault_tolerant) {
      LOG_WARNING("invalid CONECT record on line " << line_num 
                  << ": record is too long");
      return;
    }
    std::stringstream ss;
    ss << "invalid CONECT record on line ";
    ss << line_num <<": whole record is too long";
    throw IOException(ss.str());
  }
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  const int starting_atom = line.substr(6,5).trim().to_int().second;
  // map for bonds, in the form of <serial_number, bond_order>
  std::map<int, unsigned char> connected_atoms;
  for (int i=0; i<4; i++) {
    if (static_cast<int>(line.length()) < 11+i*5+5) break;
    const int connected_atom = line.substr(11+i*5, 5).trim().to_int().second;
    connected_atoms[connected_atom]+=1;
  }
  for (auto& pair : connected_atoms) {
    auto at_one_it = this->amap_.find(starting_atom);
    auto at_two_it = this->amap_.find(pair.first);
    if(at_one_it != this->amap_.end() && at_two_it != this->amap_.end()) {
      editor.Connect(at_one_it->second, at_two_it->second, pair.second);
    }
  }
}


}}
