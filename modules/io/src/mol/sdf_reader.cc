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
  Author: Tobias Schmidt
 */

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/format.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/lexical_cast.hpp>
#include <ost/mol/bond_handle.hh>
#include <ost/conop/conop.hh>
#include <ost/io/io_exception.hh>
#include <ost/log.hh>
#include <ost/mol/xcs_editor.hh>

#include "sdf_reader.hh"

namespace ost { namespace io {

using boost::format;

SDFReader::SDFReader(const String& filename)
  : infile_(filename), instream_(infile_)
{
  this->ClearState(boost::filesystem::path(filename));
}

SDFReader::SDFReader(const boost::filesystem::path& loc)
  : infile_(loc), instream_(infile_)
{
  this->ClearState(loc);
}

SDFReader::SDFReader(std::istream& instream)
  : infile_(), instream_(instream)
{
  this->ClearState(boost::filesystem::path(""));
}

// import data from provided stream
void SDFReader::Import(mol::EntityHandle& ent)
{
  String line;
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  while (std::getline(in_,line)) {
    ++line_num;

    // std::getline removes EOL character but may leave a DOS CR (\r) in Unix
    size_t cr_pos = line.find("\r");
    if (cr_pos != String::npos) {
        LOG_TRACE( "Remove CR@" << cr_pos);
        line.erase(cr_pos);
    }

    if (line_num<=4) {
      ParseAndAddHeader(line, line_num, ent, editor);
    } else if (line_num<=atom_count_+4) {
      ParseAndAddAtom(line, line_num, ent, true, editor);
    } else if (line_num<=bond_count_+atom_count_+4) {
      ParseAndAddBond(line, line_num, ent, editor);
    } else if (boost::iequals(line.substr(0,2), "> ")) {
      // parse data items
      int data_header_start = line.find('<');
      if(data_header_start==-1) {
        String msg="Bad data line %d: Could not find start or end of header";
        throw IOException(str(format(msg) % line_num));
      }
      String data_header = line.substr(data_header_start+1,line.rfind('>')-data_header_start-1);
      if(data_header.empty()) {
        String msg="Bad data line %d: Could not find data header";
        throw IOException(str(format(msg) % line_num));
      }
      String data_value="";
      while(std::getline(in_,line) && !boost::iequals(line, "")) {
        data_value.append(line);
      }
      curr_chain_.SetStringProp(data_header, data_value);
    } else if (boost::iequals(line, "$$$$")) {
      LOG_VERBOSE("MOLECULE " << curr_chain_.GetName() << " (" << chain_count_ << ") added.")
      NextMolecule();
    }
  }

  LOG_INFO("imported " << chain_count_ << " chains, " << residue_count_
               << " residues, " << atom_count_ << " atoms");
}

void SDFReader::ClearState(const boost::filesystem::path& loc)
{
  if (boost::iequals(".gz", boost::filesystem::extension(loc))) {
    in_.push(boost::iostreams::gzip_decompressor());
  }
  in_.push(instream_);
  if(!infile_) throw IOException("could not open "+loc.string());
  curr_chain_=mol::ChainHandle();
  curr_residue_=mol::ResidueHandle();
  chain_count_=0;
  residue_count_=0;
  atom_count_=0;
  bond_count_=0;
  line_num=0;
}

void SDFReader::NextMolecule()
{
  residue_count_=0;
  atom_count_=0;
  bond_count_=0;
  line_num=0;
  curr_residue_ = ost::mol::ResidueHandle();
  curr_chain_ = ost::mol::ChainHandle();
}

void SDFReader::ParseAndAddHeader(const String& line, int line_num,
                                  mol::EntityHandle& ent, mol::XCSEditor& editor)
{
  LOG_TRACE( "line: [" << line << "]" );
  format chain_fmter("%05i_%s");
  switch(line_num)
  {
    case 1:  // title line
    {
      ++chain_count_;
      String s_title=line;
      String s_chain;
      chain_fmter % chain_count_ % boost::trim_copy(s_title);
      s_chain=chain_fmter.str();
      if(s_chain.empty()) {
        String msg="Bad molecule name line %d: Line is empty";
        throw IOException(str(format(msg) % line_num));
      }
      // prepeare required variables to add new chain and residue
      // once we parse the first atom
      curr_chain_name_ = s_chain;
      curr_res_key_ = boost::trim_copy(s_title);
      break;
    }
    case 2:  // user information line
      break;
    case 3:  // comments line
      break;
    case 4:  // counts line
    {
      String s_anum=line.substr(0,3);
      try {
        atom_count_=boost::lexical_cast<int>(boost::trim_copy(s_anum));
      } catch(boost::bad_lexical_cast&) {
        String msg="Bad counts line %d: Can't convert number of atoms"
                   " '%s' to integral constant.";
        throw IOException(str(format(msg) % line_num % s_anum));
      }
      String s_bnum=line.substr(3,3);
      try {
        bond_count_=boost::lexical_cast<int>(boost::trim_copy(s_bnum));
      } catch(boost::bad_lexical_cast&) {
        String msg="Bad counts line %d: Can't convert number of bonds"
                   " '%s' to integral constant.";
        throw IOException(str(format(msg) % line_num % s_bnum));
      }
      break;
    }
  }
}

void SDFReader::ParseAndAddAtom(const String& line, int line_num,
                                mol::EntityHandle& ent, bool hetatm,
                                mol::XCSEditor& editor)
{

  LOG_TRACE( "line: [" << line << "]" );

  if(line.length()<48 || line.length()>69) {
    // Handle the case where we have trailing space characters
    if (line.length()>69 && boost::trim_copy(line.substr(69)) == "") {
      LOG_DEBUG( "Ignoring trailing space" );
    }
    else {
      String msg="Bad atom line %d: Not correct number of characters on the"
                 " line: %i (should be between 48 and 69)";
      throw IOException(str(format(msg) % line_num % line.length()));
    }
  }
  int anum = line_num-4;  // start at 1 on fifth line since first four lines are header
  String s_posx=line.substr(0,10);
  String s_posy=line.substr(10,10);
  String s_posz=line.substr(20,10);
  String s_ele=line.substr(31,3);
  String s_charge=line.substr(36,3);

  geom::Vec3 apos;
  try {
    apos=geom::Vec3(boost::lexical_cast<Real>(boost::trim_copy(s_posx)),
                    boost::lexical_cast<Real>(boost::trim_copy(s_posy)),
                    boost::lexical_cast<Real>(boost::trim_copy(s_posz)));
  } catch(boost::bad_lexical_cast&) {
    String msg="Bad atom line %d: Can't convert coordinates to "
               "floating point numbers";
    throw IOException(str(format(msg) % line_num));
  }

  String ele=boost::trim_copy(s_ele);
  String upper_ele=ele;
  std::transform(upper_ele.begin(),upper_ele.end(),upper_ele.begin(),toupper);
  String aname=boost::lexical_cast<String>(anum);
  
  Real charge=0.0;  
  try {
    charge=boost::lexical_cast<Real>(boost::trim_copy(s_charge));
    if (charge!=0) {
      charge=4-charge;
    } //4-sdf_charge=real_charge if not 0
  } catch(boost::bad_lexical_cast&) {
    String msg="Bad atom line %d: Can't convert charge"
               " '%s' to integral constant.";
    throw IOException(str(format(msg) % line_num % s_charge));
  }

  if(!curr_chain_.IsValid()) {
      curr_chain_=editor.InsertChain(curr_chain_name_);
      LOG_DEBUG("new chain " << curr_chain_name_);
  }

  if(!curr_residue_.IsValid()) {
      mol::ResNum rnum(++residue_count_);
      curr_residue_=editor.AppendResidue(curr_chain_, curr_res_key_, rnum);
      LOG_DEBUG("new residue " << curr_res_key_ << "(" << rnum << ")");
  }

  LOG_DEBUG("adding atom " << aname << " (" << s_ele << ") @" << apos);

  mol::AtomHandle atom=editor.InsertAtom(curr_residue_, aname, apos, upper_ele);
  atom.SetHetAtom(hetatm);
  atom.SetCharge(charge);
}


void SDFReader::ParseAndAddBond(const String& line, int line_num,
                                mol::EntityHandle& ent, mol::XCSEditor& editor)
{

  LOG_TRACE( "line: [" << line << "]" );

  if(line.length()<9 || line.length()>21) {
    // Handle the case where we have trailing space characters
    if (line.length()>21 && boost::trim_copy(line.substr(21)) == "") {
      LOG_DEBUG( "Ignoring trailing space" );
    }
    else {
      String msg="Bad bond line %d: Not correct number of characters on the"
                 " line: %i (should be between 9 and 21)";
      throw IOException(str(format(msg) % line_num % line.length()));
    }
  }

  String s_first_name=line.substr(0,3);
  String s_second_name=line.substr(3,3);
  String s_type=line.substr(6,3);
  String first_name, second_name;
  unsigned char type;
  mol::BondHandle bond;

  first_name=boost::trim_copy(s_first_name);
  second_name=boost::trim_copy(s_second_name);

  try {
    type=boost::lexical_cast<int>(boost::trim_copy(s_type));
    if (type<1 || type>8) {
      String msg="Bad bond line %d: Bond type number"
                      " '%s' not within accepted range (1-8).";
      throw IOException(str(format(msg) % line_num % s_type));
    }
  } catch(boost::bad_lexical_cast&) {
    String msg="Bad bond line %d: Can't convert bond type number"
                " '%s' to integral constant.";
    throw IOException(str(format(msg) % line_num % s_type));
  }

  mol::AtomHandle first,second;

  first = ent.FindAtom(curr_chain_.GetName(), mol::ResNum(residue_count_), 
                       first_name);
  second = ent.FindAtom(curr_chain_.GetName(), mol::ResNum(residue_count_), 
                        second_name);

  if (first.IsValid() && second.IsValid()) {
    bond = editor.Connect(first, second);
    bond.SetBondOrder(type);
  } else {
    String msg="Bad bond line %d: Can't find the atom names '%s', '%s'"
               " in entity.";
    throw IOException(str(format(msg) % line_num % first % second));
  }

  LOG_DEBUG("adding bond " << s_first_name << " " << s_second_name << " (" 
            << s_type << ") ");
}

}}
