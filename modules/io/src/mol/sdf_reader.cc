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

boost::iostreams::filtering_stream<boost::iostreams::input>& SDFReader::GetLine(
  boost::iostreams::filtering_stream<boost::iostreams::input>& in,
  String& line)
  // Read next line from in and place it in line.
  // Remove trailing \r characters.
{
  std::getline(in, line);
  size_t cr_pos = line.find("\r");
  if (cr_pos != String::npos) {
      LOG_TRACE( "Remove CR@" << cr_pos);
      line.erase(cr_pos);
  }
  return in;
}

// import data from provided stream
void SDFReader::Import(mol::EntityHandle& ent)
{
  String line;
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  while (GetLine(in_,line)) {
    ++line_num;

    if (line_num<=4) {
      ParseHeader(line, line_num, ent, editor);
    } else if (version_ == "V2000" && line_num<=atom_count_+4) {
      AddAtom(ParseAtom(line, line_num), line_num, ent, true, editor);
    } else if (version_ == "V2000" && line_num<=bond_count_+atom_count_+4) {
      AddBond(ParseBond(line, line_num), line_num, ent, editor);
    } else if (version_ == "V2000" &&  boost::iequals(line.substr(0,6), "M  CHG")) {
      AddCharge(ParseMCharge(line, line_num), line_num, ent, editor);
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
      while(GetLine(in_,line) && !boost::iequals(line, "")) {
        data_value.append(line);
      }
      curr_chain_.SetStringProp(data_header, data_value);
    } else if (boost::iequals(line, "$$$$")) {
      LOG_VERBOSE("MOLECULE " << curr_chain_.GetName() << " (" << chain_count_ << ") added.")
      NextMolecule();
    } else if (version_ == "V3000") {
      ProcessV3000Line(line, ent, editor);
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
  version_="";
  v3000_bond_block_=false;
  v3000_atom_block_=false;
  charges_reset_=false;
}

void SDFReader::NextMolecule()
{
  residue_count_=0;
  atom_count_=0;
  bond_count_=0;
  line_num=0;
  version_="";
  v3000_bond_block_=false;
  v3000_atom_block_=false;
  charges_reset_=false;
  curr_residue_ = ost::mol::ResidueHandle();
  curr_chain_ = ost::mol::ChainHandle();
}

void SDFReader::ParseHeader(const String& line, int line_num,
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
      if (line.length() < 39) {
        String msg="Bad counts line %d: Not correct number of characters on "
                   "the line: %i (should be at least 39)";
        throw IOException(str(format(msg) % line_num % line.length()));
      }
      String version_str=line.substr(34, 5);
      if (version_str == "V2000" || version_str == "V3000") {
        version_=version_str;
      }
      else {
        String msg="Unsupported SDF version: %s.";
        throw IOException(str(format(msg) % version_str));
      }
      // Counts will be overridden in V3000
      String s_anum=line.substr(0,3);
      String s_bnum=line.substr(3,3);
      SetCounts(s_anum, s_bnum, line_num);
      break;
    }
  }
}

void SDFReader::SetCounts(const String& anum, const String bnum, int line_num)
{
  try {
    atom_count_=boost::lexical_cast<int>(boost::trim_copy(anum));
  } catch(boost::bad_lexical_cast&) {
    String msg="Bad counts line %d: Can't convert number of atoms"
               " '%s' to integral constant.";
    throw IOException(str(format(msg) % line_num % anum));
  }
  try {
    bond_count_=boost::lexical_cast<int>(boost::trim_copy(bnum));
  } catch(boost::bad_lexical_cast&) {
    String msg="Bad counts line %d: Can't convert number of bonds"
               " '%s' to integral constant.";
    throw IOException(str(format(msg) % line_num % bnum));
  }
}

SDFReader::atom_data SDFReader::ParseAtom(const String& line, int line_num)
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

  return std::make_tuple(anum, s_posx, s_posy, s_posz, s_ele, s_charge);
}

void SDFReader::AddAtom(const atom_data& atom_tuple, int line_num, mol::EntityHandle& ent,
                        bool hetatm, mol::XCSEditor& editor)
{
  int anum;
  String s_posx, s_posy, s_posz, s_ele, s_charge;
  tie(anum, s_posx, s_posy, s_posz, s_ele, s_charge) = atom_tuple;

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


SDFReader::bond_data SDFReader::ParseBond(const String& line, int line_num)
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

  return std::make_tuple(s_first_name, s_second_name, s_type);
}

void SDFReader::AddBond(const bond_data& bond_tuple, int line_num, mol::EntityHandle& ent,
                        mol::XCSEditor& editor)
{
  String s_first_name, s_second_name, s_type;
  tie(s_first_name, s_second_name, s_type) = bond_tuple;

  unsigned char type;
  mol::BondHandle bond;

  String first_name, second_name;
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


void SDFReader::ResetCharges()
// from doc of V2000 Atom Block:
// > Retained for compatibility with older Ctabs, M CHG and M RAD lines take
// > precedence.
// Therefore we must reset all charges of the residue if we encounter an
// M  CHG line.
{
  LOG_DEBUG("Resetting all charges to 0.");
  for (mol::AtomHandle & atom : curr_residue_.GetAtomList()) {
    atom.SetCharge(0.0);
  }
  charges_reset_=true;
}


SDFReader::charge_data SDFReader::ParseMCharge(const String& line, int line_num)
{

  LOG_TRACE( "line: [" << line << "]" );

  if (!charges_reset_) {
    ResetCharges();
  }

  if(line.length()<15 || line.length()>17) {
    // Handle the case where we have trailing space characters
    if (line.length()>17 && boost::trim_copy(line.substr(17)) == "") {
      LOG_DEBUG( "Ignoring trailing space" );
    }
    else {
      String msg="Bad Charge line %d: Not correct number of characters on the"
                 " line: %i (should be between 15 and 17)";
      throw IOException(str(format(msg) % line_num % line.length()));
    }
  }

  String atom_index=line.substr(10,3);
  String charge=line.substr(14,3);

  return std::make_tuple(atom_index, charge);
}

  void SDFReader::AddCharge(const charge_data& charge_tuple, int line_num, mol::EntityHandle& ent,
                        mol::XCSEditor& editor)
{
  String s_atom_index, s_charge;
  tie(s_atom_index, s_charge) = charge_tuple;

  int atom_index;
  Real charge;

  try {
    atom_index=boost::lexical_cast<int>(boost::trim_copy(s_atom_index));
    if (atom_index > atom_count_) {
      String msg="Bad charge line %d: Atom index"
                      " '%d' greater than number of atoms in the molecule (%d).";
      throw IOException(str(format(msg) % line_num % atom_index % atom_count_));
    } else if (atom_index < 1) {
      String msg="Bad charge line %d: Atom index %d < 1.";
      throw IOException(str(format(msg) % line_num % atom_index));
    }
  } catch(boost::bad_lexical_cast&) {
    String msg="Bad charge line %d: Can't convert atom index"
                " '%s' to integral constant.";
    throw IOException(str(format(msg) % line_num % s_atom_index));
  }

  try {
    charge=boost::lexical_cast<Real>(boost::trim_copy(s_charge));
  } catch(boost::bad_lexical_cast&) {
    String msg="Bad charge line %d: Can't convert charge"
                " '%s' to real number.";
    throw IOException(str(format(msg) % line_num % s_charge));
  }

  curr_residue_.GetAtomList()[atom_index - 1].SetCharge(charge);

  LOG_DEBUG("Setting charge of atom " << atom_index - 1 << " to " << charge);
}

SDFReader::v3000_line_tokens SDFReader::TokenizeV3000Line(const String& line,
                                                          int line_num,
                                                          int num_posval)
// Read whitespace-separated tokens from a V3000 line.
// Tokens can be separated by any amount of whitespace.
// The function is guaranteed to return exactly num_posval positional elements,
// or throws an error. It returns any number of keyword elements with only
// syntax checks (ie no checks if the keywords are correct, only well-formed).
{
  std::istringstream v30_stream(line);
  std::vector<String> positional;
  std::map<String, String> keywords;
  String token;
  bool keywords_reached = false;
  size_t kw_equal_pos;
  positional.reserve(num_posval);

  while (v30_stream.tellg() != -1) {
    std::getline(v30_stream, token, ' ');
    if (token.empty()) {
      continue;
    }
    kw_equal_pos = token.find('=');
    if (kw_equal_pos != String::npos) {
      keywords_reached = true;
    }
    if (keywords_reached) {
      // Token can contain a list in round brackets
      // We don't use them in OST so no fancy parsing, just capture them
      // as a string keyword
      if (token.find('(') == kw_equal_pos + 1) {
        // Search for the closing bracket
        while (token.find(')') == String::npos) {
          String next_token;
          std::getline(v30_stream, next_token, ' ');
          token = token + " " + next_token;
        }
      }

      // Check if keyword is well formed
      if (token.size() < 3 // too short
          || kw_equal_pos == String::npos // no =
          || kw_equal_pos == 0 // no key (starts with =)
          || kw_equal_pos == token.size() - 1 // no value (ends with =)
          ) {
        String msg="Bad V3000 keyword on line %d: '%s'.";
        throw IOException(str(format(msg) % line_num % token));
      }
      String key = token.substr(0, kw_equal_pos);
      String value = token.substr(kw_equal_pos + 1);
      keywords.insert({key, value});
    }
    else {
      positional.push_back(token);
    }
  }

  int obtained_posval = positional.size();
  if (obtained_posval != num_posval) {
    String msg="Bad V3000 line %d: expected %d positional values, got %d.";
    throw IOException(str(format(msg) % line_num % num_posval %
                          obtained_posval));
  }

  return std::make_tuple(positional, keywords);
}

String SDFReader::CleanupV3000Line(const String& line)
// String cleanup and aggregation for V3000
// Return a string with no "M  V30 " and not ending with -
{
  String v30_line = line;
  if (v30_line.substr(0, 7) != "M  V30 ") {
    String msg="Bad V3000 line %d: starts with '%s'.";
    throw IOException(str(format(msg) % line_num % line.substr(0, 6)));
  }

  // Handle line continuation character -
  while (v30_line.find("-") == v30_line.length()-1) {
    // Read and append the next line
    String next_line;
    GetLine(in_,next_line);
    ++line_num; // Update class member

    // Ensure we have a valid next_line
    if (next_line.substr(0, 7) != "M  V30 ") {
      String msg="Bad V3000 line %d: starts with '%s'.";
      throw IOException(str(format(msg) % line_num % next_line.substr(0, 6)));
    }
    // All clear, add data
    v30_line = v30_line.erase(v30_line.find("-")) + next_line.substr(7);
    LOG_TRACE( "V3000 line: [" << v30_line << "]" );
  }

  // Cleanup the line
  return v30_line.substr(7); // We previously ensured it starts with M  V30
}

SDFReader::atom_data SDFReader::ParseV3000Atom(const String& line, int line_num)
{
  v3000_line_tokens tokens = TokenizeV3000Line(line, line_num, 6);
  std::vector<String> posval;
  std::map<String, String> keywords;
  tie(posval, keywords) = tokens;

  String s_anum = posval[0];
  String atype = posval[1];
  String posx = posval[2];
  String posy = posval[3];
  String posz = posval[4];

  String chg;
  try {
    chg = keywords.at("CHG");
  } catch(std::out_of_range&) {
    chg = "0";
  }

  int anum;
  try {
    anum=boost::lexical_cast<int>(boost::trim_copy(s_anum));
  } catch(boost::bad_lexical_cast&) {
    String msg="Bad atom index '%s' on line %d.";
    throw IOException(str(format(msg) % s_anum % line_num));
  }

  return std::make_tuple(anum, posx, posy, posz, atype, chg);
}

SDFReader::bond_data SDFReader::ParseV3000Bond(const String& line, int line_num)
{
  v3000_line_tokens tokens = TokenizeV3000Line(line, line_num, 4);
  std::vector<String> posval;
  tie(posval, std::ignore) = tokens;

  String btype = posval[1];
  String s_first_name = posval[2];
  String s_second_name = posval[3];

  return std::make_tuple(s_first_name, s_second_name, btype);
}

std::tuple<String, String> SDFReader::ParseV3000Counts(const String& line, int line_num)
{
  v3000_line_tokens tokens = TokenizeV3000Line(line, line_num, 5);
  std::vector<String> posval;
  tie(posval, std::ignore) = tokens;

  String anum = posval[0];
  String bnum = posval[1];

  return std::make_tuple(anum, bnum);
}

void SDFReader::VerifyV3000Counts()
{
  int actual_atom_count = curr_residue_.GetAtomCount();
  int actual_bond_count = curr_residue_.GetBondCount();
  if (actual_atom_count != atom_count_) {
      String msg="Bad counts for molecule ending on line %d: "
                 "expected %d atoms, got %d.";
      throw IOException(str(format(msg) % line_num % atom_count_ %
                        actual_atom_count));
    }
  if (actual_bond_count != bond_count_) {
      String msg="Bad counts for molecule ending on line %d: "
                 "expected %d bonds, got %d.";
      throw IOException(str(format(msg) % line_num % bond_count_ %
                        actual_bond_count));
    }
}

void SDFReader::ProcessV3000Line(const String& line,
                                 mol::EntityHandle& ent,
                                 mol::XCSEditor& editor)
{
  if (line.substr(0, 6) == "M  END") {
    VerifyV3000Counts();
    return;
  }
  String v30_line = CleanupV3000Line(line);

  if (v30_line.substr(0, 6) == "COUNTS") {
    String anum, bnum;
    std::tie(anum, bnum) = ParseV3000Counts(v30_line.substr(7), line_num);
    SetCounts(anum, bnum, line_num);
  }
  else if (v30_line.substr(0, 10) == "BEGIN ATOM") {
    v3000_atom_block_=true;
  }
  else if (v30_line.substr(0, 8) == "END ATOM") {
    v3000_atom_block_=false;
  }
  else if (v30_line.substr(0, 10) == "BEGIN BOND") {
    v3000_bond_block_=true;
  }
  else if (v30_line.substr(0, 8) == "END BOND") {
    v3000_bond_block_=false;
  }
  else if (v3000_atom_block_) {
    AddAtom(ParseV3000Atom(v30_line, line_num), line_num, ent, true, editor);
  }
  else if (v3000_bond_block_) {
    AddBond(ParseV3000Bond(v30_line, line_num), line_num, ent, editor);
  }
  else {
    LOG_TRACE( "ignoring line: [" << v30_line << "]" );
  }
}

}}
