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
  Author: Marco Biasini
 */
#include <locale>

#include <string.h>

#include <boost/format.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/algorithm/string.hpp>

#include <ost/base.hh>
#include <ost/boost_filesystem_helper.hh>
#include <ost/io/io_exception.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/mol/entity_visitor.hh>
#include <ost/mol/bond_handle.hh>
#include "pdb_writer.hh"

using boost::format;

namespace ost { namespace io {


namespace {

// determine whether the element name has to be shifted to to the left by one
// column.
bool shift_left(const String& atom_name, bool is_hetatm, 
                const String& element, float mass) 
{
  if (atom_name.length()==4) {
    return true;
  }
  if (is_hetatm==false) {
    return false;
  }

  if (isdigit(atom_name[0]) || atom_name=="UNK") {
    return true;
  }
  if (mass>34) {
    if (element=="W" || element=="V") {
      return false;
    }    
    return true;
  }
  return (element=="K" || element=="CA" || element=="NA" || 
          element=="MG" || element=="LI");
}


bool atom_pos_ok(geom::Vec3 p) {
  for (int i=0; i<3; ++i) {
    if (p[i]<=-1000 || p[i] >= 10000) {
      return false;
    }
  }
  return true;
}

// write the serial field for ATOM and TER records, use '*****' for 
// indexes > 99999 (borrowed from charmm)
void write_serial(int serial, FormattedLine& line) {
  // Avoid writing out atomnumbers larger than 5 digits
  if (serial > 99999) {
    line( 6, 5)=fmt::LPadded("*****");
  } else {
    line( 6, 5)=fmt::LPaddedInt(serial);
  }
}

void write_atom(std::ostream& ostr, FormattedLine& line, 
                const mol::AtomHandle& atom, int atomnum, 
                bool is_pqr, bool charmm_style)
{
  mol::ResidueHandle res=atom.GetResidue();
  char ins_code=res.GetNumber().GetInsCode();
  StringRef record_name(atom.IsHetAtom() ? "HETATM" : "ATOM  ", 6);
  std::vector<String> names=atom.GetAltGroupNames();
  
  geom::Vec3 p=atom.GetPos();
  line( 0, 6)=record_name;
  write_serial(atomnum, line);
  String atom_name=atom.GetName();
  if (atom_name.size()>4) {
    throw IOException("Atom name '"+atom.GetQualifiedName()+
                      "' is too long for PDB output. At most 4 character "
                      "are allowed");
  }
  if (shift_left(atom_name, atom.IsHetAtom(), atom.GetElement(), 
                 atom.GetMass())) {
    line(12, 4)=fmt::RPadded(atom_name);
  } else {
    line(13, 3)=fmt::RPadded(atom_name);
  }
  if (charmm_style) {
    if (res.GetKey().size()>4) {
      throw IOException("Residue name '"+res.GetName()+
                        "' is too long for CHARMM-PDB output. At most 4 "
                        "characters are allowed");
    }
    line(17, 4)=fmt::RPadded(res.GetKey());
  } else {
    if (res.GetKey().size()>3) {
      throw IOException("Residue name '"+res.GetName()+
                        "' is too long for PDB output. At most 3 characters "
                        "are allowed");
    }
    line(17, 3)=fmt::LPadded(res.GetKey());    
  }

  
  String chain_name=res.GetChain().GetName();
  if (!charmm_style) {
      if (! chain_name.empty()) {
      if (chain_name.size()==1) {
        line[21]=chain_name[0];
      } else {
        throw IOException("PDB format only support single-letter chain names");
      }
    }    
  }

  int num = res.GetNumber().GetNum();
  if(num < -999 || num > 9999) {
    throw IOException("Residue number from " + res.GetQualifiedName() +
                      " is out of range supported by the PDB format " 
                      "(-999 to 9999)");
  }
  line(22, 4)=fmt::LPaddedInt(num);
  if (ins_code!=0) {
    line[26]=ins_code;
  }
  
  // deal with alternative atom locations
  if (names.empty()) {
    if (!atom_pos_ok(p)) {
      throw IOException("Atom position outside of bounds supported by PDB format");
    }
    line(30, 8)=fmt::LPaddedFloat(p[0],  3);
    line(38, 8)=fmt::LPaddedFloat(p[1],  3);
    line(46, 8)=fmt::LPaddedFloat(p[2],  3);
    
    if (is_pqr) {
      line(55, 7)=fmt::LPaddedFloat(atom.GetCharge(), 4);
      line(63, 6)=fmt::LPaddedFloat(atom.GetRadius(), 4);
    } else {
      line(54, 6)=fmt::LPaddedFloat(atom.GetOccupancy(), 2);
      Real bfac=atom.GetBFactor();
      if (bfac>999.99) {
        line(60, 6)=fmt::LPaddedFloat(999.99, 2);
      } else {
        line(60, 6)=fmt::LPaddedFloat(bfac, 2);
      }
    }
    if (charmm_style) {
      line(72, 4)=fmt::RPadded(chain_name);
    }
    line(76, 2)=fmt::LPadded(atom.GetElement());
    ostr << line;
  } else {
    for (std::vector<String>::const_iterator
         i=names.begin(), e=names.end(); i!=e; ++i) {
      // GetAltPos always return orig pos, i.e. does not honor
      // transformations like GetPos does - so apply it here
      // seems like a FIXME...
      p=atom.GetEntity().GetTransform().Apply(atom.GetAltPos(*i));
      line(30, 50).Clear();

      if (!atom_pos_ok(p)) {
        throw IOException("Atom position outside of bounds supported by PDB format");
      }
      if (i->size()>1) {
        throw IOException("Alternative atom indicator '"+atom.GetQualifiedName()+
                          "("+*i+")' too long for PDB output. At most 1 "
                          "character is allowed");
      }
      line[16]=(*i)[0];
      line(30, 8)=fmt::LPaddedFloat(p[0],  3);
      line(38, 8)=fmt::LPaddedFloat(p[1],  3);
      line(46, 8)=fmt::LPaddedFloat(p[2],  3);

      if (is_pqr) {
       line(54, 6)=fmt::LPaddedFloat(atom.GetCharge(), 2);
       line(60, 6)=fmt::LPaddedFloat(atom.GetRadius(), 2);
      } else {
       line(54, 6)=fmt::LPaddedFloat(atom.GetAltOcc(*i), 2);
       Real bfac=atom.GetAltBFactor(*i);
       if (bfac>999.99) {
         line(60, 6)=fmt::LPaddedFloat(999.99, 2);
       } else {
         line(60, 6)=fmt::LPaddedFloat(bfac, 2);
       }
      }
      if (charmm_style) {
        line(72, 4)=fmt::RPadded(chain_name);
      }
      line(76, 2)=fmt::LPadded(atom.GetElement());
      ostr << line;
    }
  }

  line.Clear();
}

void write_conect(std::ostream& ostr, int atom_index,
                  std::list<int>& partner_indices)
{
  int counter=0;
  partner_indices.sort();
  for (std::list<int>::const_iterator i=partner_indices.begin();
       i!=partner_indices.end(); ++i, counter++) {
    if (counter>3) {
      ostr << std::endl;
      counter=0;
    }
    if (counter==0) ostr << "CONECT" << format("%5i") % atom_index;
    ostr << format("%5i") % (*i);
  }
  ostr << std::endl;
}

class PDBWriterImpl : public mol::EntityVisitor {
public:
  PDBWriterImpl(std::ostream& ostream, FormattedLine& line, 
                std::map<long,int>& atom_indices, bool charmm_style)
    : ostr_(ostream), counter_(0), is_pqr_(false),
      atom_indices_(atom_indices), line_(line), peptide_(false),
      charmm_style_(charmm_style) {
  }
private:
public:
  virtual bool VisitAtom(const mol::AtomHandle& atom) {
    counter_++;
    write_atom(ostr_, line_, atom, counter_, is_pqr_, charmm_style_);
    if (atom.IsHetAtom()) {
      atom_indices_[atom.GetHashCode()]=counter_;
    }
    return true;
  }
  
  virtual bool VisitResidue(const mol::ResidueHandle& res)
  {
    if (res.IsPeptideLinking()) {
      peptide_=true;
    } else {
      if (peptide_) {
        this->WriteTer(prev_);
      }
      peptide_=false;
    }
    prev_=res;
    return true;
  }
  
  virtual bool VisitChain(const mol::ChainHandle& chain)
  {
    if (chain.GetResidueCount()!=0) {
      if (peptide_) {
        this->WriteTer(prev_);
      }
      peptide_=false;
    }
    return true;
  }

  virtual void OnExit()
  {
    if (peptide_) {
      this->WriteTer(prev_);
    }
  }

  void WriteTer(mol::ResidueHandle res)
  {
    counter_++;
    line_(0, 6)=StringRef("TER   ", 6);
    write_serial(counter_, line_);
    if(charmm_style_){
      line_(17, 4)=fmt::RPadded(res.GetKey());
    }
    else{
      line_(17, 3)=fmt::LPadded(res.GetKey()); 
    }
    if (!res.GetChain().GetName().empty() && !charmm_style_) {
      line_[21]=res.GetChain().GetName()[0];
    }

    line_(22, 4)=fmt::LPaddedInt(res.GetNumber().GetNum());
    char ins_code=res.GetNumber().GetInsCode();
    if (ins_code!=0) {
      line_[26]=ins_code;
    }    
    ostr_ << line_;
    line_.Clear();
  }
  
  void SetIsPQR(bool t) {
    is_pqr_=t;
  }
private:
  std::ostream&       ostr_;
  int                 counter_;
  bool                is_pqr_;
  std::map<long,int>& atom_indices_;
  FormattedLine&      line_;
  mol::ResidueHandle  prev_;
  bool                peptide_;
  bool                charmm_style_;
};

class PDBConectWriterImpl : public mol::EntityVisitor {
public:
  PDBConectWriterImpl(std::ostream& ostream, std::map<long,int>& atom_indices)
    : ostr_(ostream), atom_indices_(atom_indices) {
  }
private:
public:
virtual bool VisitAtom(const mol::AtomHandle& atom) {
    if (atom.IsHetAtom()) {
      bool has_partner=false;
      int atom_index=atom_indices_[atom.GetHashCode()];
      mol::AtomHandleList partners=atom.GetBondPartners();
      std::list<int> partner_indices;
        for (auto partner : partners){
          mol::BondHandle bond = atom.FindBondToAtom(partner);
          int pind=atom_indices_[partner.GetHashCode()];
          if (pind!=0) {   
            for (int i=0; i < int(bond.GetBondOrder()); i++)
              {partner_indices.push_back(pind);}
            has_partner=true;
          }
        }
      if (has_partner) {
        write_conect(ostr_, atom_index, partner_indices);
      }
    }
    return true;
  }
private:
  std::ostream&       ostr_;
  std::map<long,int>& atom_indices_;
};

struct ForcePOSIX {
  std::locale old_locale;
  ForcePOSIX(){
    old_locale=std::locale();
    setlocale(LC_NUMERIC, "POSIX");
  }
  ~ForcePOSIX(){
    setlocale(LC_NUMERIC, old_locale.name().c_str());
  }
};

}

PDBWriter::PDBWriter(std::ostream& stream, const IOProfile& profile):
  outfile_(), outstream_(stream), mol_count_(0), line_(80),
  multi_model_(false), charmm_style_(profile.dialect=="CHARMM"), is_pqr_(false),
  profile_(profile)
{
  out_.push(outstream_);
}

PDBWriter::PDBWriter(const boost::filesystem::path& filename, 
                     const IOProfile& profile):
  outfile_(BFPathToString(filename).c_str()), outstream_(outfile_),
  mol_count_(0), line_(80), multi_model_(false), 
  charmm_style_(profile.dialect=="CHARMM"), is_pqr_(false),
  profile_(profile), filename_("")
{
  if (boost::iequals(".pqr", boost::filesystem::extension(filename))) {
    is_pqr_=true;
  }
  if (boost::iequals(".gz", boost::filesystem::extension(filename))) {
    out_.push(boost::iostreams::gzip_compressor());
  }
  out_.push(outstream_);
}

PDBWriter::PDBWriter(const String& filename, const IOProfile& profile):
  outfile_(filename.c_str()), outstream_(outfile_), mol_count_(0), line_(80), 
  multi_model_(false), charmm_style_(profile.dialect=="CHARMM"), 
  is_pqr_(false), profile_(profile), filename_(filename)
{
  if (!outfile_.is_open()) {
    throw IOException("Failed to open: " + filename);
  }
  if (boost::iequals(".pqr", boost::filesystem::extension(filename))) {
    is_pqr_=true;
  }
  if (boost::iequals(".gz", boost::filesystem::extension(filename))) {
    out_.push(boost::iostreams::gzip_compressor());
  }
  out_.push(outstream_);
}

void PDBWriter::WriteModelLeader()
{
  ++mol_count_;
  if (multi_model_) {
    out_ << "MODEL     " << mol_count_ << std::endl;
  } else if (mol_count_>1) {
    throw IOException("Trying to write several models into one file without "
                      "multi model mode enabled!");
  }
}

void PDBWriter::WriteModelTrailer()
{
  if (multi_model_) {
    out_ << "ENDMDL" << std::endl;
  }
}

template <typename H>
void PDBWriter::WriteModel(H ent)
{
  if (!out_) {
    if (!filename_.empty()) {
      std::stringstream ss;
      ss << "Can't write PDB to file '" << filename_ << "'";
      throw IOException(ss.str());
    }
    throw IOException("Can't write PDB. Bad stream");
  }
  ForcePOSIX posix;
  this->WriteModelLeader();
  PDBWriterImpl writer(out_, line_, atom_indices_, charmm_style_);
  writer.SetIsPQR(is_pqr_);
  ent.Apply(writer);
  PDBConectWriterImpl con_writer(out_, atom_indices_);
  ent.Apply(con_writer);
  this->WriteModelTrailer();
}

void PDBWriter::Write(const mol::EntityView& ent)
{
  this->WriteModel(ent);
}

void PDBWriter::Write(const mol::EntityHandle& ent)
{
  this->WriteModel(ent);
}

void PDBWriter::Write(const mol::AtomHandleList& atoms)
{
  ForcePOSIX posix = ForcePOSIX();
  this->WriteModelLeader();
  int counter=1;
  mol::ChainHandle last_chain;
  for (mol::AtomHandleList::const_iterator i=atoms.begin(),
       e=atoms.end(); i!=e; ++i, ++counter) {
    write_atom(out_, line_, *i, counter, is_pqr_, charmm_style_);
  }
  this->WriteModelTrailer();
}

  
PDBWriter::~PDBWriter()
{
  out_ << "END   " << std::endl;
}

}}
