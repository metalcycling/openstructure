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

#include "sdf_writer.hh"

#include <ost/boost_filesystem_helper.hh>
#include <ost/mol/atom_view.hh>
#include <ost/mol/residue_view.hh>
#include <ost/mol/chain_view.hh>
#include <ost/mol/bond_handle.hh>
#include <boost/regex.hpp>
#include <boost/bind/bind.hpp>

namespace ost { namespace io {

using boost::format;

namespace {

  class SDFAtomWriter : public mol::EntityViewVisitor {
    public:
      SDFAtomWriter(std::ostream& ostream, std::map<long, int>& atom_indices)
      : ostr_(ostream), atom_indices_(atom_indices), counter_(0) {
        atom_indices_.clear();
      }
    private:
    public:
      virtual bool VisitAtom(const mol::AtomView& atom) {
        atom_indices_[atom.GetHashCode()] = ++counter_;
        ostr_ << format("%10.4f") % atom.GetPos()[0]
              << format("%10.4f") % atom.GetPos()[1]
              << format("%10.4f ") % atom.GetPos()[2]
              << format("%-3s") % SDFAtomWriter::FormatEle(atom.GetElement())
              << " 0  0  0  0  0  0"
              << std::endl;
        return true;
      }

      static String FormatEle(const String& ele) {
        // OpenStructure has no strict requirements on lower or upper case
        // for elements. However, some sdf readers (read: OpenBabel) want the
        // first character to upper case, the rest in lower case.
        String return_ele = ele;
        if(!return_ele.empty()) return_ele[0] = toupper(return_ele[0]);
        for(size_t i = 1; i < return_ele.size(); ++i) {
          return_ele[i] = tolower(return_ele[i]);
        }
        return return_ele;
      }
    private:
      std::ostream&      ostr_;
      std::map<long, int>& atom_indices_;
      int counter_;
  };

  class SDFBondWriter : public mol::EntityViewVisitor {
  public:
    SDFBondWriter(std::ostream& ostream,
                  const std::map<long, int>& atom_indices)
      : ostr_(ostream), atom_indices_(atom_indices), counter_(0) {
    }
  private:
    // compare two atoms according to their indices (used for sorting)
    bool CompareAtomIdx(const mol::AtomView& first,
                        const mol::AtomView& second) {
      std::map<long, int>::const_iterator aidx_first(
                                  atom_indices_.find(first.GetHashCode()));
      std::map<long, int>::const_iterator aidx_second(
                                  atom_indices_.find(second.GetHashCode()));

      if(aidx_first==atom_indices_.end() || aidx_second==atom_indices_.end()) {
        throw IOException("Cannot write bond: atom idx not found for sorting");
      }
      return (aidx_first->second < aidx_second->second);
    }

  public:
    virtual bool VisitAtom(const mol::AtomView& atom) {
      ++counter_; // current atom index

      // get all neighboring atoms and sort them according to their atom index
      mol::AtomViewList atoms = atom.GetBondPartners();
      std::sort(atoms.begin(), atoms.end(), bind(&SDFBondWriter::CompareAtomIdx,
                                                 this, boost::placeholders::_1,
                                                 boost::placeholders::_2));

      // iterate all neighboring atoms and print bonds to all atoms with index
      // larger than current atom index
      for(mol::AtomViewList::iterator atom_iter = atoms.begin();
          atom_iter != atoms.end(); ++atom_iter) {
        std::map<long, int>::const_iterator aidx(
                               atom_indices_.find((*atom_iter).GetHashCode()));

        // check if index was found
        if(aidx==atom_indices_.end()) {
          throw IOException("Cannot write bond between " +
                            atom.GetQualifiedName() + " and " +
                            atom_iter->GetQualifiedName() +
                            ": atom index not found");
        }

        // only print bonds to atoms with larger index than current index
        if(aidx->second > counter_) {
          mol::BondHandle bond(atom.GetHandle().FindBondToAtom(
                                                   atom_iter->GetHandle()));
          if(!bond.IsValid()) {
            throw IOException("Bond is invalid between " +
                              atom.GetQualifiedName() + " and " +
                              atom_iter->GetQualifiedName());
          }
          int type = bond.GetBondOrder();
          ostr_ << format("%3i") % counter_
                << format("%3i") % aidx->second
                << format("%3i") % type
                << "  0  0  0"
                << std::endl;
        }
      }
      return true;
    }

  private:
    std::ostream&      ostr_;
    const std::map<long, int>& atom_indices_;
    int counter_;
  };
}

SDFWriter::SDFWriter(std::ostream& ostream)
  : outfile_(), ostr_(ostream), counter_(0), atom_indices_() {
}

SDFWriter::SDFWriter(const String& filename)
  : outfile_(filename.c_str()), ostr_(outfile_), counter_(0), atom_indices_() {
}

SDFWriter::SDFWriter(const boost::filesystem::path& filename): 
  outfile_(BFPathToString(filename).c_str()),
  ostr_(outfile_), counter_(0), atom_indices_() {}

void SDFWriter::Write(const mol::EntityView& ent) {
  if (!ostr_) {
    throw IOException("Can't write SDF file. Bad output stream");
  }
  mol::EntityView non_const_view = ent;
  non_const_view.Apply(*this);
}

void SDFWriter::Write(const mol::EntityHandle& ent) {
  if (!ostr_) {
    throw IOException("Can't write SDF file. Bad output stream");
  }
  mol::EntityView non_const_view = ent.CreateFullView();
  non_const_view.Apply(*this);
}

bool SDFWriter::VisitChain(const mol::ChainView& chain) {
  // print end of molecule line
  if(counter_ != 0) {
    ostr_ << "$$$$" << std::endl;
    counter_ = 0;
    atom_indices_.clear();
  }

  // remove chain number if added when reading from sdf file
  String cname = chain.GetName();
  if (cname.length()>6) {
    boost::regex pattern = boost::regex("^[0-9]{5}_");
    if (boost::regex_search(cname, pattern)) {
      cname = cname.substr(6);
    }
  }

  // print header lines
  ostr_ << cname << std::endl;
  ostr_ << std::endl;
  ostr_ << std::endl;
  // print counts line
  ostr_ << format("%3d") % chain.GetAtomCount()
        << format("%3d") % chain.GetBondCount()
        << "  0  0  0  0            999 V2000"
        << std::endl;

  // write atom block
  SDFAtomWriter atom_writer(ostr_, atom_indices_);
  mol::ChainView non_const_chain = chain;
  non_const_chain.Apply(atom_writer);

  // write bond block
  SDFBondWriter bond_writer(ostr_, atom_indices_);
  non_const_chain.Apply(bond_writer);

  // write property block
  //TODO: write property block
  ostr_ << "M  END" << std::endl;

  // write data block
  std::map<String,GenericPropValue> prop_map = non_const_chain.GetPropMap();
  std::map<String,GenericPropValue>::iterator iter;
  for(iter = prop_map.begin(); iter != prop_map.end(); ++iter) {
    ostr_ << "> <" << (*iter).first << ">" << std::endl;
    ostr_ << (*iter).second << std::endl;
    ostr_ << std::endl;
  }

  // write molecule endline
  ostr_ << "$$$$" << std::endl;

  return true;
}

}}
