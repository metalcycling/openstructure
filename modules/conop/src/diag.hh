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
#ifndef OST_CONOP_DIAG_HH
#define OST_CONOP_DIAG_HH

#include <ost/mol/atom_handle.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_handle.hh>
#include <ost/conop/module_config.hh>

namespace ost { namespace conop {

typedef enum {
  DIAG_ARG_TYPE_ATOM,
  DIAG_ARG_TYPE_RESIDUE,
  DIAG_ARG_TYPE_CHAIN,
  DIAG_ARG_TYPE_STRING,
  DIAG_ARG_TYPE_INT
} DiagArgType;

typedef enum {
  DIAG_UNK_ATOM,
  DIAG_UNK_RESIDUE,
  DIAG_MISSING_ATOM,
  DIAG_NONSTD_RESIDUE
} DiagType;

class DLLEXPORT_OST_CONOP Diag {
public:
  Diag(DiagType typ, const char* fmt): type_(typ), format_(fmt) {}
  DiagType GetType() const { return type_; }
  Diag& AddAtom(mol::AtomHandle atom)
  {
    atoms_.push_back(atom);
    args_.push_back(ArgDesc(atoms_.size()-1, DIAG_ARG_TYPE_ATOM));
    return *this;
  }

  Diag& AddResidue(mol::ResidueHandle res)
  {
    residues_.push_back(res);
    args_.push_back(ArgDesc(residues_.size()-1, DIAG_ARG_TYPE_RESIDUE));
    return *this;
  }
  Diag& AddChain(mol::ChainHandle chain)
  {
    chains_.push_back(chain);
    args_.push_back(ArgDesc(chains_.size()-1, DIAG_ARG_TYPE_CHAIN));
    return *this;
  }
  Diag& AddInt(int int_val)
  {
    ints_.push_back(int_val);
    args_.push_back(ArgDesc(ints_.size()-1, DIAG_ARG_TYPE_INT));
    return *this;
  }
  Diag& AddString(const String& str)
  {
    strings_.push_back(str);
    args_.push_back(ArgDesc(strings_.size()-1, DIAG_ARG_TYPE_STRING));
    return *this;
  }
  mol::AtomHandle GetAtom(size_t index) const
  {
    assert(index<args_.size());
    return atoms_[args_[index].index];
  }
  mol::ResidueHandle GetResidue(size_t index) const
  {
    assert(index<args_.size());
    return residues_[args_[index].index];
  }
  mol::ChainHandle GetChain(size_t index) const
  {
    assert(index<args_.size());
    return chains_[args_[index].index];
  }
  String Format(bool colored=true) const;
private:
  struct ArgDesc {
    ArgDesc(size_t i, DiagArgType t): index(i), type(t) { }
    size_t      index;
    DiagArgType type;
  };
  DiagType               type_;
  String                 format_;
  mol::AtomHandleList    atoms_;
  mol::ResidueHandleList residues_;
  mol::ChainHandleList   chains_;
  std::vector<String>    strings_;
  std::vector<int>       ints_;
  std::vector<ArgDesc>   args_;
};

class Diagnostics;

typedef boost::shared_ptr<Diagnostics> DiagnosticsPtr;


class DLLEXPORT DiagError : public Error {
public:
  DiagError(const Diag& diag) : Error(diag.Format(false)) {}
};

class DLLEXPORT_OST_CONOP Diagnostics {
public:
  typedef std::vector<Diag*>::iterator diag_iterator;
  typedef std::vector<Diag*>::const_iterator const_diag_iterator;
  Diagnostics() {}

  ~Diagnostics()
  {
    for(std::vector<Diag*>::iterator
        i=diags_.begin(), e=diags_.end(); i!=e;++i) {
      delete *i;
    }
  }

  Diag& AddDiag(DiagType type, const char* fmt)
  {
    diags_.push_back(new Diag(type, fmt));
    return *diags_.back();
  }
  diag_iterator diags_begin() { return diags_.begin(); } 
  diag_iterator diags_end() { return diags_.end(); } 
  const_diag_iterator diags_begin() const { return diags_.begin(); } 
  const_diag_iterator diags_end() const { return diags_.end(); } 
  size_t diag_count() const { return diags_.size(); }
private:
  std::vector<Diag*> diags_;
};

}} /* ost::conop */
#endif
