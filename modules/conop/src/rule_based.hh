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
#ifndef OST_CONOP_RULE_BASED_HH
#define OST_CONOP_RULE_BASED_HH

#include <ost/mol/entity_handle.hh>
#include "compound_lib.hh"
#include "diag.hh"
#include "processor.hh"

namespace ost { namespace conop {

/// \brief returns all atoms not listed in  the specifictaion of compound
mol::AtomHandleList DLLEXPORT_OST_CONOP GetUnknownAtoms(mol::ResidueHandle res, 
                                                        CompoundPtr compound);

class RuleBasedProcessor;

typedef boost::shared_ptr<RuleBasedProcessor> RuleBasedProcessorPtr;

class DLLEXPORT_OST_CONOP RuleBasedProcessor  : public Processor {
public:
  RuleBasedProcessor(CompoundLibPtr compound_lib): 
    lib_(compound_lib), fix_element_(true), strict_hydrogens_(false), 
    unk_res_treatment_(CONOP_WARN), unk_atom_treatment_(CONOP_WARN)
  {
    _CheckLib();
  }

  RuleBasedProcessor(CompoundLibPtr compound_lib, bool fe, bool sh,
                     ConopAction ur, ConopAction ua, bool bf, bool at, bool cn,
                     bool aa, bool ch, ConopAction zo): 
    Processor(bf, at, cn, aa, ch, zo), lib_(compound_lib), fix_element_(fe), 
    strict_hydrogens_(sh), unk_res_treatment_(ur), 
    unk_atom_treatment_(ua)
  {
    _CheckLib();
  }

  ConopAction GetUnkResidueTreatment() const {
    return unk_res_treatment_;
  }

  ConopAction GetUnkAtomTreatment() const {
    return unk_atom_treatment_;
  }

  bool GetFixElement() const {
    return fix_element_;
  }
  void SetFixElement(bool flag) {
    fix_element_ = flag;
  }
  bool GetStrictHydrogens() const {
    return strict_hydrogens_;
  }

  void SetStrictHydrogens(bool flag) {
    strict_hydrogens_ = flag;
  }
  void SetUnkResidueTreatment(ConopAction action) {
    unk_res_treatment_ = action;
  }

  void SetUnkAtomTreatment(ConopAction action) {
    unk_atom_treatment_ = action;
  }
  virtual ProcessorPtr Copy() const {
    return ProcessorPtr(new RuleBasedProcessor(*this));
  }

  virtual String ToString() const;
protected:
  void ProcessUnkResidue(DiagnosticsPtr diags,
                         mol::ResidueHandle res, 
                         mol::AtomHandleList& remaining) const;
  void ProcessUnkAtoms(DiagnosticsPtr diags,
                       mol::ResidueHandle res,
                       mol::AtomHandleList unks,
                       mol::AtomHandleList& remaining) const;
  virtual void DoProcess(DiagnosticsPtr diags, 
                         mol::EntityHandle ent) const;
private:
  void _CheckLib() const;

  CompoundLibPtr lib_;
  bool fix_element_;
  bool strict_hydrogens_;
  ConopAction unk_res_treatment_;
  ConopAction unk_atom_treatment_;
};



}}
#endif

