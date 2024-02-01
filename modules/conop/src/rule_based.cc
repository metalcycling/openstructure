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
#include <limits>
#include <ost/log.hh>
#include <ost/profile.hh>
#include <ost/message.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/mol/torsion_handle.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/residue_handle.hh>
#include "rule_based.hh"
#include "conop.hh"

namespace ost { namespace conop {

void RuleBasedProcessor::DoProcess(DiagnosticsPtr diags, 
                                   mol::EntityHandle ent) const
{
  Profile prof("RuleBasedProcessor::Process");
  mol::ChainHandleList chains=ent.GetChainList();
  for (mol::ChainHandleList::iterator 
       i = chains.begin(), e = chains.end(); i!= e; ++i) {
    mol::ChainHandle& chain = *i;
    mol::ResidueHandleList residues = chain.GetResidueList();
    mol::ResidueHandle prev;
    for (mol::ResidueHandleList::iterator 
         j = residues.begin(), e2 = residues.end(); j != e2; ++j) {
      mol::ResidueHandle residue = *j;
      mol::AtomHandleList atoms_to_connect;
      CompoundPtr compound = lib_->FindCompound(residue.GetName(), Compound::PDB);
      if (!compound && this->GetConnect()) {
        // process unknown residue...
        this->ProcessUnkResidue(diags, residue, atoms_to_connect);
        for (mol::AtomHandleList::iterator k = atoms_to_connect.begin(),
             e3=atoms_to_connect.end(); k!= e3; ++k) {
          this->DistanceBasedConnect(*k);
        }
        continue;
      }
      this->ReorderAtoms(residue, compound, this->GetFixElement());  
      bool unks = this->HasUnknownAtoms(residue, compound, 
                                        this->GetStrictHydrogens());
      if (unks) {
        mol::AtomHandleList unk_atoms;
        unk_atoms = GetUnknownAtomsOfResidue(residue, compound, 
                                             this->GetStrictHydrogens());
        this->ProcessUnkAtoms(diags, residue, unk_atoms, atoms_to_connect);
        residue.SetChemClass(mol::ChemClass(mol::ChemClass::UNKNOWN));
        residue.SetChemType(mol::ChemType(mol::ChemType::UNKNOWN));
        residue.SetOneLetterCode('?');
      } else {
        this->FillResidueProps(residue, compound);
      }
      if (this->GetConnect()) {
        this->ConnectAtomsOfResidue(residue, compound, 
                                    this->GetStrictHydrogens());
        if (this->GetConnectAminoAcids())
          this->ConnectResidues(prev, residue);
        for (mol::AtomHandleList::iterator k = atoms_to_connect.begin(),
             e3=atoms_to_connect.end(); k!= e3; ++k) {
          this->DistanceBasedConnect(*k);
        }
        if (!this->GetStrictHydrogens()) {
          mol::AtomHandleList atoms = residue.GetAtomList();
          for (mol::AtomHandleList::iterator k = atoms.begin(),
               e3 = atoms.end(); k != e3; ++k) {
            const String& ele = k->GetElement();
            if ((ele == "D" || ele == "H") && k->GetBondCount() == 0) {
              this->DistanceBasedConnect(*k);
            }
          }
        }
      }
      prev = residue;
    }
    if (residues.empty() || !this->GetAssignTorsions()) {
      continue;
    }
    AssignBackboneTorsions(residues);
  }
}

void RuleBasedProcessor::ProcessUnkResidue(DiagnosticsPtr diags,
                                           mol::ResidueHandle res,
                                           mol::AtomHandleList& remaining_atoms) const
{
  mol::AtomHandleList atoms = res.GetAtomList();    
  if (this->GetUnkResidueTreatment() == CONOP_SILENT ||
      this->GetUnkResidueTreatment() == CONOP_WARN) {
    res.SetOneLetterCode('?');
    res.SetChemClass(mol::ChemClass(mol::ChemClass::UNKNOWN));
    for (mol::AtomHandleList::iterator i = atoms.begin(), 
         e = atoms.end(); i !=e; ++i) {
      remaining_atoms.push_back(*i);
      if (!Conopology::Instance().IsValidElement(i->GetElement()))
        i->SetElement(GuessAtomElement(i->GetName(), i->IsHetAtom(), atoms.size()));
    }
  }
  // Don't do anything if treatment is set to SILENT
  if (this->GetUnkResidueTreatment() == CONOP_SILENT) {
    return;
  }
  switch (this->GetUnkResidueTreatment()) {
    case CONOP_WARN:
      diags->AddDiag(DIAG_UNK_RESIDUE, "unknown residue %0")
            .AddResidue(res);
      break;
    case CONOP_FATAL:
      throw DiagError(Diag(DIAG_UNK_RESIDUE, "unknown residue %0").AddResidue(res));
      break;
    case CONOP_REMOVE_RESIDUE:
    case CONOP_REMOVE:
      diags->AddDiag(DIAG_UNK_RESIDUE, "removed unknown residue %0")
           .AddResidue(res);
      res.GetEntity().EditXCS().DeleteResidue(res);
      return;
    default:
      assert(0 && "unhandled switch");
  }
}

void RuleBasedProcessor::ProcessUnkAtoms(DiagnosticsPtr diags,
                                         mol::ResidueHandle res,
                                         mol::AtomHandleList unks,
                                         mol::AtomHandleList& remaining_atoms) const
{
  if (this->GetUnkAtomTreatment() == CONOP_SILENT ||
      this->GetUnkAtomTreatment() == CONOP_WARN) {
    for (mol::AtomHandleList::iterator i = unks.begin(), 
         e = unks.end(); i !=e; ++i) {
      remaining_atoms.push_back(*i);
      if (!Conopology::Instance().IsValidElement(i->GetElement()))
        i->SetElement(GuessAtomElement(i->GetName(), i->IsHetAtom(), res.GetAtomCount()));
    }
  }

  // Don't do anything if treatment is set to SILENT
  if (this->GetUnkAtomTreatment() == CONOP_SILENT) {
    return;
  }
 
  assert(!unks.empty() && "empty unk list");
  mol::XCSEditor edi = unks.front().GetEntity().EditXCS();
  for (mol::AtomHandleList::iterator 
       i = unks.begin(), e = unks.end(); i != e; ++i) {
    switch (this->GetUnkAtomTreatment()) {
      case CONOP_REMOVE_ATOM:
      case CONOP_REMOVE:
        edi.DeleteAtom(*i);
        diags->AddDiag(DIAG_UNK_ATOM, "removed unknown atom %0 from residue %1")
             .AddString(i->GetName()).AddResidue(i->GetResidue());
        break;
      case CONOP_WARN:
        diags->AddDiag(DIAG_UNK_ATOM, "residue %0 contains unknown atom %1")
             .AddResidue(i->GetResidue()).AddString(i->GetName());
        break;
      case CONOP_FATAL:
        throw DiagError(Diag(DIAG_UNK_ATOM, "unknown atom %0").AddAtom(*i));
        break;
      case CONOP_REMOVE_RESIDUE:
        diags->AddDiag(DIAG_UNK_ATOM, "removed residue %0 with unknown atom %1")
             .AddString(i->GetResidue().GetQualifiedName())
             .AddString(i->GetName());
        edi.DeleteResidue(i->GetResidue());
        return;
      default:
        assert(0 && "unhandled switch");
    }
  }
}

String RuleBasedProcessor::ToString() const {
  std::stringstream ss;
  ss << "RuleBasedProcessor(" << this->OptionsToString() 
     << ", fix_element=" << (fix_element_ ? "True" : "False")
     << ", strict_hydrogens=" << (strict_hydrogens_ ? "True" : "False")
     << ", unk_res_treatment='" << StringFromConopAction(unk_res_treatment_)
     << "', unk_atom_treatment='" << StringFromConopAction(unk_atom_treatment_)
     << "')";
  return ss.str();
}

void RuleBasedProcessor::_CheckLib() const {
  if (!lib_) {
    throw Error("Cannot initialize RuleBasedProcessor without a valid "
                "CompoundLib object!");
  }
}

}}
