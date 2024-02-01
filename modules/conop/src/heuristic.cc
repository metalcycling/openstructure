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

#include <ost/log.hh>
#include <ost/profile.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/mol/torsion_handle.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/residue_handle.hh>
#include "heuristic.hh"
#include "conop.hh"

namespace ost { namespace conop {

void HeuristicProcessor::ProcessUnkResidue(DiagnosticsPtr diags, 
                                           mol::ResidueHandle res) const {
  mol::AtomHandleList atoms = res.GetAtomList(); 
  for (mol::AtomHandleList::iterator 
       i = atoms.begin(), e = atoms.end(); i != e; ++i) {
    mol::AtomHandle a = *i;
    if (!Conopology::Instance().IsValidElement(a.GetElement()))
      a.SetElement(GuessAtomElement(a.GetName(), a.IsHetAtom(), atoms.size()));
  }

  res.SetChemClass(GuessChemClass(res));

  if (!this->GetConnect()) { return; }
  for (mol::AtomHandleList::iterator 
       i = atoms.begin(), e = atoms.end(); i != e; ++i) {
    this->DistanceBasedConnect(*i);
  }
}

void HeuristicProcessor::DoProcess(DiagnosticsPtr diags, 
                                   mol::EntityHandle ent) const
{
  Profile prof("HeuristicProcessor::Process");
  mol::ChainHandleList chains=ent.GetChainList();
  for (mol::ChainHandleList::iterator 
       i = chains.begin(), e = chains.end(); i!= e; ++i) {
    mol::ChainHandle& chain = *i;
    mol::ResidueHandleList residues = chain.GetResidueList();
    mol::ResidueHandle prev;
    for (mol::ResidueHandleList::iterator 
         j = residues.begin(), e2 = residues.end(); j != e2; ++j) {
      mol::ResidueHandle residue = *j;
      CompoundPtr compound = lib_.FindCompound(residue.GetName(), Compound::PDB);
      if (!compound) {
        // process unknown residue...
        this->ProcessUnkResidue(diags, residue);
        if (this->GetConnectAminoAcids()) 
          this->ConnectResidues(prev, residue);
        prev = residue;
        continue;
      }
      this->ReorderAtoms(residue, compound, true);  
      this->FillResidueProps(residue, compound);
      if (this->GetConnect()) {
        this->ConnectAtomsOfResidue(residue, compound, false);
        if (this->GetConnectAminoAcids()) 
          this->ConnectResidues(prev, residue);
      }
      prev = residue;
      if (!this->GetConnect()) { continue; }
      mol::AtomHandleList atoms = residue.GetAtomList();
      for (mol::AtomHandleList::iterator 
          i = atoms.begin(), e = atoms.end(); i != e; ++i) {
        if (i->GetBondCount() == 0)
          this->DistanceBasedConnect(*i);
      }
    }
    if (residues.empty() || !this->GetAssignTorsions()) {
      continue;
    }
    AssignBackboneTorsions(residues);
  }
}

String HeuristicProcessor::ToString() const {
  std::stringstream ss;
  ss << "HeuristicProcessor(" << this->OptionsToString() << ")";
  return ss.str();
}

}}
