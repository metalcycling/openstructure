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
#include <ost/mol/xcs_editor.hh>
#include <ost/mol/bond_handle.hh>
#include <ost/mol/torsion_handle.hh>
#include <ost/mol/impl/residue_impl.hh>
#include <ost/mol/impl/atom_impl.hh>
#include <ost/mol/residue_handle.hh>
#include "processor.hh"

namespace ost { namespace conop {

struct OrdinalAtomComp {
  bool operator()(const mol::impl::AtomImplPtr& a,
                  const mol::impl::AtomImplPtr& b) const {
    return a->GetState()<b->GetState();
  }
};
DiagnosticsPtr Processor::Process(mol::EntityHandle ent, 
                                  bool log_diags) const
{
  DiagnosticsPtr diags(new Diagnostics);
  if (!this->BeginProcessing(diags, ent)) {
    return diags;
  }
  this->DoProcess(diags, ent);
  
  this->EndProcessing(diags, ent);
  if (log_diags) {
    for (Diagnostics::diag_iterator i = diags->diags_begin(),
         e = diags->diags_end(); i != e; ++i) {
      LOG_WARNING((*i)->Format(false));
    }
  }
  return diags;
}

void AssignBackboneTorsions(mol::ChainHandle chain)
{
  AssignBackboneTorsions(chain.GetResidueList());
}

void AssignBackboneTorsions(mol::ResidueHandleList residues)
{
  if (residues.empty()) { return; }
  mol::ResidueHandle r1;
  mol::ResidueHandle r2;
  mol::ResidueHandle r3 = residues.front();
  for (mol::ResidueHandleList::iterator 
        j = residues.begin()+1, e2 = residues.end(); j != e2; ++j) {
    r1 = r2;
    r2 = r3;
    r3 = *j;
    AssignBackboneTorsions(r1, r2, r3); 
  }
  AssignBackboneTorsions(r2, r3, mol::ResidueHandle());
}


void AssignBackboneTorsions(mol::ResidueHandle prev, 
                            mol::ResidueHandle res, 
                            mol::ResidueHandle next)
{

  mol::XCSEditor e=res.GetEntity().EditXCS(mol::BUFFERED_EDIT);
  //psi
  if (next.IsValid() && next.IsPeptideLinking()){
    mol::AtomHandle ca_this=res.FindAtom("CA");
    mol::AtomHandle n_this=res.FindAtom("N");
    mol::AtomHandle c_this=res.FindAtom("C");
    mol::AtomHandle n_next=next.FindAtom("N");
    if ((ca_this && n_this && c_this && n_next &&  BondExists(c_this, n_next))
        && !res.GetPsiTorsion()) {
      e.AddTorsion("PSI", n_this, ca_this, c_this, n_next);
    }
  };
  //phi
  if (prev.IsValid() && prev.IsPeptideLinking()) {
    mol::AtomHandle c_prev=prev.FindAtom("C");
    mol::AtomHandle n_this=res.FindAtom("N");
    mol::AtomHandle ca_this=res.FindAtom("CA");
    mol::AtomHandle c_this=res.FindAtom("C");
    if ((c_prev && n_this && ca_this && c_this && BondExists(c_prev, n_this))
        && !res.GetPhiTorsion()) {
      e.AddTorsion("PHI", c_prev, n_this, ca_this, c_this);
    }
  }
  //omega
  if (prev.IsValid() && prev.IsPeptideLinking()) {
    mol::AtomHandle ca_prev=prev.FindAtom("CA");
    mol::AtomHandle c_prev=prev.FindAtom("C");
    mol::AtomHandle n=res.FindAtom("N");
    mol::AtomHandle ca=res.FindAtom("CA");
    if ((ca_prev && c_prev && n && ca && BondExists(c_prev, n))
        && !res.GetOmegaTorsion()) {
      e.AddTorsion("OMEGA",ca_prev , c_prev, n, ca);
    }
  }
} 

mol::ChemClass GuessChemClass(mol::ResidueHandle res)
{
  mol::AtomHandle ca=res.FindAtom("CA");
  if (!ca.IsValid() || ca.GetElement()!="C") return mol::ChemClass();
  mol::AtomHandle n=res.FindAtom("N");
  if (!n.IsValid() || n.GetElement()!="N") return mol::ChemClass();
  mol::AtomHandle c=res.FindAtom("C");
  if (!c.IsValid() || c.GetElement()!="C") return mol::ChemClass();
  mol::AtomHandle o=res.FindAtom("O");
  if (!o.IsValid() || o.GetElement()!="O") return mol::ChemClass();
  if (IsBondFeasible(n, ca) && IsBondFeasible(ca, c) &&
      IsBondFeasible(c, o)) {
    return mol::ChemClass(mol::ChemClass::PEPTIDE_LINKING);
  }
  return mol::ChemClass();
}

String GuessAtomElement(const String& aname, bool hetatm, int atom_count)
{
  static String l1[] = {
    "H","C","N","O","P","S","K"
  };
  static int l1c=7;
  static String l2[] = {
    "NA","MG","AL","SI","CL",
    "CA","CR","MN","FE","CO","NI","CU","ZN","AS","SE","BR","MO"
  };
  static int l2c=17;
  static String l3[] = {
    "B","F","V"
  };
  static int l3c=3;
  static String l4[] = {
    "HE","LI","BE","AR","SC","TI","GA","GE","KR"
  };
  static int l4c=9;

  String ele=aname.substr(0,2);
  // hydrogen hack
  if(ele[0]=='H') {
    return "H";
  }
  if (!hetatm || atom_count > 1) {
    if (ele=="CA" || ele=="CB") {
      return "C";
    }
  }

  // two characters
  if(aname.size()==2) {
    for(int i=0;i<l2c;i++) {
      if(ele==l2[i]) return ele;
    }
    // check second character for match
    for(int i=0;i<l1c;i++) {
      if(ele[1]==l1[i][0]) {
        return l1[i];
      }
    }
    // still no match, repeat with less likely tables
    for(int i=0;i<l4c;i++) {
      if(ele==l4[i]) return ele;
    }
    // check second character for match
    for(int i=0;i<l3c;i++) {
      if(ele[1]==l3[i][0]) {
        return l3[i];
      }
    }
    // check second character for match
    for(int i=0;i<l1c;i++) {
      if(ele[0]==l1[i][0]) {
        return l1[i];
      }
    }
  } else {
    for(int i=0;i<l1c;i++) {
      if(ele==l1[i]) return ele;
    }
    for(int i=0;i<l3c;i++) {
      if(ele==l3[i]) return ele;
    }
  }
  size_t i=0;
  while (i<aname.size() && isdigit(aname[i])) {
    ++i;
  }
  return i<aname.size() ? String(1, aname[i]) : "";
}

bool IsBondFeasible(const mol::AtomHandle& atom_a,
                    const mol::AtomHandle& atom_b)
{
  Real radii=0.0;
  if (atom_a.GetRadius()>0.0) {
    radii=atom_a.GetRadius();
  } else {
    return false;
  }
  if (atom_b.GetRadius()>0.0) {
    radii+=atom_b.GetRadius();
  } else {
    return false;
  } 
  Real len=geom::Length2(atom_a.GetPos()-atom_b.GetPos());
  Real lower_bound=radii*radii*0.0625;
  Real upper_bound=lower_bound*6.0;
  return (len<=upper_bound && len>=lower_bound);
}

void Processor::ConnectResidues(mol::ResidueHandle rh,
                                mol::ResidueHandle next) const
{
  if (!next.IsValid() || !rh.IsValid()) {
    return;
  }

  mol::XCSEditor e=rh.GetEntity().EditXCS(mol::BUFFERED_EDIT);

  // check if both of the residues are able to form a peptide bond.
  if (rh.IsPeptideLinking() && next.IsPeptideLinking()) {
    // If we have an OXT then there is no peptide bond connecting the two
    // residues.
    if (rh.FindAtom("OXT"))
      return;
    mol::AtomHandle c=rh.FindAtom("C");
    mol::AtomHandle n=next.FindAtom("N");
    // Give subclasses a chance to give us their opinions on the feasibility of
    // the peptide bond.
    if (c.IsValid() && n.IsValid() && IsBondFeasible(c, n)) {
      e.Connect(c, n, 1);
      rh.SetIsProtein(true);
      next.SetIsProtein(true);
    }
  } else if (rh.IsNucleotideLinking() && next.IsNucleotideLinking()) {
    mol::AtomHandle c=rh.FindAtom("O3'");
    mol::AtomHandle n=next.FindAtom("P");
    if (c.IsValid() && n.IsValid() && IsBondFeasible(c, n)) {
      e.Connect(c, n, 1);
    }
  }
}


mol::AtomHandle Processor::LocateAtom(const mol::AtomHandleList& ahl, 
                                      int ordinal) const
{
  if (ahl.empty())
    return mol::AtomHandle();
  const mol::AtomHandle* r_it=&ahl.back();
  if (static_cast<int>(ahl.size())>ordinal) {
    r_it=&ahl.front()+ordinal;
  }
  while ((r_it>=&ahl.front()) &&
         (static_cast<int>(r_it->Impl()->GetState())>ordinal)) {
    --r_it;
  }
  bool not_found=(r_it<&ahl.front() ||
                 static_cast<int>(r_it->Impl()->GetState())!=ordinal);
  return  not_found ? mol::AtomHandle() : *r_it;
}

void Processor::FillResidueProps(mol::ResidueHandle residue,
                                 CompoundPtr compound) const
{
  residue.SetChemClass(compound->GetChemClass());
  residue.SetChemType(compound->GetChemType());
  residue.SetOneLetterCode(compound->GetOneLetterCode());
}

void Processor::ConnectAtomsOfResidue(mol::ResidueHandle rh, 
                                      CompoundPtr compound,
                                      bool strict_hydrogens) const
{
  mol::XCSEditor e=rh.GetEntity().EditXCS(mol::BUFFERED_EDIT);
  BondSpecList::const_iterator j=compound->GetBondSpecs().begin();
  mol::AtomHandleList atoms=rh.GetAtomList();
  for(; j!=compound->GetBondSpecs().end(); ++j) {
      const BondSpec& bond=*j;
      mol::AtomHandle a1=this->LocateAtom(atoms, bond.atom_one);
      mol::AtomHandle a2=this->LocateAtom(atoms, bond.atom_two);
      if (a1.IsValid() && a2.IsValid()) { 
        if (!strict_hydrogens && (a1.GetElement()=="H" ||
                                  a1.GetElement()=="D" ||
                                  a2.GetElement()=="H" ||
                                  a2.GetElement()=="D")) {
          continue;
        }
        if (!connect_hetatm_ && a1.IsHetAtom() && a2.IsHetAtom()) {
          continue;
        }
        if (!this->GetCheckBondFeasibility()) {
          e.Connect(a1, a2, bond.order);
        } else { 
          if (IsBondFeasible(a1, a2)) {
            e.Connect(a1, a2, bond.order);
          }
        }
      }
  }
}


void Processor::ReorderAtoms(mol::ResidueHandle residue, 
                             CompoundPtr compound, 
                             bool fix_element) const
{
  mol::impl::ResidueImplPtr impl=residue.Impl();
  mol::impl::AtomImplList::iterator i=impl->GetAtomList().begin();
  for (; i!=impl->GetAtomList().end(); ++i) {
    mol::impl::AtomImplPtr atom=*i;
    atom->SetState(std::numeric_limits<int>::max());
    int index=compound->GetAtomSpecIndex(atom->GetName());
    if (index==-1) {
      atom->SetState(std::numeric_limits<int>::max());
      continue;
    }
    atom->SetState((compound->GetAtomSpecs())[index].ordinal);
    // override element
    if (fix_element) {
      atom->SetElement((compound->GetAtomSpecs())[index].element);
    }
  }
  std::sort(impl->GetAtomList().begin(), impl->GetAtomList().end(),
            OrdinalAtomComp());
}

bool Processor::HasUnknownAtoms(mol::ResidueHandle res, 
                                CompoundPtr compound,
                                bool strict_hydrogens) const 
{
  mol::AtomHandleList atoms=res.GetAtomList();
  for (mol::AtomHandleList::iterator 
       i=atoms.begin(), e=atoms.end(); i!=e; ++i) {
    if ((*i).Impl()->GetState()==std::numeric_limits<int>::max()) {
      if (((*i).GetElement()=="H" || (*i).GetElement()=="D") && 
          !strict_hydrogens) {
        continue;
      }
      return true;
    }
  }
  return false;
}

mol::AtomHandleList GetUnknownAtomsOfResidue(mol::ResidueHandle res, 
                                             CompoundPtr compound,
                                             bool strict_hydrogens)
{
  mol::AtomHandleList atoms=res.GetAtomList();
  mol::AtomHandleList unks;
  const AtomSpecList& atom_specs=compound->GetAtomSpecs();
  for (mol::AtomHandleList::const_iterator
        j=atoms.begin(), e2=atoms.end(); j!=e2; ++j) {
    bool found=false;
    for (AtomSpecList::const_iterator
          k=atom_specs.begin(), e3=atom_specs.end(); k!=e3; ++k) {
      if (k->name==j->GetName() || k->alt_name==j->GetName()) {
        found=true;
        break;
      }
    }
    if (!found) {
      unks.push_back(*j);
    }
  }
  return unks;
}

void Processor::DistanceBasedConnect(mol::AtomHandle atom) const
{
  mol::EntityHandle ent=atom.GetEntity();
  mol::XCSEditor editor=ent.EditXCS(mol::BUFFERED_EDIT);
  mol::AtomHandleList alist = ent.FindWithin(atom.GetPos(),4.0);
  mol::ResidueHandle res_a=atom.GetResidue();
  for (mol::AtomHandleList::const_iterator it=alist.begin(),
       e=alist.end();it!=e;++it) {
    if (*it!=atom) {
      if (!connect_hetatm_ && it->IsHetAtom() && atom.IsHetAtom()) {
        continue;
      }
      if (IsBondFeasible(atom, *it)) {
        if (it->GetResidue()==res_a ||
            Processor::AreResiduesConsecutive(res_a, it->GetResidue())) {
            editor.Connect(*it, atom);
        }
      }
    }
  }
}

String StringFromConopAction(ConopAction action) {
  switch (action) {
    case CONOP_FATAL:
      return "fatal";
    case CONOP_SILENT:
      return "silent";
    case CONOP_REMOVE:
      return "rm";
    case CONOP_REMOVE_RESIDUE:
      return "rm_residue";
    case CONOP_REMOVE_ATOM:
      return "rm_atom";
    case CONOP_WARN:
      return "warn";
  }
  assert(0 && "unhandled conop action");
  return "unhandled conop action";
}

bool Processor::AreResiduesConsecutive(mol::ResidueHandle r1, 
                                       mol::ResidueHandle r2)
{
  if (!r1 || !r2) return false; // protect against invalid handles
  if(r1.GetChain() != r2.GetChain()) return false;
  return r2.GetNumber().GetNum()==r1.GetNumber().GetNum()+1 ||
         r2.GetNumber().GetInsCode()==r1.GetNumber().NextInsertionCode().GetInsCode();
}


String Processor::OptionsToString() const {
  std::stringstream ss;
  ss << "check_bond_feasibility=" << (check_bond_feasibility_ ? "True" : "False") << ", "
     << "assign_torsions=" << (assign_torsions_ ? "True" : "False") << ", "
     << "connect=" << (connect_ ? "True" : "False") << ", "
     << "connect_peptides=" << (connect_aa_ ? "True" : "False")  << ", "
     << "zero_occ_treatment='" << StringFromConopAction(zero_occ_treatment_) << "'";
  return ss.str();
}


}}
