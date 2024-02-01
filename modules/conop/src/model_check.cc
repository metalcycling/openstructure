#include "model_check.hh"
#include "amino_acids.hh"

using namespace ost::mol;

namespace ost { namespace conop {

void Checker::CheckForCompleteness(bool require_hydrogens)
{
  for (ResidueHandleList::const_iterator i=residues_.begin(), 
       e = residues_.end(); i!=e; ++i) {
    ResidueHandle res=*i;
    String anames="";
    CompoundPtr compound=lib_->FindCompound(res.GetName(),Compound::PDB);
    if (!compound) {
      if (checked_unk_res_) {
        continue;
      }
      diags_.AddDiag(DIAG_UNK_RESIDUE, "unknown residue %0")
            .AddResidue(res);
      continue;
    }
    int missing_atoms=0;
    const AtomSpecList& atom_specs=compound->GetAtomSpecs();
    for (AtomSpecList::const_iterator
         j=atom_specs.begin(), e2=atom_specs.end(); j!=e2; ++j) {
      if (j->is_leaving) {
        continue;
      }
      if (!require_hydrogens && (j->element=="H" || j->element=="D")) {
         continue;
      }
      if (!res.FindAtom(j->name) && !res.FindAtom(j->alt_name)) {
        if (!anames.empty()) {
          anames+=", ";
        }
        anames+="'"+j->name+"'";
        missing_atoms+=1;
      }
    }
    if (missing_atoms>0) {
      diags_.AddDiag(DIAG_MISSING_ATOM, "residue %0 is missing %1 atom%s1: %2")
            .AddResidue(res)
            .AddInt(missing_atoms)
            .AddString(anames);
    }
  }
  checked_unk_res_=true;
}

mol::AtomHandleList Checker::GetHydrogens()
{
  AtomHandleList hydlist;
  for (ResidueHandleList::const_iterator i=residues_.begin(), 
       e = residues_.end(); i!=e; ++i) {
    ResidueHandle res=*i;
    CompoundPtr compound=lib_->FindCompound(res.GetName(),Compound::PDB);
    if (!compound) {
      continue;
    }
    AtomHandleList atoms=res.GetAtomList();
    for (AtomHandleList::const_iterator
         j=atoms.begin(), e2=atoms.end(); j!=e2; ++j) {
      int specindx=compound->GetAtomSpecIndex(j->GetName());
      if (specindx!=-1) {  
        if (compound->GetAtomSpecs()[specindx].element=="H" || compound->GetAtomSpecs()[specindx].element=="D") {
           hydlist.push_back(*j);
        }
      }
    }
  }
  return hydlist;
}

mol::AtomHandleList Checker::GetZeroOccupancy()
{
  AtomHandleList zerolist;
  for (ResidueHandleList::const_iterator i=residues_.begin(), 
       e = residues_.end(); i!=e; ++i) {
    ResidueHandle res=*i;
    AtomHandleList atoms=res.GetAtomList();
    for (AtomHandleList::const_iterator j=atoms.begin(), e2=atoms.end(); j!=e2; ++j) {
      if (j->GetOccupancy()==0.0) {
        zerolist.push_back(*j);
      }
    }
  }
  return zerolist;
}

void Checker::CheckForNonStandard()
{
  for (ResidueHandleList::const_iterator i=residues_.begin(), 
       e = residues_.end(); i!=e; ++i) {
    ResidueHandle res=*i;
    CompoundPtr compound=lib_->FindCompound(res.GetName(),Compound::PDB);
    if (!compound) {
      if (checked_unk_res_) {
        continue;
      }
      diags_.AddDiag(DIAG_UNK_RESIDUE, "unknown residue %0")            
            .AddResidue(res);
      continue;
    }
    if (ResidueToAminoAcid(res)==XXX) {
      diags_.AddDiag(DIAG_NONSTD_RESIDUE, "%0 is not a standard amino acid")            
            .AddResidue(res);
    }
  }
  checked_unk_res_=true;
}

void Checker::CheckForUnknownAtoms()
{
  for (ResidueHandleList::const_iterator i=residues_.begin(), 
       e = residues_.end(); i!=e; ++i) {
    ResidueHandle res=*i;
    String anames="";
    CompoundPtr compound=lib_->FindCompound(res.GetName(),Compound::PDB);
    if (!compound) {
      if (checked_unk_res_) {
        continue;
      }
      diags_.AddDiag(DIAG_UNK_RESIDUE, "unknown residue %0")            
            .AddResidue(res);
      continue;
    }
    AtomHandleList atoms=res.GetAtomList();
    const AtomSpecList& atom_specs=compound->GetAtomSpecs();
    for (AtomHandleList::const_iterator
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
        diags_.AddDiag(DIAG_UNK_ATOM, "residue %0 contains unknown atom %1")
              .AddResidue(res)
              .AddAtom(*j);
      }
    }
  }
  checked_unk_res_=true;
}
}} /* ost::conop */
