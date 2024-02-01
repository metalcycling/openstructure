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
  Author: Marco Biasini, Juergen Haas
 */

#include <ost/message.hh>
#include <ost/log.hh>
#include <ost/dyn_cast.hh>
#include <ost/conop/conop.hh>
#include <ost/mol/mol.hh>
#include <ost/mol/alg/construct_cbeta.hh>
#include <ost/conop/compound_lib.hh>
#include "nonstandard.hh"
using namespace ost::mol;
using namespace ost;
using namespace ost::conop;

namespace ost { namespace mol { namespace alg {


bool CopyResidue(ResidueHandle src_res, ResidueHandle dst_res, XCSEditor& edi)
{
  // first let's get our hands on the component library
  conop::CompoundLibPtr comp_lib=conop::Conopology::Instance().GetDefaultLib();
  return CopyResidue(src_res, dst_res, edi, comp_lib); 
}

bool CopyResidue(ResidueHandle src_res, ResidueHandle dst_res, 
                 XCSEditor& edi, CompoundLibPtr comp_lib)
{
  bool has_cbeta = false;
  bool ret;
  char parent_src = (comp_lib->FindCompound(src_res.GetName(),
                                            Compound::PDB))->GetOneLetterCode();
  char parent_dst = (comp_lib->FindCompound(dst_res.GetName(),
                                            Compound::PDB))->GetOneLetterCode();
  if (parent_src == parent_dst) {
    ret = CopyConserved(src_res, dst_res, edi, has_cbeta, comp_lib);
  } else {
    ret = CopyNonConserved(src_res, dst_res, edi, has_cbeta);
  }
  // insert Cbeta, unless dst residue is glycine.
  if (!has_cbeta && dst_res.GetName() != "GLY") {
    try {
      geom::Vec3 cbeta_pos = mol::alg::CBetaPosition(dst_res);
      edi.InsertAtom(dst_res, "CB", cbeta_pos, "C");
    } catch (ost::Error& e) {
      LOG_WARNING("Issue in CB reconstruction while copying residue "
                  << src_res.GetQualifiedName() << " to "
                  << dst_res.GetQualifiedName() << ": " << e.what()
                  << ". Skipped reconstruction.");
    }
  }
  return ret;
}


bool CopyConserved(ResidueHandle src_res, ResidueHandle dst_res, XCSEditor& edi,
                   bool& has_cbeta)
{
  // first let's get our hands on the component library
  conop::CompoundLibPtr comp_lib=conop::Conopology::Instance().GetDefaultLib();
  return CopyConserved(src_res,dst_res,edi,has_cbeta,comp_lib);
}

bool CopyConserved(ResidueHandle src_res, ResidueHandle dst_res, XCSEditor& edi,
                   bool& has_cbeta, CompoundLibPtr comp_lib)
{
  // check if the residue name of dst and src are the same. In the easy 
  // case, the two residue names match and we just copy over all atoms.
  // If they don't match, we are dealing with modified residues.

  //~ return copy_conserved(src_res, dst_res, edi, has_cbeta);

  if (dst_res.GetName()==src_res.GetName()) {
    return CopyIdentical(src_res, dst_res, edi, has_cbeta);
  } else if (src_res.GetName()=="MSE") {
    return CopyMSE(src_res, dst_res, edi, has_cbeta);
  } else {
    return CopyModified(src_res, dst_res, edi, has_cbeta, comp_lib);
  }
}


bool CopyIdentical(ResidueHandle src_res, ResidueHandle dst_res, XCSEditor& edi, 
                   bool& has_cbeta)
{
  //~ return copy_identical();
  AtomHandleList atoms=src_res.GetAtomList();
  for (AtomHandleList::const_iterator 
       a=atoms.begin(), e=atoms.end(); a!=e; ++a) {
    if (a->GetName()=="CB") {
      has_cbeta=true;
    }
    if (a->GetElement()=="D" || a->GetElement()=="H") {
      continue;
    }
    edi.InsertAtom(dst_res, a->GetName(), a->GetPos(), a->GetElement(), 
                   a->GetOccupancy(), a->GetBFactor());
  }
  return true;
}


bool CopyMSE(ResidueHandle src_res, ResidueHandle dst_res, XCSEditor& edi, 
             bool& has_cbeta)
{
  // selenium methionine is easy. We copy all atoms and substitute the SE 
  // with SD
  AtomHandleList atoms=src_res.GetAtomList();
  for (AtomHandleList::const_iterator 
       a=atoms.begin(), e=atoms.end(); a!=e; ++a) {
    if (a->GetName()=="CB") {
      has_cbeta=true;
    }
    if (a->GetElement()=="D" || a->GetElement()=="H") {
      continue;
    }
    if (a->GetName()=="SE") {
      edi.InsertAtom(dst_res, "SD", a->GetPos(), "S", 
                     a->GetOccupancy(), a->GetBFactor());             
    } else {
      edi.InsertAtom(dst_res, a->GetName(), a->GetPos(), a->GetElement(), 
                     a->GetOccupancy(), a->GetBFactor());            
    }
  }
  return true;
}

bool CopyModified(ResidueHandle src_res, ResidueHandle dst_res, 
                  XCSEditor& edi, bool& has_cbeta, CompoundLibPtr comp_lib)
{
  // FIXME: for now this functions ignores chirality changes of sidechain 
  //        chiral atoms. To detect those changes, we would need to store the 
  //        chirality of those centers in the compound library.
  
  // For now, this function just handles cases where the src_res contains 
  // additional atoms, but the dst_res doesn't contain any atoms the src_res
  // doesn't have. If these two requirements are not met, we fall back to
  // CopyNonConserved.
  
  CompoundPtr src_compound=comp_lib->FindCompound(src_res.GetName(), 
                                                  Compound::PDB);
  if (!src_compound) {
    // OK, that's bad. Let's copy the backbone and be done with it!
    return CopyNonConserved(src_res, dst_res, edi, has_cbeta);
  }
  // since dst_res is one of the 20 amino acids, we don't have to check for 
  // existence of the compound. We know it is there!
  CompoundPtr dst_compound=comp_lib->FindCompound(dst_res.GetName(),
                                                  Compound::PDB);
  std::set<String> dst_atoms;
  std::set<String> src_atoms;
  // to compare the atoms of dst_res with those of src_res, let's create two 
  // sets containing all heavy atoms.
  for (AtomSpecList::const_iterator i=dst_compound->GetAtomSpecs().begin(), 
       e=dst_compound->GetAtomSpecs().end(); i!=e; ++i) {
    if (i->element=="H" || i->element=="D") {
      continue;
    }
    dst_atoms.insert(i->name);
  }
  for (AtomSpecList::const_iterator i=src_compound->GetAtomSpecs().begin(), 
       e=src_compound->GetAtomSpecs().end(); i!=e; ++i) {
    if (i->element=="H" || i->element=="D") {
      continue;
    }
    src_atoms.insert(i->name);
  }
  for (std::set<String>::const_iterator i=dst_atoms.begin(), 
       e=dst_atoms.end(); i!=e; ++i) {
    if (src_atoms.find(*i)==src_atoms.end()) {
      return CopyNonConserved(src_res, dst_res, edi, has_cbeta);
    }
  }
  // Muahaha, all is good. Let's copy the atoms. Muahaha
  AtomHandleList atoms=src_res.GetAtomList();
  for (AtomHandleList::const_iterator 
       a=atoms.begin(), e=atoms.end(); a!=e; ++a) {
    if (a->GetName()=="CB") {
      if (dst_res.GetName()=="GLY") {
        continue;
      }
      has_cbeta=true;
    }
    if (a->GetElement()=="D" || a->GetElement()=="H") {
      continue;
    }
    if (dst_atoms.find(a->GetName())==dst_atoms.end()) {
      continue;
    }
    edi.InsertAtom(dst_res, a->GetName(), a->GetPos(), a->GetElement(), 
                   a->GetOccupancy(), a->GetBFactor());
  }
  return true;
}


bool CopyNonConserved(ResidueHandle src_res, ResidueHandle dst_res, 
                      XCSEditor& edi, bool& has_cbeta) 
{  
  AtomHandleList atoms=src_res.GetAtomList();
  for (AtomHandleList::const_iterator 
       a=atoms.begin(), e=atoms.end(); a!=e; ++a) {
    String name=a->GetName();
    if (name=="CA" || name=="N" || name=="C" || name=="O" || name=="CB") {
      if (name=="CB") {
        if (dst_res.GetName()=="GLY") {
          continue;
        }
        has_cbeta=true;
      }
      if (a->GetElement()=="D" || a->GetElement()=="H") {
        continue;
      }
      edi.InsertAtom(dst_res, a->GetName(), a->GetPos(), a->GetElement(), 
                     a->GetOccupancy(), a->GetBFactor());     
    }
  }
  return false;
}





}}}
