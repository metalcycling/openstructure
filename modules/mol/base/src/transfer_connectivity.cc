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

#include "transfer_connectivity.hh"
#include "entity_visitor.hh"
#include "residue_handle.hh"
#include "atom_handle.hh"
#include "bond_handle.hh"
#include "xcs_editor.hh"


namespace ost { namespace mol {

namespace {

class ConnectivityTransferVisitor : public EntityVisitor {
public:
  ConnectivityTransferVisitor(EntityHandle dst,
                              const std::map<ResidueHandle, ResidueHandle>* to_from):
    edi_(dst.EditXCS(mol::BUFFERED_EDIT)), to_from_(to_from) {
  }

  bool VisitResidue(const ResidueHandle& dest_res) {
    std::map<ResidueHandle, ResidueHandle>::const_iterator i = to_from_->find(dest_res);
    if (i == to_from_->end()) {
      return false;
    }
    const std::pair<ResidueHandle, ResidueHandle>& ft = *i;
    if (!ft.second.IsValid()) {
      return false;
    }

    // define a map from atom name to atom handle for the dest residue.
    std::map<String, AtomHandle> name_to_handle;
    AtomHandleList dest_atoms = dest_res.GetAtomList();
    for (AtomHandleList::iterator j = dest_atoms.begin(), e = dest_atoms.end(); j !=e ; ++j) {
      name_to_handle[j->GetName()] = *j;
    }

    AtomHandleList src_atoms =  ft.second.GetAtomList();
    for (AtomHandleList::iterator j = src_atoms.begin(), e = src_atoms.end(); j !=e ; ++j) {
      BondHandleList bonds = j->GetBondList();
      for (BondHandleList::iterator k = bonds.begin(), e2 = bonds.end(); k !=e2; ++k) {
        AtomHandle first = k->GetFirst();
        AtomHandle second = k->GetSecond();
        AtomHandle dst_frst = this->GetDestAtomForSrcAtom(first, ft.second, ft.first,
                                                          name_to_handle);
        AtomHandle dst_scnd = this->GetDestAtomForSrcAtom(second, ft.second, ft.first,
                                                          name_to_handle);
        // add bond...
        if (!dst_frst || !dst_scnd) {
          continue;
        }
        edi_.Connect(dst_frst, dst_scnd);
      }
    }
    return false;
  }

  bool CheckInsertionCode(ResidueHandle chk_res, ResidueHandle src_res,
                          ResidueHandle dst_res) {
    // This is a hack to make GetDestAtomForSrcAtom work with mol::alg::PDBize,
    // maybe with insertion codes in general...
    // Depending on the value of min_polymer_size, PDBize puts small polymers in
    // the ligand chain called '_'. The original chain is annotated by all
    // residues having the same residue number but different insertion codes.
    // That can lead to the same issue described further down for branched
    // entities. Basically the problem is that in branched entities the residues
    // do not need to be connected with their immediate neighbours breaking the
    // original mechanics of GetDestAtomForSrcAtom.
    // check if residue has an insertion code
    if (chk_res.GetNumber().GetInsCode() != '\0') {
      // check if the original residue also had no inscode
      if (src_res.GetNumber().GetInsCode() == '\0') {
        // There must be another residue with the same number and since
        // this function is only used by PDBize, that residue originates
        // from the same chain as dst_res and may be the correct partner.
        if (chk_res.GetNumber().GetNum() != dst_res.GetNumber().GetNum()) {
          return false;
        }
      }
    }
    return true;
  }

  AtomHandle GetDestAtomForSrcAtom(AtomHandle src_atom, ResidueHandle src_res,
                                   ResidueHandle dst_res,
                                   const std::map<String, AtomHandle>& name_to_atom) {
    ResidueHandle r = src_atom.GetResidue();
    if (r == src_res) {
      // fast track...
      std::map<String, AtomHandle>::const_iterator i = name_to_atom.find(src_atom.GetName());
      return i == name_to_atom.end() ? AtomHandle() : i->second;
    }
    // educated guess: we are trying to connect to the previous/next residue
    std::map<ResidueHandle, ResidueHandle>::const_iterator j;
    j = to_from_->find(dst_res.GetPrev());
    if (j != to_from_->end()) {
      if (j->second == r) {
        if (CheckInsertionCode(
                              j->first.FindAtom(src_atom.GetName()).GetResidue(),
                              r, dst_res)) {
          return j->first.FindAtom(src_atom.GetName());
        }
      }
    }
    j = to_from_->find(dst_res.GetNext());
    if (j != to_from_->end()) {
      if (j->second == r) {
        if (CheckInsertionCode(
                              j->first.FindAtom(src_atom.GetName()).GetResidue(),
                              r, dst_res)) {
            return j->first.FindAtom(src_atom.GetName());
          }
      }
    }
    // still nothing. scan linearly through all residues.
    for ( j = to_from_->begin(); j != to_from_->end(); ++j) {
      if (j->second == r) {
        // Check that we are connecting to the same chain, otherwise we need to
        // check that the found residue also connects outwards.
        if (j->first.GetChain() == dst_res.GetChain()) {
          if (CheckInsertionCode(
                              j->first.FindAtom(src_atom.GetName()).GetResidue(),
                              r, dst_res)) {
            return j->first.FindAtom(src_atom.GetName());
          }
        }
        // Using the found residue would mean connecting to a different chain
        // which is unusual. So we make sure that in the residue we transfer
        // from, there is also a bond to a different chain. Otherwise the
        // search continues.
        // For branched mmCIF entities where connected residues may be not
        // direct neighbours since branched entities are non-linear, that gave
        // a problem when a bio unit doubled an oligosaccharide and the search
        // found the bond of the copied oligosaccharide first. In that case the
        // copied oligosaccharide is a different chain and the residue from the
        // same chain comes later in the list.
        BondHandleList jbonds = src_atom.GetBondList();
        for (BondHandleList::iterator k = jbonds.begin(), e2 = jbonds.end();
             k !=e2; ++k) {
          // determine which "side" the src_res sits on
          if (k->GetFirst() == src_atom) {
            if (src_atom.GetResidue().GetChain()
                != k->GetSecond().GetResidue().GetChain()) {
              return j->first.FindAtom(src_atom.GetName());
            }
            continue;
          }
          if (src_atom.GetResidue().GetChain()
              != k->GetFirst().GetResidue().GetChain()) {
            return j->first.FindAtom(src_atom.GetName());
          }
        }
      }
    }
    return AtomHandle();
  }
private:
  XCSEditor              edi_;
  const std::map<ResidueHandle, ResidueHandle>* to_from_;
};

}
void TransferConnectivity(EntityHandle dest,
                          const std::map<ResidueHandle, ResidueHandle>& to_from) {
  ConnectivityTransferVisitor visitor(dest, &to_from);
  dest.Apply(visitor);
}

}}

//  LocalWords:  PDBize ligand
