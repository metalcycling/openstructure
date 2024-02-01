
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
#ifndef OST_MOL_BASE_BUILDER_HH
#define OST_MOL_BASE_BUILDER_HH

#include "module_config.hh"

#include "atom_handle.hh"
#include "residue_handle.hh"
#include "chain_handle.hh"
#include "entity_handle.hh"
#include "bond_handle.hh"
#include "xcs_editor.hh"


namespace ost { namespace mol {

// a helper class to easily create entities
class DLLEXPORT Builder {
public:
   Builder(): ent_(CreateEntity()), edi_(ent_.EditXCS(BUFFERED_EDIT)) {
   }

   // conversion to entity handle
   operator EntityHandle() { return ent_; }

   Builder& Chain(const String& name) {
     chain_ = edi_.InsertChain(name);
     res_ = ResidueHandle();
     atom_ = AtomHandle();
     return *this;
   }

   Builder& Residue(const String& name) {
     res_ = edi_.AppendResidue(chain_, name);
     atom_ = AtomHandle();
     return *this;
   }

   Builder& OneLetterCode(char olc) {
     res_.SetOneLetterCode(olc);
     return *this;
   }

   Builder& Atom(const String& name, const geom::Vec3& pos=geom::Vec3()) {
     edi_.InsertAtom(res_, name, pos);
     return *this;
   }
   Builder& Gly(bool connect=true) {
     this->Residue("GLY");
     this->Atom("N");
     this->Atom("CA");
     this->Atom("C");
     this->Atom("O");
     if (connect) {
       edi_.Connect(res_.FindAtom("N"), res_.FindAtom("CA"));
       edi_.Connect(res_.FindAtom("CA"), res_.FindAtom("C"));
       edi_.Connect(res_.FindAtom("C"), res_.FindAtom("O"));
     }
     return *this;
   }
   Builder& ConnectToPrev() {
     AtomHandle ah = res_.FindAtom("N");
     AtomHandle pa = res_.GetPrev().FindAtom("C");
     if (pa) {
       edi_.Connect(pa, ah);
     }
     return *this;
   }
private:
   EntityHandle  ent_;
   ChainHandle   chain_;
   ResidueHandle res_;
   AtomHandle    atom_;
   XCSEditor     edi_;


};

}}
#endif

