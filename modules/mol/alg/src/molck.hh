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
#ifndef OST_MOL_ALG_MOLCK_HH
#define OST_MOL_ALG_MOLCK_HH

#include <string>
#include <ost/mol/entity_handle.hh>
#include <ost/conop/compound_lib.hh>

namespace {
  inline std::string BoolToString(bool b)
    {
      return b ? "True" : "False";
    }
}

namespace ost { namespace mol{ namespace alg {

struct MolckSettings;

struct MolckSettings{
  
  bool rm_unk_atoms;
  bool rm_non_std;
  bool rm_hyd_atoms;
  bool rm_oxt_atoms;
  bool rm_zero_occ_atoms;
  bool colored;
  bool map_nonstd_res;
  bool assign_elem;

  MolckSettings(bool init_rm_unk_atoms=false,
                bool init_rm_non_std=false,
                bool init_rm_hyd_atoms=true,
                bool init_rm_oxt_atoms=false,
                bool init_rm_zero_occ_atoms=false,
                bool init_colored=false,
                bool init_map_nonstd_res=true,
                bool init_assign_elem=true):
      rm_unk_atoms(init_rm_unk_atoms), // Remove unknown and atoms not following the nomenclature
      rm_non_std(init_rm_non_std), // Remove all residues not one of the 20 standard amino acids
      rm_hyd_atoms(init_rm_hyd_atoms), // Remove hydrogen atoms
      rm_oxt_atoms(init_rm_oxt_atoms), // Remove terminal oxygens
      rm_zero_occ_atoms(init_rm_zero_occ_atoms), // Remove atoms with zero occupancy
      colored(init_colored), // Whether the output should be colored
      map_nonstd_res(init_map_nonstd_res), // Map non standard residues back to standard ones (e.g.: MSE->MET,SEP->SER,etc.)
      assign_elem(init_assign_elem){} // Clean up element column

  public:
    std::string ToString(){
      std::string rep = "MolckSettings(rm_unk_atoms=" + BoolToString(rm_unk_atoms) +
        ", rm_unk_atoms=" + BoolToString(rm_unk_atoms) +
        ", rm_non_std=" + BoolToString(rm_non_std) +
        ", rm_hyd_atoms=" + BoolToString(rm_hyd_atoms) +
        ", rm_oxt_atoms=" + BoolToString(rm_oxt_atoms) +
        ", rm_zero_occ_atoms=" + BoolToString(rm_zero_occ_atoms) +
        ", colored=" + BoolToString(colored) +
        ", map_nonstd_res=" + BoolToString(map_nonstd_res) +
        ", assign_elem=" + BoolToString(assign_elem) +
        ")";
      return rep;
    }

};

void MapNonStandardResidues(ost::mol::EntityHandle& ent,
                            ost::conop::CompoundLibPtr lib,
                            bool reprocess=true);

void RemoveAtoms(ost::mol::EntityHandle& ent,
                 ost::conop::CompoundLibPtr lib,
                 bool rm_unk_atoms,
                 bool rm_non_std,
                 bool rm_hyd_atoms,
                 bool rm_oxt_atoms,
                 bool rm_zero_occ_atoms,
                 bool colored=true,
                 bool reprocess=true);

void CleanUpElementColumn(ost::mol::EntityHandle& ent,
                          ost::conop::CompoundLibPtr lib);

void Molck(ost::mol::EntityHandle& ent,
           ost::conop::CompoundLibPtr lib,
           const MolckSettings& settings,
           bool prune = true);


}}} // namespace

#endif
