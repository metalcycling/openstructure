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
#ifndef OST_IO_OMF_HH
#define OST_IO_OMF_HH

#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <boost/iostreams/filtering_stream.hpp>
#include <ost/mol/entity_handle.hh>
#include <ost/geom/mat4.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/mmcif_info.hh>

namespace ost { namespace io {

const int OMF_VERSION = 3;

class ChainData;
class OMF;
typedef boost::shared_ptr<OMF> OMFPtr;
typedef boost::shared_ptr<ChainData> ChainDataPtr;

struct SidechainAtomRule {
  int sidechain_atom_idx;
  int anchor_idx[3];
  Real bond_length;
  Real angle;
  // 0: chi1, 1: chi2, 2: chi3, 3: chi4, 4: 0.0
  int dihedral_idx;
  // the value of the dihedral above will be added to base_dihedral to get
  // the final diheral angle. If you want to have the effect of chi3 + M_PI
  // you define dihedral_idx as 2 and base_dihedral = M_PI.
  Real base_dihedral;
};

struct ChiDefinition{
  int idx_one;
  int idx_two;
  int idx_three;
  int idx_four;
};

struct ResidueDefinition {

  ResidueDefinition() { };

  ResidueDefinition(const ost::mol::ResidueHandle& res);

  bool operator==(const ResidueDefinition& other) const { 
    return (name == other.name && 
            olc == other.olc &&
            chem_type == other.chem_type &&
            chem_class == other.chem_class &&
            anames == other.anames &&
            elements == other.elements &&
            is_hetatm == other.is_hetatm &&
            bonds == other.bonds &&
            bond_orders == other.bond_orders);
  }

  bool operator!=(const ResidueDefinition& other) const { 
    return !(*this == other);
  }

  void ToStream(std::ostream& stream) const;

  void FromStream(std::istream& stream);

  int GetIdx(const String& aname) const;

  const std::set<int>& GetRotamericAtoms() const;

  const std::vector<ChiDefinition>& GetChiDefinitions() const;

  const std::vector<SidechainAtomRule>& GetSidechainAtomRules() const;

  int GetNChiAngles() const;

  void _InitIdxMapper() const;

  void _AddChiDefinition(int idx_one, int idx_two, int idx_three,
                         int idx_four);

  void _AddAtomRule(int a_idx, int anch_one_idx,
                    int anch_two_idx, int anch_three_idx, 
                    Real bond_length, Real angle, int dihedral_idx, 
                    Real base_dihedral);

  String name;
  char olc;
  char chem_type;
  char chem_class;
  std::vector<String> anames;
  std::vector<String> elements;
  std::vector<bool> is_hetatm;
  std::vector<int> bonds;
  std::vector<int> bond_orders;
  mutable std::map<String, int> idx_mapper;
  std::set<int> rotameric_atoms;
  std::vector<ChiDefinition> chi_definitions;
  std::vector<SidechainAtomRule> sidechain_atom_rules;
  std::set<int> critical_sidechain_angles;
};

struct ChainData {

  ChainData(): ch_name(""), chain_type(ost::mol::CHAINTYPE_UNKNOWN) { }

  ChainData(const ost::mol::ChainHandle& chain,
            const std::vector<ResidueDefinition>& residue_definitions,
            const std::unordered_map<unsigned long, int>& res_idx_map,
            const std::vector<std::pair<unsigned long, unsigned long> >& 
            inter_residue_bonds,
            const std::vector<int>& inter_residue_bond_orders,
            std::unordered_map<long, int>& atom_idx_mapper);

  void ToStream(std::ostream& stream,
                const std::vector<ResidueDefinition>& res_def,
                Real max_error, bool avg_bfactors, bool round_bfactors,
                bool skip_ss) const;

  void FromStream(std::istream& stream,
                  const std::vector<ResidueDefinition>& res_def,
                  int version, Real max_error, bool avg_bfactors,
                  bool round_bfactors, bool skip_ss);

  // chain features
  String ch_name;
  ost::mol::ChainType chain_type;

  // residue features
  std::vector<int> res_def_indices;
  std::vector<int> rnums;
  std::vector<char> insertion_codes;
  std::vector<char> sec_structures;

  // atom features    
  std::vector<Real> occupancies;
  std::vector<Real> bfactors;
  geom::Vec3List positions;

  // bond features - only for bonds that are inter-residue
  // e.g. peptide bonds
  std::vector<int> bonds;
  std::vector<int> bond_orders;
};


class DefaultPepLib{
public:
  static DefaultPepLib& Instance() {
    static DefaultPepLib instance;
    return instance;
  }
  std::vector<ResidueDefinition> residue_definitions;

private:
  DefaultPepLib();
  DefaultPepLib(DefaultPepLib const& copy); 
  DefaultPepLib& operator=(DefaultPepLib const& copy);
};

class OMF {

public:

  enum OMFOption {DEFAULT_PEPLIB = 1, AVG_BFACTORS = 2, ROUND_BFACTORS = 4,
                  SKIP_SS = 8, INFER_PEP_BONDS = 16};

  bool OptionSet(OMFOption opt) const {
    return (opt & options_) == opt;
  }

  static OMFPtr FromEntity(const ost::mol::EntityHandle& ent,
                           Real max_error = 0.0,
                           uint8_t options = 0);

  static OMFPtr FromFile(const String& fn);

  static OMFPtr FromString(const String& s);

  void ToFile(const String& fn) const;

  String ToString() const;

  ost::mol::EntityHandle GetAU() const;

  ost::mol::EntityHandle GetEntity() const {
    return this->GetAU();
  }

  ost::mol::EntityHandle GetAUChain(const String& name) const;

  ost::mol::EntityHandle GetEntityChain(const String& name) const {
    return this->GetAUChain(name);
  }

  int GetVersion() const { return version_; }

  static int GetCurrentOMFVersion() { return OMF_VERSION; }

  Real GetMaxError() const { return 0.001 * max_error_; }

  // data access without requirement of generating a full
  // OpenStructure entity

  String GetName() const { return name_; }

  std::vector<String> GetChainNames() const;

  const geom::Vec3List& GetPositions(const String& cname) const;

  const std::vector<Real>& GetBFactors(const String& cname) const;

  std::vector<Real> GetAvgBFactors(const String& cname) const;

  String GetSequence(const String& cname) const;

private:
  // only construct with static functions
  OMF(): options_(0) { }

  void ToStream(std::ostream& stream) const;

  void FromStream(std::istream& stream);

  void FillChain(const ChainDataPtr data, ost::mol::XCSEditor& ed,
                 ost::mol::ChainHandle& chain) const;

  String name_;
  uint16_t max_error_;
  std::vector<ResidueDefinition> residue_definitions_;
  std::map<String, ChainDataPtr> chain_data_;

  // bond features - only for bonds that are inter-chain
  // given n bonds, bond_chain_names_ and bond_atoms_ have length 2*n and are
  // organized as follows: 
  // [bond1_at1_x, bond1_at2_x, ..., bondn_at1_x, bondn_at2_x]
  // bond_orders_ on the other hand has length n
  std::vector<String> bond_chain_names_;
  std::vector<int> bond_atoms_;
  std::vector<int> bond_orders_;

  // bitfield with options
  uint8_t options_;

  int version_;
};

}} //ns

#endif
