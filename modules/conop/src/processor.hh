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
#ifndef OST_CONOP_PROCESSOR_HH
#define OST_CONOP_PROCESSOR_HH

#include <ost/mol/entity_handle.hh>
#include "module_config.hh"
#include "diag.hh"
#include "compound.hh"
namespace ost { namespace conop {

typedef enum {
  PDB_DIALECT,
  CHARMM_DIALECT
} Dialect;

enum ConopAction {
  CONOP_WARN = 0,
  CONOP_SILENT,
  CONOP_REMOVE,
  CONOP_REMOVE_ATOM,
  CONOP_REMOVE_RESIDUE,
  CONOP_FATAL
};

class Processor;
typedef boost::shared_ptr<Processor> ProcessorPtr;
// the base class for all options
class DLLEXPORT_OST_CONOP Processor {
public:
  DiagnosticsPtr Process(mol::EntityHandle ent, bool log_diags=true) const;
  virtual ProcessorPtr Copy() const = 0;
  virtual ~Processor() {}
protected:
  virtual void DoProcess(DiagnosticsPtr diags, 
                         mol::EntityHandle ent) const = 0;
  virtual bool BeginProcessing(DiagnosticsPtr diags, 
                               mol::EntityHandle ent) const { return true; }  
  virtual bool EndProcessing(DiagnosticsPtr diags, 
                             mol::EntityHandle ent) const { return true; }
  bool HasUnknownAtoms(mol::ResidueHandle residue, CompoundPtr compound, 
                       bool strict_hydrogens) const;
  void ReorderAtoms(mol::ResidueHandle residue, CompoundPtr compound, 
                    bool fix_element) const;
  void FillResidueProps(mol::ResidueHandle residue, CompoundPtr compound) const;
  void ConnectAtomsOfResidue(mol::ResidueHandle residue, 
                             CompoundPtr compound, bool strict_hydrogens) const;
  void ConnectResidues(mol::ResidueHandle residue, mol::ResidueHandle next) const;
  static bool AreResiduesConsecutive(mol::ResidueHandle a, mol::ResidueHandle b);
  void DistanceBasedConnect(mol::AtomHandle atom) const;
  mol::AtomHandle LocateAtom(const mol::AtomHandleList&, int ordinal) const;
public:
  Processor(bool bf, bool at, bool cn, bool aa, bool ch, ConopAction zo):
    check_bond_feasibility_(bf), assign_torsions_(at), connect_(cn),
    connect_aa_(aa), connect_hetatm_(ch), zero_occ_treatment_(zo) {}
  Processor(): check_bond_feasibility_(false),
    assign_torsions_(true), connect_(true), connect_aa_(true),
    connect_hetatm_(true), zero_occ_treatment_(CONOP_SILENT) {}

  void SetConnect(bool connect) {
    connect_ = connect;
  }
  bool GetConnect() const {
    return connect_;
  }

  void SetAssignTorsions(bool flag) {
    assign_torsions_ = flag;
  }
  bool GetAssignTorsions() const {
    return assign_torsions_;
  }

  bool GetConnectAminoAcids() const {
    return connect_aa_;
  }
  void SetConnectAminoAcids(bool c) {
    connect_aa_ = c;
  }

  bool GetConnectHetatm() const {
    return connect_hetatm_;
  }
  void SetConnectHetatm(bool c) {
    connect_hetatm_ = c;
  }

  bool GetCheckBondFeasibility() const {
    return check_bond_feasibility_;
  }

  void SetCheckBondFeasibility(bool flag) {
    check_bond_feasibility_ = flag;
  }

  ConopAction GetZeroOccTreatment() const {
    return zero_occ_treatment_;
  }
  void SetZeroOccTreatment(ConopAction action) {
    zero_occ_treatment_ = action;
  }

  virtual String ToString() const = 0;
protected:
  String OptionsToString() const; 
private:
  bool check_bond_feasibility_;
  bool assign_torsions_;
  bool connect_;
  bool connect_aa_;
  bool connect_hetatm_;
  ConopAction zero_occ_treatment_;
};

ConopAction DLLEXPORT_OST_CONOP ConopActionFromString(const String& name);

String DLLEXPORT_OST_CONOP StringFromConopAction(ConopAction action);


/// \brief guess element of atom based on name and hetatm flag
String DLLEXPORT_OST_CONOP GuessAtomElement(const String& atom_name, bool hetatm,
                                            int atom_count);

/// \brief guess chemclass based on atoms of residue
mol::ChemClass DLLEXPORT_OST_CONOP GuessChemClass(mol::ResidueHandle res);

/// \brief assigns phi/psi/omega to all residues marked peptide-linking of 
///     the chain
/// 
/// Requires the atoms to be connected
void DLLEXPORT_OST_CONOP AssignBackboneTorsions(mol::ChainHandle chain);
void DLLEXPORT_OST_CONOP AssignBackboneTorsions(mol::ResidueHandleList residues);
void DLLEXPORT_OST_CONOP AssignBackboneTorsions(mol::ResidueHandle prev, 
                                                mol::ResidueHandle res, 
                                                mol::ResidueHandle next);

bool DLLEXPORT_OST_CONOP IsBondFeasible(const mol::AtomHandle&, 
                                        const mol::AtomHandle&);
mol::AtomHandleList DLLEXPORT_OST_CONOP 
GetUnknownAtomsOfResidue(mol::ResidueHandle residue, 
                         CompoundPtr compound,
                         bool strict_hydrogens=false);
}}

#endif
