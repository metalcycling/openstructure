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

#ifndef OST_MM_INTERACTION_HH
#define OST_MM_INTERACTION_HH

#include <vector>
#include <fstream>

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#include <ost/io/binary_data_source.hh>
#include <ost/io/binary_data_sink.hh>
#include <ost/io/io_exception.hh>
#include <ost/message.hh>
#include <ost/mol/atom_handle.hh>
#include <ost/mol/residue_handle.hh>


namespace ost { namespace mol{ namespace mm{

class Interaction;
typedef boost::shared_ptr<Interaction> InteractionPtr;

enum FuncType{
  Invalid,
  HarmonicBond,
  HarmonicAngle,
  UreyBradleyAngle,
  PeriodicDihedral,
  PeriodicImproper,
  HarmonicImproper,
  CMap,
  LJ,
  LJPair,
  GBSA,
  DistanceConstraint,
  Exclusion,
  HarmonicPositionRestraint,
  HarmonicDistanceRestraint
};

class Interaction{

public:
  Interaction(FuncType func_type);

  void SetTypes(std::vector<String> types);

  void SetNames(std::vector<String> names);

  void SetParam(std::vector<Real>& parameters);

  std::vector<String> GetTypes() const { return atom_types_; }

  std::vector<String> GetNames() const { return atom_names_; }

  std::vector<Real> GetParam() const { return parameters_; }

  ost::mol::AtomHandleList GetAtoms(const ost::mol::ResidueHandle& res) const;

  FuncType GetFuncType() const { return func_type_; }

  bool ReplaceAtom(const String& name, const String& new_name, const String& new_type);

  bool MatchTypes(const std::vector<String>& atom_types) const;

  bool MatchNames(const std::vector<String>& atom_names) const;

  bool HasName(const String& name) const;

  bool HasType(const String& type) const;

  bool IsParametrized() const { return set_parameters_; }

  bool HasTypeWildcard() const { return has_type_wildcard_; }

  bool HasNameWildcard() const { return has_name_wildcard_; }

  template <typename DS>
  void Serialize(DS& ds){
    ds & set_parameters_;
    ds & has_type_wildcard_;
    ds & has_name_wildcard_;

    if(ds.IsSource()){
      int num_types = 0;
      int num_names = 0;
      int num_param = 0;
      ds & num_types;
      ds & num_names;
      ds & num_param;
      for(int i = 0; i < num_types; ++i){
        String type;
        ds & type;
        atom_types_.push_back(type);
      }
      for(int i = 0; i < num_names; ++i){
        String name;
        ds & name;
        atom_names_.push_back(name);
      }
      for(int i = 0; i < num_param; ++i){
        Real param;
        ds & param;
        parameters_.push_back(param);
      }
    }
    else{
      int atom_types_size = atom_types_.size();
      int atom_names_size = atom_names_.size();
      int parameters_size = parameters_.size();
      ds & atom_types_size;
      ds & atom_names_size;
      ds & parameters_size;

      for(std::vector<String>::iterator i = atom_types_.begin();
          i != atom_types_.end(); ++i){
        ds & *i;
      }
      for(std::vector<String>::iterator i = atom_names_.begin();
          i != atom_names_.end(); ++i){
        ds & *i;
      }
      for(std::vector<Real>::iterator i = parameters_.begin();
          i != parameters_.end(); ++i){
        ds & *i;
      }
    }
  }

private:

  bool CheckSetNamesTypes(std::vector<String>& types);
  bool CheckSetParam(std::vector<Real>& param);

  FuncType func_type_;
  bool set_parameters_;
  bool has_type_wildcard_;
  bool has_name_wildcard_;
  std::vector<Real> parameters_;
  std::vector<String> atom_types_;
  std::vector<String> atom_names_;
};

}}} //ns

#endif
