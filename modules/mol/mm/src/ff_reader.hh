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

#ifndef OST_MM_FORCEFIELD_READER_HH
#define OST_MM_FORCEFIELD_READER_HH

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>

#include <ost/base.hh>
#include <ost/io/io_exception.hh>
#include <ost/mol/mm/forcefield.hh>
#include <ost/mol/mm/interaction.hh>
#include <ost/mol/mm/gromacs_block_modifiers.hh>


namespace ost { namespace mol{ namespace mm{

class GromacsData;
class CHARMMData;
class FFReader;
typedef boost::shared_ptr<GromacsData> GromacsDataPtr;
typedef boost::shared_ptr<FFReader> FFReaderPtr;
typedef boost::shared_ptr<CHARMMData> CHARMMDataPtr;

class GromacsData{
public:
  static GromacsDataPtr Instance();
  int GetKeywordIndex(const String& keyword);
  String ConvertToStandard(const String& res_name, const String& atom_name);
  bool ConversionExists(const String& res_name);

private:
  GromacsData();
  GromacsData(const GromacsData&);
  GromacsDataPtr operator=(const GromacsDataPtr&);
  static GromacsDataPtr instance_;
  boost::unordered_map<String,int> keyword_map_;
  boost::unordered_map<String, std::vector<std::pair<String,String> > > renaming_to_standard_;

};

class CHARMMData{
public:
  static CHARMMDataPtr Instance();
  int GetKeywordIndex(const String& keyword);

private:
  CHARMMData();
  CHARMMData(const CHARMMData&);
  CHARMMDataPtr operator=(const CHARMMDataPtr&);
  static CHARMMDataPtr instance_;
  boost::unordered_map<String,int> keyword_map_;
};

class MMPreprocessor{

public:

  MMPreprocessor(const String& basepath): basepath_(basepath) { }

  std::vector<std::vector<String> > Process(const String& filename);

  void SetDefinition(const String& def) { defines_.insert(def); }

  boost::filesystem::path GetBasedir() { return basepath_; }

private:

  //function, that can recursively resolve ifdef / ifndef statements
  void ResolveIFDEF(std::vector<std::vector<String> >& file_content, int line_counter);

  //simply reads a file, cuts it into pieces and removes comments marked by '*' and ';'
  std::vector<std::vector<String> > ReadFile(const String& filename);

  std::map<String,std::vector<String> > definitions_;
  std::set<String> defines_;
  boost::filesystem::path basepath_;
};


class FFReader {
public:

  FFReader(const String& base_dir);

  void SetPreprocessorDefinition(const String& def) { preprocessor_.SetDefinition(def); }

  void ReadGromacsForcefield();

  ForcefieldPtr GetForcefield() { return ff_;}

  void SetForcefield(ForcefieldPtr ff) { ff_ = ff; }

  void ReadResidueDatabase(const String& basename);

  void ReadITP(const String& basename);

  void ReadCHARMMPRM(const String& basename);

  void ReadCHARMMRTF(const String& basename);

private:

  InteractionPtr ParseBond(const std::vector<String>& data, 
                             bool type_definition, 
                             FuncType functype = Invalid);

  InteractionPtr ParseAngle(const std::vector<String>& data, 
                              bool type_definition, 
                              FuncType functype = Invalid);

  InteractionPtr ParseDihedral(const std::vector<String>& data, 
                                 bool type_definition, 
                                 FuncType functype = Invalid);

  InteractionPtr ParseCMap(const std::vector<String>& data, 
                             bool type_definition, 
                             FuncType functype = Invalid);

  InteractionPtr ParseLJ(const std::vector<String>& data, 
                           bool type_definition, 
                           FuncType functype = Invalid);

  InteractionPtr ParseLJPair(const std::vector<String>& data, 
                               bool type_definition, 
                               FuncType functype = Invalid);

  InteractionPtr ParseConstraint(const std::vector<String>& data, 
                                   bool type_definition, 
                                   FuncType functype = Invalid);

  InteractionPtr ParseGenborn(const std::vector<String>& data, 
                                bool type_definition, 
                                FuncType functype = Invalid);

  InteractionPtr ParseExclusion(const std::vector<String>& data,
                                  bool type_definition,
                                  FuncType functype = Invalid);

  BuildingBlockPtr BlockFromRTP(const std::vector<std::vector<String> >& data);

  BuildingBlockPtr BlockFromITP(const std::vector<std::vector<String> >& data);

  TerminiConstructorPtr ParseTermini(const std::vector<std::vector<String> >& data);

  BlockModifierPtr ParseBlockModifier(const std::vector<std::vector<String> >& data);

  void ParseHydrogenRule(const std::vector<String>& data, GromacsHydrogenConstructor& constructor);

  void ParseTerminiReplaceRule(const std::vector<String>& data, GromacsBlockModifier& constructor);

  void ParseTerminiAddRule(const std::vector<String>& data1, const std::vector<String>& data2,
                           GromacsBlockModifier& constructor);

  //Reader functions for the mandatory forcefield files
  void ParseForcefield(std::vector<std::vector<String> >& content);
  void ParseAtomTypes(std::vector<std::vector<String> >& content);

  //Reader functions for CHARMM stuff
  void ParseCHARMMPRM(std::vector<std::vector<String> >& content);

  void ParseCHARMMRTF(std::vector<std::vector<String> >& content);

  //Reader functions for all different residue database files
  void ParseRTP(std::vector<std::vector<String> >& content);
  void ParseARN(std::vector<std::vector<String> >& content);
  void ParseHDB(std::vector<std::vector<String> >& content);
  void ParseNTDB(std::vector<std::vector<String> >& content);
  void ParseCTDB(std::vector<std::vector<String> >& content);
  void ParseVSD(std::vector<std::vector<String> >& content);
  void ParseRtoB(std::vector<std::vector<String> >& content);
  //Reader function for single molecule itp files
  void ParseITP(std::vector<std::vector<String> >& content);

  boost::unordered_map<String, std::vector<std::pair<String,String> > > atom_renaming_ff_specific_;
  boost::unordered_map<String, ResidueNamesPtr> res_renaming_ff_specific_;

  std::vector<FuncType> ff_bonded_types_;

  MMPreprocessor preprocessor_;
  ForcefieldPtr ff_;

  //following part is ugly...
  //data that is read during the residue datababase parsing process gets stored in there
  std::vector<FuncType> bonded_types_;
  std::vector<String> read_residues_;
};


}}}//ns

#endif
