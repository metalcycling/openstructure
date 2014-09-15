#ifndef OST_GROMACS_READER_HH
#define OST_GROMACS_READER_HH

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/unordered_map.hpp>

#include <ost/base.hh>
#include <ost/io/io_exception.hh>
#include <ost/mol/mm/forcefield.hh>
#include <ost/mol/mm/mm_interaction.hh>
#include <ost/mol/mm/gromacs_block_modifiers.hh>


namespace ost { namespace mol{ namespace mm{

class GromacsData;
class CHARMMData;
class GromacsReader;
typedef boost::shared_ptr<GromacsData> GromacsDataPtr;
typedef boost::shared_ptr<GromacsReader> GromacsReaderPtr;
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


class GromacsReader {
public:

  GromacsReader(const String& base_dir);

  void SetPreprocessorDefinition(const String& def) { preprocessor_.SetDefinition(def); }

  void ReadGromacsForcefield();

  ForcefieldPtr GetForcefield() { return ff_;}

  void SetForcefield(ForcefieldPtr ff) { ff_ = ff; }

  void ReadResidueDatabase(const String& basename);

  void ReadITP(const String& basename);

  void ReadCHARMMPRM(const String& basename);

  void ReadCHARMMRTF(const String& basename);

private:

  MMInteractionPtr ParseBond(const std::vector<String>& data, 
                             bool type_definition, 
                             FuncType functype = Invalid);

  MMInteractionPtr ParseAngle(const std::vector<String>& data, 
                              bool type_definition, 
                              FuncType functype = Invalid);

  MMInteractionPtr ParseDihedral(const std::vector<String>& data, 
                                 bool type_definition, 
                                 FuncType functype = Invalid);

  MMInteractionPtr ParseCMap(const std::vector<String>& data, 
                             bool type_definition, 
                             FuncType functype = Invalid);

  MMInteractionPtr ParseLJ(const std::vector<String>& data, 
                           bool type_definition, 
                           FuncType functype = Invalid);

  MMInteractionPtr ParseLJPair(const std::vector<String>& data, 
                               bool type_definition, 
                               FuncType functype = Invalid);

  MMInteractionPtr ParseConstraint(const std::vector<String>& data, 
                                   bool type_definition, 
                                   FuncType functype = Invalid);

  MMInteractionPtr ParseGenborn(const std::vector<String>& data, 
                                bool type_definition, 
                                FuncType functype = Invalid);

  MMInteractionPtr ParseExclusion(const std::vector<String>& data,
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
