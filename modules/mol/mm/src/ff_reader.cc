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

#include <ost/mol/mm/ff_reader.hh>


namespace ost { namespace mol{ namespace mm {


GromacsData::GromacsData(){
  keyword_map_["atoms"] = 0;
  keyword_map_["bonds"] = 1;
  keyword_map_["angles"] = 2;
  keyword_map_["impropers"] = 3;
  keyword_map_["dihedrals"] = 4;
  keyword_map_["exclusions"] = 5;
  keyword_map_["cmap"] = 6;
  keyword_map_["replace"] = 7;
  keyword_map_["moleculetype"] = 8;
  keyword_map_["add"] = 9;
  keyword_map_["delete"] = 10;
  keyword_map_["bondtypes"] = 11;
  keyword_map_["constrainttypes"] = 12;
  keyword_map_["angletypes"] = 13;
  keyword_map_["dihedraltypes"] = 14;
  keyword_map_["atomtypes"] = 15;
  keyword_map_["implicit_genborn_params"] = 16;
  keyword_map_["pairtypes"] = 17;
  keyword_map_["nonbond_params"] = 18;
  keyword_map_["defaults"] = 19;
  keyword_map_["cmaptypes"] = 20;
  keyword_map_["pairs"] = 21;
  keyword_map_["settles"] = 22;


  renaming_to_standard_["ILE"] = std::vector<std::pair<String,String> >();
  renaming_to_standard_["HOH"] = std::vector<std::pair<String,String> >();
  renaming_to_standard_["HO4"] = std::vector<std::pair<String,String> >();
  renaming_to_standard_["HO5"] = std::vector<std::pair<String,String> >();
  renaming_to_standard_["HEME"] = std::vector<std::pair<String,String> >();
  renaming_to_standard_["protein"] = std::vector<std::pair<String,String> >();

  renaming_to_standard_["ILE"].push_back(std::make_pair("CD1","CD"));
  renaming_to_standard_["ILE"].push_back(std::make_pair("HD11","HD1"));
  renaming_to_standard_["ILE"].push_back(std::make_pair("HD12","HD2"));
  renaming_to_standard_["ILE"].push_back(std::make_pair("HD13","HD3"));

  renaming_to_standard_["HOH"].push_back(std::make_pair("O","OW"));
  renaming_to_standard_["HOH"].push_back(std::make_pair("OW1","OW"));
  renaming_to_standard_["HO4"].push_back(std::make_pair("O","OW"));
  renaming_to_standard_["HO4"].push_back(std::make_pair("OW1","OW"));
  renaming_to_standard_["HO5"].push_back(std::make_pair("O","OW"));
  renaming_to_standard_["HO5"].push_back(std::make_pair("OW1","OW"));

  renaming_to_standard_["HEME"].push_back(std::make_pair("N_A","NA"));
  renaming_to_standard_["HEME"].push_back(std::make_pair("N_B","NB"));
  renaming_to_standard_["HEME"].push_back(std::make_pair("N_C","NC"));
  renaming_to_standard_["HEME"].push_back(std::make_pair("N_D","ND"));

  renaming_to_standard_["protein"].push_back(std::make_pair("O1","O"));
  renaming_to_standard_["protein"].push_back(std::make_pair("O2","OXT"));
  renaming_to_standard_["protein"].push_back(std::make_pair("OT1","O"));
  renaming_to_standard_["protein"].push_back(std::make_pair("OT2","OXT"));
  renaming_to_standard_["protein"].push_back(std::make_pair("OT","OXT"));
  renaming_to_standard_["protein"].push_back(std::make_pair("O'","O"));
  renaming_to_standard_["protein"].push_back(std::make_pair("O''","OXT"));
  renaming_to_standard_["protein"].push_back(std::make_pair("OC1","O"));
  renaming_to_standard_["protein"].push_back(std::make_pair("OC2","OXT"));
  renaming_to_standard_["protein"].push_back(std::make_pair("HN","H"));
  renaming_to_standard_["protein"].push_back(std::make_pair("HT1","H1"));
  renaming_to_standard_["protein"].push_back(std::make_pair("HT2","H2"));
  renaming_to_standard_["protein"].push_back(std::make_pair("HT3","H3"));
}

int GromacsData::GetKeywordIndex(const String& keyword){
  String lowercase_keyword = keyword;
  boost::algorithm::to_lower(lowercase_keyword);
  if(keyword_map_.find(lowercase_keyword) == keyword_map_.end()) return -1;
  else return keyword_map_[lowercase_keyword];
}

String GromacsData::ConvertToStandard(const String& res_name, const String& atom_name){
  
  if(renaming_to_standard_.find(res_name) == renaming_to_standard_.end()) return atom_name;
  for(std::vector<std::pair<String,String> >::iterator i = renaming_to_standard_[res_name].begin();
      i != renaming_to_standard_[res_name].end(); ++i){
    if(i->first == atom_name) return i->second;
  }
  return  atom_name;
}

bool GromacsData::ConversionExists(const String& res_name){
  return renaming_to_standard_.find(res_name) != renaming_to_standard_.end();
}

GromacsDataPtr GromacsData::Instance(){
  if(!instance_){
    instance_ = GromacsDataPtr(new GromacsData);
  }
  return instance_;
}

GromacsDataPtr GromacsData::instance_ = GromacsDataPtr();

CHARMMData::CHARMMData() {

  keyword_map_["atoms"] = 0;
  keyword_map_["bonds"] = 1;
  keyword_map_["angles"] = 2;
  keyword_map_["dihedrals"] = 3;
  keyword_map_["improper"] = 4;
  keyword_map_["cmap"] = 5;
  keyword_map_["nonbonded"] = 6;
  keyword_map_["hbond"] = 7;
  keyword_map_["end"] = 8;
  keyword_map_["mass"] = 9;
  keyword_map_["decl"] = 10;
  keyword_map_["defa"] = 11;
  keyword_map_["auto"] = 12;
  keyword_map_["resi"] = 13;
  keyword_map_["pres"] = 14;
  keyword_map_["atom"] = 15;
  keyword_map_["group"] = 16;
  keyword_map_["bond"] = 17;
  keyword_map_["angle"] = 18;
  keyword_map_["dihe"] = 19;
  keyword_map_["impr"] = 20;
  keyword_map_["donor"] = 21;
  keyword_map_["acceptor"] = 22;
  keyword_map_["ic"] = 23;
  keyword_map_["patching"] = 24;
  keyword_map_["print"] = 25;
  keyword_map_["double"] = 26;
}

int CHARMMData::GetKeywordIndex(const String& keyword){
  String lowercase_keyword = keyword;
  boost::algorithm::to_lower(lowercase_keyword);
  if(keyword_map_.find(lowercase_keyword) == keyword_map_.end()) return -1;
  else return keyword_map_[lowercase_keyword];
}


CHARMMDataPtr CHARMMData::Instance(){
  if(!instance_){
    instance_ = CHARMMDataPtr(new CHARMMData);
  }
  return instance_;
}

CHARMMDataPtr CHARMMData::instance_ = CHARMMDataPtr();

std::vector<std::vector<String> > MMPreprocessor::Process(const String& filename){

  std::vector<std::vector<String> > file_content = this->ReadFile(filename);
  std::vector<std::vector<String> > temp;
  std::vector<std::vector<String> > return_file;

  //we first resolve all include statements
  String include_file_name;
  for(uint i = 0; i<file_content.size(); ++i){
    if(file_content[i][0] ==  "#include"){
      if(file_content[i].size() != 2){
        std::stringstream ss;
        ss<<"Can only evaluate include statement in combination with a filename.";
        ss<<"Also more stuff in the same line blows the thing up...";
        throw ost::io::IOException(ss.str());
      }
      String include_filename = file_content[i][1].substr(1,file_content[i][1].size()-2);
      temp = this->ReadFile(include_filename);
      //directly delete the line with the
      //include statement and fill in the full content
      file_content.insert(file_content.erase(file_content.begin()+i), temp.begin(), temp.end());
    }
  }

  //let's resolve all definitions
  std::map<String,std::vector<String> >::iterator definitions_it;
  for(uint i=0; i<file_content.size(); ++i){
    //resolve definition
    if(file_content[i][0] == "#define"){
      if(file_content[i].size()==2){
        defines_.insert(file_content[i][1]);
        file_content.erase(file_content.begin()+i);
        if(i == file_content.size()) break;
        --i;
        continue;
      }
      definitions_[file_content[i][1]] = std::vector<String>(file_content[i].begin()+2,file_content[i].end());
      file_content.erase(file_content.begin()+i);
      if(i == file_content.size()) break;
      --i;
      continue;
    }
    //resolve ifdef
    if(file_content[i][0] == "#ifdef" || file_content[i][0] == "#ifndef"){
      this->ResolveIFDEF(file_content,i);
      --i;
      continue;
    }
    //check whether a word in current line corresponds to a found definition
    for(uint j = 0; j<file_content[i].size(); ++j){
      definitions_it = definitions_.find(file_content[i][j]);
      if(definitions_it!=definitions_.end()){
        file_content[i].insert(file_content[i].erase(file_content[i].begin()+j),definitions_it->second.begin(),definitions_it->second.end());
      }
    }
  }
  return file_content;
}

void MMPreprocessor::ResolveIFDEF(std::vector<std::vector<String> >& file_content, int line_counter){

  String ifdef_statement = file_content[line_counter][1];
  bool logical_outcome;
  if(file_content[line_counter][0] == "#ifdef"){
    logical_outcome = defines_.find(ifdef_statement) != defines_.end();
  }
  else if(file_content[line_counter][0] == "#ifndef"){
    logical_outcome = defines_.find(ifdef_statement) == defines_.end();
  }
  else{
    throw ost::Error("Invalid preprocessor statement observed in ResolveIFDEF function!");
  }

  bool has_else = true;
  uint else_position = line_counter+1;
  while(true){
    if(else_position >= file_content.size()){
      throw ost::io::IOException("Invalid ifdef statement encountered!");
    }
    if(file_content[else_position][0] == "#ifdef" || file_content[else_position][0] == "#ifndef"){
      this->ResolveIFDEF(file_content,else_position);
    }
    if(file_content[else_position][0] == "#else") break;
    if(file_content[else_position][0] == "#endif"){
      has_else = false;
      break;
    }
    ++else_position;
  }

  uint endif_position = else_position+1;
  while(true && has_else){
    if(endif_position >= file_content.size()){
      throw ost::io::IOException("Invalid ifdef statement encountered!");
    }
    if(file_content[endif_position][0] == "#ifdef" || file_content[endif_position][0] == "#ifndef"){
      this->ResolveIFDEF(file_content,endif_position);
    }
    if(file_content[endif_position][0] == "#endif") break;
    ++endif_position;
  }

  if(has_else){
    if(logical_outcome){
      //delete part after else
      file_content.erase(file_content.begin()+else_position, file_content.begin()+endif_position+1);
      file_content.erase(file_content.begin()+line_counter);
    }
    else{
      //delete part before else
      file_content.erase(file_content.begin()+line_counter,file_content.begin()+else_position+1);
      file_content.erase(file_content.begin()+line_counter+(endif_position-else_position)-1);
    }
  }
  else{
    if(logical_outcome){//only delete ifdef and endif 
      file_content.erase(file_content.begin()+else_position);
      file_content.erase(file_content.begin()+line_counter);
    }
    else{//it's not defined => delete whole block
      file_content.erase(file_content.begin()+line_counter,file_content.begin()+else_position+1);
    }
  }
}

std::vector<std::vector<String> > MMPreprocessor::ReadFile(const String& filename){

  boost::filesystem::path file_path = basepath_ / filename;
  std::ifstream file;
  file.open((file_path.string()).c_str());

  if(!file){
    std::stringstream ss;
    ss<<"Could not read file. File '";
    ss<<filename<<"' does not exist!";
    throw ost::io::IOException(ss.str());
  }

  std::vector<std::vector<String> > split_file_content;
  String string_line;
  ost::StringRef stringref_line;
  std::vector<ost::StringRef> data;
  std::vector<String> fill_data;

  while(file.good()){
    std::getline(file,string_line);
    stringref_line = ost::StringRef(string_line.c_str(), string_line.length());
    stringref_line = stringref_line.trim();
    if(stringref_line.size() != 0 && stringref_line[0] != ';') {
      // Everything after ; is considered to be comment, so we
      // already filter that out
      data = stringref_line.split(';');
      if(!data.empty()){
        // lets split the remaining stuff
        data = data[0].split();
        fill_data.clear();
        for(std::vector<ost::StringRef>::iterator it = data.begin();
            it!=data.end(); ++it){
          fill_data.push_back(it->str());
        } 

        split_file_content.push_back(fill_data);
      }
    }
  }

  //handling linebreaks ("\")
  for(uint i = 0; i != split_file_content.size(); ++i){
    if(*split_file_content[i].back().rbegin() == '\\'){
      if(i < split_file_content.size() - 1){
        if(split_file_content[i].back() == "\\"){
          split_file_content[i].pop_back();
        }
        else{
          split_file_content[i].back() = split_file_content[i].back().substr(0,split_file_content[i].back().size()-1);
        }
        if(split_file_content[i].empty()){
          throw io::IOException("Observed single linebreak (\\) in one line of input file!");
        }
        for(std::vector<String>::iterator j = split_file_content[i+1].begin();
            j != split_file_content[i+1].end(); ++j){
          split_file_content[i].push_back(*j);
        }
        split_file_content.erase(split_file_content.begin()+i+1);
        --i;
      }
    }
  }
  
  return split_file_content;
}

FFReader::FFReader(const String& base_dir): preprocessor_(base_dir),
                                                      ff_(new Forcefield)
{ 
  for(int i = 0; i < 4; ++i){
    ff_bonded_types_.push_back(Invalid);
    bonded_types_.push_back(Invalid);
  }
}

void FFReader::ReadGromacsForcefield(){

  //read in the forcefield
  if(!boost::filesystem::exists(preprocessor_.GetBasedir() / "forcefield.itp")){
    throw io::IOException("directory must contain forcefield.itp file!");
  }
  if(!boost::filesystem::exists(preprocessor_.GetBasedir() / "atomtypes.atp")){
    throw io::IOException("directory must contain atomtypes.atp file!");
  }
  std::vector<std::vector<String> > content = preprocessor_.Process("forcefield.itp");
  this->ParseForcefield(content);
  content = preprocessor_.Process("atomtypes.atp");
  this->ParseAtomTypes(content);
}

void FFReader::ReadResidueDatabase(const String& basename){
  std::vector<std::vector<String> > data = preprocessor_.Process(basename+".rtp");
  this->ParseRTP(data);

  try{
    data = preprocessor_.Process(basename+".arn");
    this->ParseARN(data);
  }catch(ost::io::IOException& e) { }

  try{
    data = preprocessor_.Process(basename+".hdb");
    this->ParseHDB(data);
  }catch(ost::io::IOException& e) { }

  try{
    data = preprocessor_.Process(basename+".n.tdb");
    this->ParseNTDB(data);
  }catch(ost::io::IOException& e) { }

  try{
    data = preprocessor_.Process(basename+".c.tdb");
    this->ParseCTDB(data);
  }catch(ost::io::IOException& e) { }

  try{
    data = preprocessor_.Process(basename+".vsd");
    this->ParseVSD(data);
  }catch(ost::io::IOException& e) { }

  try{
    data = preprocessor_.Process(basename+".r2b");
    this->ParseRtoB(data);
  }catch(ost::io::IOException& e) { }
}

void FFReader::ReadCHARMMPRM(const String& filename){
  std::vector<std::vector<String> > data = preprocessor_.Process(filename);
  this->ParseCHARMMPRM(data);
}

void FFReader::ReadCHARMMRTF(const String& filename){
  std::vector<std::vector<String> > data = preprocessor_.Process(filename);
  this->ParseCHARMMRTF(data);
}

void FFReader::ParseForcefield(std::vector<std::vector<String> >& content){

  int keyword_id = -1;
  std::vector<String> current_line;

  for(uint line_counter = 0; line_counter<content.size(); ++line_counter){
    current_line = content[line_counter];
    if(current_line[0] == "["){
      keyword_id = GromacsData::Instance()->GetKeywordIndex(current_line[1]);
      if(keyword_id == -1){
        std::stringstream ss;
        ss << "Unknown keyword '"<<current_line[1];
        ss << "' in forcefield file."<<std::endl;
        throw io::IOException(ss.str());
      }
      continue;
    }

    switch(keyword_id){
      case 11:{
        InteractionPtr bond = this->ParseBond(current_line,true);
        ff_->AddBond(bond);
        ff_bonded_types_[0] = bond->GetFuncType();
        break;
      }
      case 12:{
        InteractionPtr constraint = this->ParseConstraint(current_line,true);
        ff_->AddConstraint(constraint);
        break;
      }
      case 13:{
        InteractionPtr angle = this->ParseAngle(current_line,true);
        ff_->AddAngle(angle);
        ff_bonded_types_[1] = angle->GetFuncType();
        break;
      }
      case 14:{

        InteractionPtr p = this->ParseDihedral(current_line,true);
        FuncType functype = p->GetFuncType();

        if(functype == PeriodicDihedral){
          ff_->AddDihedral(p);
          ff_bonded_types_[2] = functype;
        }
        else if(functype == PeriodicImproper || functype == HarmonicImproper){
          ff_->AddImproper(p);
          ff_bonded_types_[3] = functype;
        }
        break;
      }
      case 15:{
        InteractionPtr lj = this->ParseLJ(current_line,true);
        ff_->AddLJ(lj);
        break;
      }
      case 16:{
        InteractionPtr genborn = this->ParseGenborn(current_line,true);
        ff_->AddImplicitGenborn(genborn);
        break;
      }
      case 17:{
        InteractionPtr lj_pair = this->ParseLJPair(current_line,true);
        ff_->AddLJPair(lj_pair);
        break;
      }
      case 18:{
        std::stringstream ss;
        ss << "Observerd nonbond_params in forcefield file. Only forcefields with ";
        ss << "combining rules for sigma epsilon calculations are supported!"<<std::endl;
        throw ost::Error(ss.str());
        break;
      }
      case 19:{
        bool gen_pairs;
        Real fudge_lj, fudge_qq;
        if(current_line[2] == "yes") gen_pairs = true;
        else if(current_line[2] == "no") gen_pairs = false;
        else throw io::IOException("Only \"yes\" and \"no\" are supported for gen_pairs parameter!");
      
        fudge_lj = boost::lexical_cast<Real>(current_line[3]);
        fudge_qq = boost::lexical_cast<Real>(current_line[4]);

        ff_->SetGenPairs(gen_pairs);
        ff_->SetFudgeLJ(fudge_lj);
        ff_->SetFudgeQQ(fudge_qq);
        break;
      }
      case 20:{
        InteractionPtr cmap = this->ParseCMap(current_line,true);
        ff_->AddCMap(cmap);
        break;
      }
      default: break;
    }
  }
}

void FFReader::ParseAtomTypes(std::vector<std::vector<String> >& content){
  for(std::vector<std::vector<String> >::iterator i = content.begin();
      i!=content.end(); ++i){
    ff_->AddMass((*i)[0],boost::lexical_cast<Real>((*i)[1]));
  }
}

void FFReader::ParseRTP(std::vector<std::vector<String> >& content){
 
  read_residues_.clear();

  //we assume, that at the beginning is the keyword 'bondedtypes'
  //we neglect anything before that...
  uint line_counter = 0;
  for(; line_counter<content.size(); ++line_counter){
    if(content[line_counter].size()>1){
      if(content[line_counter][1] == "bondedtypes"){
        ++line_counter;
        //we're interested in the next line
        if(content[line_counter].size() < 4){
          throw ost::Error("Error in parsing bondedtypes section of RTP file!");
        }
        int gromacs_functype;
        //first element defines the bonds
        gromacs_functype = boost::lexical_cast<int>(content[line_counter][0]);
        if(gromacs_functype == 1){
          bonded_types_[0] = HarmonicBond;
        }
        else{
          std::stringstream ss;
          ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
          ss << " for bond in bondedtypes section of rtp file!";
          throw ost::Error(ss.str());
        }
        //second element defines angle
        gromacs_functype = boost::lexical_cast<int>(content[line_counter][1]);
        if(gromacs_functype == 5) bonded_types_[1] = UreyBradleyAngle;
        else if(gromacs_functype == 1) bonded_types_[1] = HarmonicAngle;
        else{
          std::stringstream ss;
          ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
          ss << " for angle in bondedtypes section of rtp file!";
          throw ost::Error(ss.str());
        }
        //third element defines the dihedral
        gromacs_functype = boost::lexical_cast<int>(content[line_counter][2]);
        if(gromacs_functype == 9 || gromacs_functype == 1) bonded_types_[2] = PeriodicDihedral;
        else{
          std::stringstream ss;
          ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
          ss << " for dihedral in bondedtypes section of rtp file!";
          throw ost::Error(ss.str());
        }
        //fourth elemend defines improper dihedral
        gromacs_functype = boost::lexical_cast<int>(content[line_counter][3]);
        if(gromacs_functype == 4) bonded_types_[3] = PeriodicImproper;
        else if(gromacs_functype == 2) bonded_types_[3] = HarmonicImproper;
        else{
          std::stringstream ss;
          ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
          ss << " for dihedral!";
          throw ost::Error(ss.str());
        }
        ++line_counter;
        break;
      }
    }
  }

  //data for specific residue will be filled in here and passed
  //in the residue_building block Parse function
  std::vector<std::vector<String> > residue_data;
  String residue_name;
  for(; line_counter < content.size(); ++line_counter){
    if(content[line_counter][0] == "["){
      //I there is a keyword (start with [) not defined above,
      //we assume, that it is a new residue
      if(GromacsData::Instance()->GetKeywordIndex(content[line_counter][1]) == -1){
        //parse the old residue and set new residue name
        if(!residue_data.empty()){
          ff_->AddBuildingBlock(residue_name, this->BlockFromRTP(residue_data));
          read_residues_.push_back(residue_name);
        }
        residue_data.clear();
        residue_name = content[line_counter][1];
        continue;
      }
    }
    residue_data.push_back(content[line_counter]);
  }
  //we probably missed the last residue
  if(!residue_data.empty()){
    ff_->AddBuildingBlock(residue_name, this->BlockFromRTP(residue_data));
    read_residues_.push_back(residue_name);
  }
  

  //It could be, that certain bondedtypes in the forcefield are not defined yet,
  //since all forces are defined in the building blocks...
  //We assume, that the whole stuff is consistent across all forcefield and rtp
  //files... If not, you'll get errors when setting up your system
  for(int i = 0; i < 4; ++i){
    if(ff_bonded_types_[i] == Invalid) ff_bonded_types_[i] = bonded_types_[i];
  }

}

void FFReader::ParseARN(std::vector<std::vector<String> >& content){
  String res_name, a_one, a_two;

  for(std::vector<std::vector<String> >::iterator i = content.begin();
      i!=content.end(); ++i){
    res_name = (*i)[0];
    a_one = (*i)[1];
    a_two = (*i)[2];
    if(res_name == "*"){
      // * is wildcard... valid for all residues in the previously read rtp file
      for(std::vector<String>::iterator j = read_residues_.begin();
          j != read_residues_.end(); ++j){
        ff_->AddAtomRenamingRule(*j,a_one,a_two);
      }
    }
    else{
      ff_->AddAtomRenamingRule(res_name, a_one, a_two); 
    }
  }
}

void FFReader::ParseHDB(std::vector<std::vector<String> >& content){

  std::vector<std::vector<String> > residue_data;
  String residue_name;
  for(uint line_counter = 0; line_counter<content.size(); ++line_counter){
    if(content[line_counter].size() == 2){
      if(!residue_data.empty()){
        GromacsHydrogenConstructor* p = new GromacsHydrogenConstructor;
        for(std::vector<std::vector<String> >::iterator i = residue_data.begin();
            i != residue_data.end(); ++i){
          this->ParseHydrogenRule(*i,*p);
        }
        ff_->AddHydrogenConstructor(residue_name,HydrogenConstructorPtr(p));
      }
      residue_data.clear();
      residue_name = content[line_counter][0];
      continue;
    }
    residue_data.push_back(content[line_counter]);
  }
  //we probably missed the last residue
  if(!residue_data.empty()){
    GromacsHydrogenConstructor* p = new GromacsHydrogenConstructor;
    for(std::vector<std::vector<String> >::iterator i = residue_data.begin();
        i != residue_data.end(); ++i){
      this->ParseHydrogenRule(*i,*p);
    }
    ff_->AddHydrogenConstructor(residue_name,HydrogenConstructorPtr(p));
  }
}

void FFReader::ParseHydrogenRule(const std::vector<String>& data, 
                                      GromacsHydrogenConstructor& constructor){

  if(data.size() < 4){
    throw ost::Error("Require at least four items when parsing hydrogen rule!");
  }

  int number = boost::lexical_cast<int>(data[0]);
  int method = boost::lexical_cast<int>(data[1]);
  String prefix = data[2];
  std::vector<String> anchors;
  std::vector<String> hydrogen_names;

  for(uint i=3; i<data.size(); ++i){
    anchors.push_back(data[i]);
  } 

  switch(number){
    case 1:{
      hydrogen_names.push_back(prefix);
      break;
    }
    case 2:{
      hydrogen_names.push_back(prefix+"1");
      hydrogen_names.push_back(prefix+"2");
      break;
    }
    case 3:{
      hydrogen_names.push_back(prefix+"1");
      hydrogen_names.push_back(prefix+"2");
      hydrogen_names.push_back(prefix+"3");
      break;
    }
    default:{
      break;
    }
  }

  constructor.AddHydrogenRule(number,method,hydrogen_names, anchors);
}


void FFReader::ParseTerminiReplaceRule(const std::vector<String>& data, 
                                            GromacsBlockModifier& modifier){

  String name;
  String new_name;
  String new_type;
  Real new_charge;
  if(data.size() == 5 && data.back() != "0"){
    name = data[0];
    new_name = data[1];
    new_type = data[2];
    new_charge = boost::lexical_cast<Real>(data[4]);
  }
  else{
    name = data[0];
    new_name = name;
    new_type = data[1];
    new_charge = boost::lexical_cast<Real>(data[3]);
  }
 
  modifier.AddReplaceRule(name, new_name, new_type, new_charge);
  //we neglect the mass! This will be read from the forcefield itself!
}

void FFReader::ParseTerminiAddRule(const std::vector<String>& data1, 
                                        const std::vector<String>& data2,
                                        GromacsBlockModifier& modifier){

  int number = boost::lexical_cast<int>(data1[0]);
  int method = boost::lexical_cast<int>(data1[1]);
  String prefix = data1[2];

  std::vector<String> anchors;
  for(uint i = 3; i < data1.size(); ++i){
    anchors.push_back(data1[i]);
  }

  String type = data2[0];
  Real charge = boost::lexical_cast<Real>(data2[2]);
  //we again neglect the masses...

  std::vector<String> names;

  switch(number){
    case 1:{
      names.push_back(prefix);
      break;
    }
    case 2:{
      names.push_back(prefix+"1");
      names.push_back(prefix+"2");
      break;
    }
    case 3:{
      names.push_back(prefix+"1");
      names.push_back(prefix+"2");
      names.push_back(prefix+"3");
      break;
    }
    default:{
      std::stringstream ss;
      ss<<"Can only add one, two or three atoms at one position.";
      ss<<"The provided force field file seems to tell something different";
      throw ost::Error(ss.str());
    }
  }

  modifier.AddAddRule(number, method, names, anchors, type, charge);

}


void FFReader::ParseNTDB(std::vector<std::vector<String> >& content){

  //Besided general termini constructors, there are also residue 
  //specific constructors => custom constructors
  std::vector<String> termini;
  std::map<String,std::vector<String> > custom_termini;

  std::vector<std::vector<String> > data;
  String name;

  for(uint line_counter = 0; line_counter < content.size(); ++line_counter){
    if(content[line_counter][0] == "["){
      //I there is a keyword (start with [) not defined above,
      //we assume, that it is a new residue
      if(GromacsData::Instance()->GetKeywordIndex(content[line_counter][1]) == -1){
        //parse the old residue and set new residue name
        if(!data.empty()){
          BlockModifierPtr p = this->ParseBlockModifier(data);
          ff_->AddBlockModifier(name,p);
          size_t minus_pos = name.find("-");
          if(minus_pos != String::npos && minus_pos < (name.size()-1)){
            String a = name.substr(0,minus_pos);
            if(custom_termini.find(a) == custom_termini.end()){
              custom_termini[a] = std::vector<String>();
            }
            custom_termini[a].push_back(name); 
          }
          else{
            termini.push_back(name);
          }
        }
        data.clear();
        name = content[line_counter][1];
        continue;
      }
    }
    data.push_back(content[line_counter]);
  }

  if(!data.empty()){
    BlockModifierPtr p = this->ParseBlockModifier(data);
    ff_->AddBlockModifier(name,p);
    size_t minus_pos = name.find("-");
    if(minus_pos != String::npos && minus_pos < (name.size()-1)){
      String a = name.substr(0,minus_pos);
      if(custom_termini.find(a) == custom_termini.end()){
        custom_termini[a] = std::vector<String>();
      }
      custom_termini[a].push_back(name);
    }
    else{
      termini.push_back(name);
    }
  }

  for(std::vector<String>::iterator i = read_residues_.begin(); 
      i != read_residues_.end(); ++i){
    if(custom_termini.find(*i) != custom_termini.end()){
      if(!custom_termini[*i].empty()) ff_->SetStandardNTer(*i,custom_termini[*i][0]);
      continue;
    }
    if(!termini.empty()) ff_->SetStandardNTer(*i,termini[0]);
  }
}

void FFReader::ParseCTDB(std::vector<std::vector<String> >& content){

  //Besided general termini constructors, there are also residue 
  //specific constructors => custom constructors
  std::vector<String> termini;
  std::map<String,std::vector<String> > custom_termini;

  std::vector<std::vector<String> > data;
  String name;

  for(uint line_counter = 0; line_counter < content.size(); ++line_counter){
    if(content[line_counter][0] == "["){
      //I there is a keyword (start with [) not defined above,
      //we assume, that it is a new residue
      if(GromacsData::Instance()->GetKeywordIndex(content[line_counter][1]) == -1){
        //parse the old residue and set new residue name
        if(!data.empty()){
          BlockModifierPtr p = this->ParseBlockModifier(data);
          ff_->AddBlockModifier(name,p);
          size_t minus_pos = name.find("-");
          if(minus_pos != String::npos && minus_pos < (name.size()-1)){
            String a = name.substr(0,minus_pos);
            if(custom_termini.find(a) == custom_termini.end()){
              custom_termini[a] = std::vector<String>();
            }
            custom_termini[a].push_back(name); 
          }
          else{
            termini.push_back(name);
          }
        }
        data.clear();
        name = content[line_counter][1];
        continue;
      }
    }
    data.push_back(content[line_counter]);
  }

  if(!data.empty()){
    BlockModifierPtr p = this->ParseBlockModifier(data);
    ff_->AddBlockModifier(name,p);
    size_t minus_pos = name.find("-");
    if(minus_pos != String::npos && minus_pos < (name.size()-1)){
      String a = name.substr(0,minus_pos);
      if(custom_termini.find(a) == custom_termini.end()){
        custom_termini[a] = std::vector<String>();
      }
      custom_termini[a].push_back(name);
    }
    else{
      termini.push_back(name);
    }
  }

  for(std::vector<String>::iterator i = read_residues_.begin(); 
      i != read_residues_.end(); ++i){
    if(custom_termini.find(*i) != custom_termini.end()){
      if(!custom_termini[*i].empty()) ff_->SetStandardCTer(*i,custom_termini[*i][0]);
      continue;
    }
    if(!termini.empty()) ff_->SetStandardCTer(*i,termini[0]);
  }
}

void FFReader::ParseVSD(std::vector<std::vector<String> >& content){

}

void FFReader::ParseRtoB(std::vector<std::vector<String> >& content){

  //there are two versions of this file, either with two or five
  //columns further instructions can be found in the gromacs manual

  String res_name, main, n_ter, c_ter, two_ter;

  for(std::vector<std::vector<String> >::iterator i = content.begin();
      i != content.end(); ++i){
    if((*i).size() == 2){
      res_name = (*i)[0];
      main = (*i)[1];
      n_ter = main;
      c_ter = main;
      two_ter = main;
      ff_->AddResidueRenamingRule(res_name,main,n_ter,c_ter,two_ter);
      continue;
    }
    if((*i).size() == 5){
      res_name = (*i)[0];
      main = (*i)[1];
      n_ter = (*i)[2];
      c_ter = (*i)[3];
      two_ter = (*i)[4];
      ff_->AddResidueRenamingRule(res_name,main,n_ter,c_ter,two_ter);
      continue;
    }
    throw io::IOException("Entries in r2b files must contain 2 or 5 columns");
  }
}

void FFReader::ReadITP(const String& basename){
  std::vector<std::vector<String> > content = preprocessor_.Process(basename + ".itp");
  this->ParseITP(content);
}

void FFReader::ParseITP(std::vector<std::vector<String> >& content){

  std::vector<std::vector<String> > residue_data;
  String residue_name;
  for(uint line_counter = 0; line_counter < content.size(); ++line_counter){
    if(content[line_counter][0] == "["){
      if(GromacsData::Instance()->GetKeywordIndex(content[line_counter][1]) == 8){
        //parse the old residue and set new residue name
        if(!residue_data.empty()){
          BuildingBlockPtr p = this->BlockFromITP(residue_data);
          ff_->AddBuildingBlock(residue_name, p);
        }
        residue_data.clear();
        residue_name = content[++line_counter][0];
        continue;
      }
    }
    residue_data.push_back(content[line_counter]);
  }
  //we probably missed the last residue
  if(!residue_data.empty()){
    BuildingBlockPtr p = this->BlockFromITP(residue_data);
    ff_->AddBuildingBlock(residue_name, p);
  }
}

BuildingBlockPtr FFReader::BlockFromRTP(const std::vector<std::vector<String> >& data){  

  if(bonded_types_.size()<4){
    std::stringstream ss;
    ss << "There must be at least 4 parameters in the bondedtypes section of";
    ss << ".rtp file";
    throw io::IOException(ss.str());
  }

  if(data[0][1] != "atoms"){
    std::stringstream ss;
    ss << "Expect residue definition to start with the atom keyword!"; 
    throw io::IOException(ss.str());
  }

  int keyword_index = 0;
  BuildingBlockPtr p(new BuildingBlock);

  for(uint i=1;i<data.size(); ++i){
    if(data[i].size()>1){
      if(data[i][0]=="["){
        keyword_index = GromacsData::Instance()->GetKeywordIndex(data[i][1]);
        continue;
      }
    }

    switch(keyword_index){
      case 0:{
        p->AddAtom(data[i][0],
                   data[i][1],
                   boost::lexical_cast<Real>(data[i][2]));
        break;
      }
      case 1:{
        p->AddBond(this->ParseBond(data[i],false,bonded_types_[0]));
        break;
      }
      case 2:{
        p->AddAngle(this->ParseAngle(data[i],false,bonded_types_[1]));
        break;
      }
      case 3:{
        p->AddImproper(this->ParseDihedral(data[i],false,bonded_types_[3]));
        break;
      }
      case 4:{
        p->AddDihedral(this->ParseDihedral(data[i],false,bonded_types_[2]));
        break;
      }
      case 5:{
        p->AddExclusion(this->ParseExclusion(data[i],false)); 
        break;
      }
      case 6:{
        p->AddCMap(this->ParseCMap(data[i],false));
        break;
      }
      default:{
        std::stringstream ss;
        ss << "Encountered invalid keyword while parsing rtp file: " << data[i-1][1];
        throw ost::io::IOException(ss.str());
      }
    }
  }
  return p;
}

BuildingBlockPtr FFReader::BlockFromITP(const std::vector<std::vector<String> >& data){

  int keyword_index;

  if(data[0][1] != "atoms"){
    std::stringstream ss;
    ss << "Expect residue definition to start with the atom keyword!"; 
    throw io::IOException(ss.str());
  }

  std::map<int,String> index_name_mapper;
  keyword_index = 0;
  BuildingBlockPtr p(new BuildingBlock);
  std::vector<String> current_line;

  for(uint i=1;i<data.size(); ++i){
    if(data[i].size()>1){
      if(data[i][0]=="["){
        keyword_index = GromacsData::Instance()->GetKeywordIndex(data[i][1]);
        continue;
      }
    }
    current_line = data[i];
    switch(keyword_index){
      case 0:{
        p->AddAtom(current_line[4],
                   current_line[1],
                   boost::lexical_cast<Real>(current_line[6]));
        index_name_mapper[boost::lexical_cast<int>(current_line[0])] = current_line[4];
        break;
      }
      case 1:{
        current_line[0] = index_name_mapper[boost::lexical_cast<int>(current_line[0])];
        current_line[1] = index_name_mapper[boost::lexical_cast<int>(current_line[1])];
        p->AddBond(this->ParseBond(current_line,false));
        break;
      }
      case 2:{
        current_line[0] = index_name_mapper[boost::lexical_cast<int>(current_line[0])];
        current_line[1] = index_name_mapper[boost::lexical_cast<int>(current_line[1])];       
        current_line[2] = index_name_mapper[boost::lexical_cast<int>(current_line[2])];     
        p->AddAngle(this->ParseAngle(current_line,false));
        break;
      }
      case 4:{
        //there must be a function type, otherwise we can't distinguish
        //between dihedrals and impropers
        current_line[0] = index_name_mapper[boost::lexical_cast<int>(current_line[0])];
        current_line[1] = index_name_mapper[boost::lexical_cast<int>(current_line[1])];       
        current_line[2] = index_name_mapper[boost::lexical_cast<int>(current_line[2])];
        current_line[3] = index_name_mapper[boost::lexical_cast<int>(current_line[3])];

        InteractionPtr int_ptr = this->ParseDihedral(current_line,false);
        FuncType functype = int_ptr->GetFuncType();

        if(functype == PeriodicDihedral){
          p->AddDihedral(int_ptr);
        }
        else if(functype == PeriodicImproper || functype == HarmonicImproper){
          p->AddImproper(int_ptr);
        }
        break;
      }
      case 5:{
        //current_line[0] = index_name_mapper[boost::lexical_cast<int>(current_line[0])];
        //current_line[1] = index_name_mapper[boost::lexical_cast<int>(current_line[1])];
        //p->exclusions_.push_back(Interaction::ParseExclusion(current_line,false)); 
        break;
      }
      case 6:{
        current_line[0] = index_name_mapper[boost::lexical_cast<int>(current_line[0])];
        current_line[1] = index_name_mapper[boost::lexical_cast<int>(current_line[1])];       
        current_line[2] = index_name_mapper[boost::lexical_cast<int>(current_line[2])];
        current_line[3] = index_name_mapper[boost::lexical_cast<int>(current_line[3])];
        current_line[4] = index_name_mapper[boost::lexical_cast<int>(current_line[4])];
        p->AddCMap(this->ParseCMap(current_line,false));
        break;
      }
      case 21:{
        //we extract pairs manually...
        break;
      }
      case 22:{
        std::vector<String> block_atom_names = p->GetAtoms();
        if(block_atom_names.size() < 3){
          throw ost::Error("Molecule with settles keyword must contain at least 3 atoms!");
        }
        std::vector<Real> parameters;
        std::vector<String> names;


        InteractionPtr p1(new Interaction(DistanceConstraint));
        InteractionPtr p2(new Interaction(DistanceConstraint));
        InteractionPtr p3(new Interaction(DistanceConstraint));

        //first constraint dist => O - H
        parameters.push_back(boost::lexical_cast<Real>(current_line[2]));
        names.push_back(block_atom_names[0]);
        names.push_back(block_atom_names[1]);
        p1->SetNames(names);
        p1->SetParam(parameters);
        p->AddConstraint(p1);

        //second constraint dist => O - H
        names[1] = block_atom_names[2];
        p2->SetNames(names);
        p2->SetParam(parameters);
        p->AddConstraint(p2);

        //third constraint dist => H - H
        names[0] = block_atom_names[1];
        parameters[0] = boost::lexical_cast<Real>(current_line[3]);
        p3->SetNames(names);
        p3->SetParam(parameters);
        p->AddConstraint(p3);

        break;
      }
      default:{
        std::stringstream ss;
        ss << "Encountered invalid keyword while parsing itp file: " << data[i-1][1];
        throw ost::io::IOException(ss.str());
      }
    }
  } 
  return p;
}


BlockModifierPtr FFReader::ParseBlockModifier(const std::vector<std::vector<String> >& data){


  int keyword_index = -1;

  if(data[0].size() > 0){
    keyword_index = GromacsData::Instance()->GetKeywordIndex(data[0][1]);
    if(keyword_index == -1){
      throw ost::io::IOException("Expect termini data to start with a valid keyword!");
    }
  }

  GromacsBlockModifier* p = new GromacsBlockModifier; 

  for(uint i = 1; i < data.size(); ++i){
    if(data[i].size()>1){
      if(data[i][0] == "["){
        keyword_index = GromacsData::Instance()->GetKeywordIndex(data[i][1]);
        continue;
      }
    }
    switch(keyword_index){
      case 7:{
        this->ParseTerminiReplaceRule(data[i],*p);
        break;
      }
      case 9:{
        std::vector<String> line1 = data[i];
        std::vector<String> line2 = data[++i];
        this->ParseTerminiAddRule(line1, line2, *p);
        break;
      }
      case 10:{
        p->AddDeleteAtom(data[i][0]);
        break;
      }
      case 1:{
        p->AddBond(this->ParseBond(data[i],false,bonded_types_[0]));
        break;
      }
      case 2:{
        p->AddAngle(this->ParseAngle(data[i],false,bonded_types_[1]));
        break;
      }
      case 4:{
        p->AddDihedral(this->ParseDihedral(data[i],false,bonded_types_[2]));
        break;
      }
      case 3:{
        p->AddImproper(this->ParseDihedral(data[i],false,bonded_types_[3]));
        break;
      }
      case 6:{
        p->AddCMap(this->ParseCMap(data[i],false));
        break;
      }

      default:{
        std::stringstream ss;
        ss << "Encountered invalid keyword while parsing termini: " << data[i-1][1];
        throw ost::io::IOException(ss.str());
      }
    }
  }

  return BlockModifierPtr(p);
}


InteractionPtr FFReader::ParseBond(const std::vector<String>& data, 
                                                 bool type_definition,
                                                 FuncType functype){

  if(data.size() < 2){
    throw ost::io::IOException("Require at least 2 arguments in parsing bond!");
  }

  if(functype == Invalid){
    if(data.size() == 2){
      throw ost::Error("Error in determining bond function type!");
    }
    int gromacs_functype;
    try{ gromacs_functype = boost::lexical_cast<int>(data[2]); } catch(std::exception& e){
      throw ost::io::IOException("Could not cast function type into int while parsing bond!");
    }
    if(gromacs_functype == 1){
      functype = HarmonicBond;
    }
    else{
      std::stringstream ss;
      ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
      ss << " for bond!";
      throw ost::Error(ss.str());
    }
  }

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);
  atom_definition.push_back(data[1]);

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  try{
    if(data.size() > 3){
      std::vector<Real> parameters;
      parameters.push_back(boost::lexical_cast<Real>(data[data.size()-2]));         
      parameters.push_back(boost::lexical_cast<Real>(data[data.size()-1]));
     
      p->SetParam(parameters);
    }
  }catch(std::exception& e){ 
    throw ost::io::IOException("Could not cast bond parameters into real type!");
  }
  return p;
}

InteractionPtr FFReader::ParseAngle(const std::vector<String>& data, 
                                           bool type_definition,
                                           FuncType functype){

  if(data.size() < 3){
    throw ost::io::IOException("Require at least 3 arguments in parsing angle!");
  }

  if(functype == Invalid){
    if(data.size() == 3){
      throw ost::Error("Error in determining angle function type!");
    }
    int gromacs_functype;
    try{ gromacs_functype = boost::lexical_cast<int>(data[3]); } catch(std::exception& e){
      throw ost::io::IOException("Could not cast function type into int while parsing angle!");
    }
    if(gromacs_functype == 5) functype = UreyBradleyAngle;
    else if(gromacs_functype == 1) functype = HarmonicAngle;
    else{
      std::stringstream ss;
      ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
      ss << " for angle!";
      throw ost::Error(ss.str());
    }
  }

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);
  atom_definition.push_back(data[1]);
  atom_definition.push_back(data[2]);

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  try{
    if(data.size() > 4){
      std::vector<Real> parameters;
      if(functype == UreyBradleyAngle){
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-4])/360*2*M_PI);         
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-3]));
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-2]));         
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-1]));

      }
      else{
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-2])/360*2*M_PI);         
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-1]));
      }
      p->SetParam(parameters);
    }
  }catch(std::exception& e){
      throw ost::io::IOException("Could not cast angle parameters into real type!");
  }
  return p;
}

InteractionPtr FFReader::ParseDihedral(const std::vector<String>& data, 
                                              bool type_definition,
                                              FuncType functype){
  if(data.size() < 4){
    throw ost::io::IOException("Require at least 4 arguments in parsing dihedral/improper!");
  }

  if(functype == Invalid){
    if(data.size() == 4){
      throw ost::Error("Error in determining angle function type!");
    }
    int gromacs_functype;
    try{ gromacs_functype = boost::lexical_cast<int>(data[4]); } catch(std::exception& e){
      throw ost::io::IOException("Could not cast function type into int while parsing dihedral!");
    }
    if(gromacs_functype == 9 || gromacs_functype == 1) functype = PeriodicDihedral;
    else if(gromacs_functype == 4) functype = PeriodicImproper;
    else if(gromacs_functype == 2) functype = HarmonicImproper;
    else{
      std::stringstream ss;
      ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
      ss << " for dihedral!";
      throw ost::Error(ss.str());
    }
  }

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);
  atom_definition.push_back(data[1]);
  atom_definition.push_back(data[2]);
  atom_definition.push_back(data[3]);

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  try{
    if(data.size() > 5){
      std::vector<Real> parameters;
      if(functype == HarmonicImproper){
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-2])/360*2*M_PI);         
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-1]));
      }
      else{
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-1]));
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-3])/360*2*M_PI);
        parameters.push_back(boost::lexical_cast<Real>(data[data.size()-2]));         
      }
      p->SetParam(parameters);
    }
  }catch(std::exception& e){
      throw ost::io::IOException("Could not cast dihedral/improper parameters into real type!");
  }
  return p;
}

InteractionPtr FFReader::ParseCMap(const std::vector<String>& data, 
                                                   bool type_definition,
                                                   FuncType functype){
  if(data.size() < 5){
    throw ost::io::IOException("Require at least 5 arguments in parsing cmap!");
  }

  if(functype == Invalid){
    functype = CMap;  //well, there is only this possibility
  }

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);
  atom_definition.push_back(data[1]);
  atom_definition.push_back(data[2]);
  atom_definition.push_back(data[3]);
  atom_definition.push_back(data[4]);

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  if(data.size() > 6){
    uint a,b;
    try{
      a = boost::lexical_cast<uint>(data[6]);
      b = boost::lexical_cast<uint>(data[7]);
    }catch(std::exception& e){
      throw ost::io::IOException("Failed to parse size of the cmap!");
    }
    if(a != b){
      throw ost::io::IOException("Expect cmap to have equal extension into both dimensions!");
    }

    if((8+a*b) > data.size()){
      throw ost::io::IOException("Expected more data to fill cmap!");
    }
    try{
      std::vector<Real> parameters;
      parameters.push_back(Real(a));
      for(uint i = 0; i < a*a; ++i){
        parameters.push_back(boost::lexical_cast<Real>(data[8+i]));
      }
      p->SetParam(parameters);
    }catch(std::exception& e){
      throw ost::io::IOException("Failed to cast cmap values into real type!");
    }
  }
  return p;
}


InteractionPtr FFReader::ParseLJ(const std::vector<String>& data, 
                                        bool type_definition,
                                        FuncType functype){

  if(data.size() < 1){
    throw ost::io::IOException("Require at least 1 arguments in parsing LJ!");
  }

  if(functype == Invalid){
    functype = LJ;  //well, there is only this possibility
  }

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  try{
    if(data.size() > 2){
      std::vector<Real> parameters;
      parameters.push_back(boost::lexical_cast<Real>(data[data.size()-2])); 
      parameters.push_back(boost::lexical_cast<Real>(data[data.size()-1]));
      p->SetParam(parameters);
    }
  }catch(std::exception& e){ 
    throw ost::io::IOException("Could not cast atomtype parameters into real type!");
  }
  return p;
}


InteractionPtr FFReader::ParseLJPair(const std::vector<String>& data, 
                                            bool type_definition,
                                            FuncType functype){

  if(data.size() < 2){
    throw ost::io::IOException("Require at least 2 arguments in parsing nonbonded parameters!");
  }

  if(functype == Invalid){
    if(data.size() == 2){
      throw ost::Error("Error in determining LJ Pair function type!");
    }
    int gromacs_functype;
    try{ gromacs_functype = boost::lexical_cast<int>(data[2]); } catch(std::exception& e){
      throw ost::io::IOException("Could not cast function type into int while parsing lj pair!");
    }
    if(gromacs_functype == 1){
      functype = LJPair;
    }
    else{
      std::stringstream ss;
      ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
      ss << " for LJ Pair!";
      throw ost::Error(ss.str());
    }
  }

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);
  atom_definition.push_back(data[1]);

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  try{
    if(data.size() > 3){
      std::vector<Real> parameters;
      parameters.push_back(boost::lexical_cast<Real>(data[data.size()-2]));         
      parameters.push_back(boost::lexical_cast<Real>(data[data.size()-1]));
      p->SetParam(parameters);
    }
  }catch(std::exception& e){ 
    throw ost::io::IOException("Could not cast nonbonded parameters into real type!");
  }
  return p;
}

InteractionPtr FFReader::ParseConstraint(const std::vector<String>& data, 
                                                bool type_definition,
                                                FuncType functype){

  if(data.size() < 2){
    throw ost::io::IOException("Require at least 2 arguments in parsing nonbonded parameters!");
  }

  if(functype == Invalid){
    if(data.size() == 2){
      throw ost::Error("Error in determining constraint function type!");
    }
    int gromacs_functype;
    try{ gromacs_functype = boost::lexical_cast<int>(data[2]); } catch(std::exception& e){
      throw ost::io::IOException("Could not cast function type into int while parsing constraint!");
    }
    if(gromacs_functype == 1 || gromacs_functype == 2){
      functype = DistanceConstraint;
    }
    else{
      std::stringstream ss;
      ss << "Not supported gromacs function type \""<<gromacs_functype<<"\"";
      ss << " for constraint!";
      throw ost::Error(ss.str());
    }
  }

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);
  atom_definition.push_back(data[1]);

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  try{
    if(data.size() > 3){
      std::vector<Real> parameters;
      parameters.push_back(boost::lexical_cast<Real>(data[data.size()-1]));
      p->SetParam(parameters);
    }
  }catch(std::exception& e){ 
    throw ost::io::IOException("Could not cast constraint parameters into real type!");
  }
  return p;
}

InteractionPtr FFReader::ParseGenborn(const std::vector<String>& data,
                                             bool type_definition,
                                             FuncType functype){


  if(data.size() < 1){
    throw ost::io::IOException("Require at least 1 argument in parsing genborn parameters!");
  }

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);

  functype = GBSA; //bit hacky, I know, but there is only one type for gbsa anyway....

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  try{
    if(data.size() > 1){
      std::vector<Real> parameters;
      parameters.push_back(boost::lexical_cast<Real>(data[4]));
      parameters.push_back(boost::lexical_cast<Real>(data[5]));
      p->SetParam(parameters);
    }
  }catch(std::exception& e){ 
    throw ost::io::IOException("Could not cast genborn parameters into real type!");
  }
  return p;
}


InteractionPtr FFReader::ParseExclusion(const std::vector<String>& data,
                                               bool type_definition,
                                               FuncType functype){

  if(data.size() < 2){
    throw ost::io::IOException("Require at least 2 argument in parsing exclusion parameters!");
  }

  functype = Exclusion; //bit hacky, I know, but there is only one type for exclusion anyway....

  std::vector<String> atom_definition;
  atom_definition.push_back(data[0]);
  atom_definition.push_back(data[1]);

  InteractionPtr p(new Interaction(functype));

  if(type_definition) p->SetTypes(atom_definition);
  else p->SetNames(atom_definition);

  return p;
}

void FFReader::ParseCHARMMPRM(std::vector<std::vector<String> >& content){


  int keyword_id = -1;
  std::vector<String> current_line;

  //will be used to calculate all pairwise 1,4 lj parameters later on
  std::vector<String> lj_14_types;
  std::vector<Real> lj_14_sigmas;
  std::vector<Real> lj_14_epsilons;
  std::vector<String> lj_types;
  std::vector<Real> lj_sigmas;
  std::vector<Real> lj_epsilons;

  for(uint line_counter = 0; line_counter<content.size(); ++line_counter){
    current_line = content[line_counter];
    if(CHARMMData::Instance()->GetKeywordIndex(current_line[0]) != -1){
      keyword_id = CHARMMData::Instance()->GetKeywordIndex(current_line[0]);
      continue;
    }

    if(keyword_id == 8) break;

    switch(keyword_id){
      case 0:{ //atoms keyword
        ff_->AddMass(current_line[2],boost::lexical_cast<Real>(current_line[3]));
        break;
      }
      case 1:{ //bonds
        InteractionPtr p(new Interaction(HarmonicBond));
        std::vector<String> types;
        std::vector<Real> parameters;
        types.push_back(current_line[0]);
        types.push_back(current_line[1]);
        try{
          parameters.push_back(0.1 * boost::lexical_cast<Real>(current_line[3]));
          parameters.push_back(418.4*2 * boost::lexical_cast<Real>(current_line[2]));
        }catch(std::exception& e){
          throw ost::io::IOException("Could not cast bond parameters into real type!");
        }
        p->SetTypes(types);
        p->SetParam(parameters);
        ff_->AddBond(p);
        break;
      }
      case 2:{ //angles
        InteractionPtr p(new Interaction(UreyBradleyAngle));
        std::vector<String> types;
        std::vector<Real> parameters;
        types.push_back(current_line[0]);
        types.push_back(current_line[1]);
        types.push_back(current_line[2]);
        try{
          parameters.push_back(boost::lexical_cast<Real>(current_line[4])/360*2*M_PI);         
          parameters.push_back(boost::lexical_cast<Real>(current_line[3])*4.184*2);
          if(current_line.size() > 5){
            parameters.push_back(boost::lexical_cast<Real>(current_line[6])*0.1);
            parameters.push_back(boost::lexical_cast<Real>(current_line[5])*418.4*2);
          }
          else{
            parameters.push_back(0.0);
            parameters.push_back(0.0);
          }
        }catch(std::exception& e){
          throw ost::io::IOException("Could not cast angle parameters into real type!");
        }
        p->SetTypes(types);
        p->SetParam(parameters);
        ff_->AddAngle(p);
        break;
      }
      case 3:{ //diherals
        InteractionPtr p(new Interaction(PeriodicDihedral));
        std::vector<String> types;
        std::vector<Real> parameters;
        types.push_back(current_line[0]);
        types.push_back(current_line[1]);
        types.push_back(current_line[2]);
        types.push_back(current_line[3]);
        try{
          parameters.push_back(boost::lexical_cast<Real>(current_line[5]));
          parameters.push_back(boost::lexical_cast<Real>(current_line[6])/360*2*M_PI);
          parameters.push_back(boost::lexical_cast<Real>(current_line[4])*4.184);
        }catch(std::exception& e){
          throw ost::io::IOException("Could not cast bond parameters into real type!");
        }
        p->SetTypes(types);
        p->SetParam(parameters);
        ff_->AddDihedral(p);
        break;
      }
      case 4:{ //improper
        InteractionPtr p(new Interaction(HarmonicImproper));
        std::vector<String> types;
        std::vector<Real> parameters;
        types.push_back(current_line[0]);
        types.push_back(current_line[1]);
        types.push_back(current_line[2]);
        types.push_back(current_line[3]);
        try{
          parameters.push_back(boost::lexical_cast<Real>(current_line[6])/360*2*M_PI);
          parameters.push_back(boost::lexical_cast<Real>(current_line[4])*4.184*2);
        }catch(std::exception& e){
          throw ost::io::IOException("Could not cast bond parameters into real type!");
        }
        p->SetTypes(types);
        p->SetParam(parameters);
        ff_->AddImproper(p);
        break;
      }
      case 5:{ //cmap
        std::vector<String> types;
        std::vector<Real> parameters;
        types.push_back(current_line[0]);
        types.push_back(current_line[1]);
        types.push_back(current_line[2]);
        types.push_back(current_line[3]);
        types.push_back(current_line[7]);
        uint dim = boost::lexical_cast<uint>(current_line[8]);
        parameters.push_back(Real(dim));

        ++line_counter;

        while(true){
          if(line_counter >= content.size()){
            throw ost::io::IOException("Not enough data to fill cmap!");
          }
          if(parameters.size() > dim*dim+1){
            std::stringstream ss;
            ss << "Expect cmap values to fit exactly into a grid with dim ";
            ss << dim << " x " << dim;
            throw ost::io::IOException(ss.str());
          }
          if(parameters.size() == dim*dim+1){
            InteractionPtr p(new Interaction(CMap));
            p->SetTypes(types);
            p->SetParam(parameters);
            ff_->AddCMap(p);
            --line_counter;
            break;
          }

          for(std::vector<String>::iterator i = content[line_counter].begin();
              i != content[line_counter].end(); ++i){
            parameters.push_back(boost::lexical_cast<Real>(*i)*4.184);
          }
          ++line_counter;
        }

        break;
      }
      case 6:{ //nonbonded
        if(current_line[0] == "cutnb") break;
        std::vector<String> types;
        std::vector<Real> parameters;
        types.push_back(current_line[0]);
        try{
          parameters.push_back(boost::lexical_cast<Real>(current_line[3])*0.1/(std::pow(2,1.0/6.0))*2);
          parameters.push_back(std::abs(boost::lexical_cast<Real>(current_line[2])*4.184));
        }catch(std::exception& e){ 
          throw ost::io::IOException("Could not cast LJ parameters into real type!");
        }
        InteractionPtr p(new Interaction(LJ));
        p->SetTypes(types);
        p->SetParam(parameters);
        ff_->AddLJ(p);

        if(current_line.size() >= 7){
          try{
          lj_14_types.push_back(current_line[0]);
          lj_14_sigmas.push_back(boost::lexical_cast<Real>(current_line[6])*0.1/(std::pow(2,1.0/6.0))*2);
          lj_14_epsilons.push_back(std::abs(boost::lexical_cast<Real>(current_line[5])*4.184));
          }catch(std::exception& e){
            throw ost::io::IOException("Could not cast LJ parameters into real type!");
          }
        }else{
          lj_types.push_back(types[0]);
          lj_sigmas.push_back(parameters[0]);
          lj_epsilons.push_back(parameters[1]);
        }

        break;
      }
      case 7:{ //hbond
        //we don't implement any things regarding hbonds in CHARMM
        break;
      }
    }
  }
  //let's add all possible pairwise 1,4 interactions, where both interaction partners
  //have explicitely defined 1,4 interactions
  for(uint i = 0; i < lj_14_types.size(); ++i){
    for(uint j = i; j < lj_14_types.size(); ++j){
      InteractionPtr p(new Interaction(LJPair));
      std::vector<String> types;
      std::vector<Real> parameters;
      types.push_back(lj_14_types[i]);
      types.push_back(lj_14_types[j]);
      parameters.push_back(0.5 * (lj_14_sigmas[i] + lj_14_sigmas[j]));
      parameters.push_back(sqrt(lj_14_epsilons[i] * lj_14_epsilons[j]));
      p->SetTypes(types);
      p->SetParam(parameters);
      ff_->AddLJPair(p);
    }
  }
  //let's add all pairwise 1,4 interactions, where only one partner has explicitely defined 1,4 interactions
  for(uint i = 0; i < lj_14_types.size(); ++i){
    for(uint j = 0; j < lj_types.size(); ++j){
      InteractionPtr p(new Interaction(LJPair));
      std::vector<String> types;
      std::vector<Real> parameters;
      types.push_back(lj_14_types[i]);
      types.push_back(lj_types[j]);
      parameters.push_back(0.5 * (lj_14_sigmas[i] + lj_sigmas[j]));
      parameters.push_back(sqrt(lj_14_epsilons[i] * lj_epsilons[j]));
      p->SetTypes(types);
      p->SetParam(parameters);
      ff_->AddLJPair(p);
    }
  }


}

void FFReader::ParseCHARMMRTF(std::vector<std::vector<String> >& content){

  bool in_block = false;

  BuildingBlockPtr actual_block;

  String actual_block_name = "";

  int keyword_id = -1;
  std::vector<String> current_line;
  bool line_read = false;

  for(uint line_counter = 0; line_counter < content.size(); ++line_counter){
    current_line = content[line_counter];
    std::cerr<<current_line[0]<<std::endl;
    if(CHARMMData::Instance()->GetKeywordIndex(current_line[0]) != -1){
      keyword_id = CHARMMData::Instance()->GetKeywordIndex(current_line[0]);
    }
    if(keyword_id == -1) continue;

    

    if(in_block){
      line_read = false;
      switch(keyword_id){
        case 15:{ //atom
          actual_block->AddAtom(current_line[1],current_line[2],boost::lexical_cast<Real>(current_line[3]));
          line_read=true;
          break;
        }
        case 16:{//group
          line_read=true;
          break;
        }
        case 17:{//bond
          int num_bonds = (current_line.size()-1) / 2;
          for(int i = 0; i < num_bonds; ++i){
            InteractionPtr p(new Interaction(HarmonicBond));
            std::vector<String> names;
            names.push_back(current_line[i*2+1]);
            names.push_back(current_line[i*2+2]);
            p->SetNames(names);
            actual_block->AddBond(p);
          }
          line_read=true;
          break;
        }
        case 18:{//angle
          int num_angles = (current_line.size()-1) / 3;
          for(int i = 0; i < num_angles; ++i){
            InteractionPtr p(new Interaction(UreyBradleyAngle));
            std::vector<String> names;
            names.push_back(current_line[i*3+1]);
            names.push_back(current_line[i*3+2]);
            names.push_back(current_line[i*3+3]);
            p->SetNames(names);
            actual_block->AddAngle(p);
          }
          line_read=true;
          break;
        }
        case 19:{//dihedral
          int num_dihedrals = (current_line.size()-1) / 4;
          for(int i = 0; i < num_dihedrals; ++i){
            InteractionPtr p(new Interaction(PeriodicDihedral));
            std::vector<String> names;
            names.push_back(current_line[i*4+1]);
            names.push_back(current_line[i*4+2]);
            names.push_back(current_line[i*4+3]);
            names.push_back(current_line[i*4+4]);
            p->SetNames(names);
            actual_block->AddDihedral(p);
          }
          line_read=true;
          break;
        }
        case 20:{//improper
          int num_impropers = (current_line.size()-1) / 4;
          for(int i = 0; i < num_impropers; ++i){
            InteractionPtr p(new Interaction(HarmonicImproper));
            std::vector<String> names;
            names.push_back(current_line[i*4+1]);
            names.push_back(current_line[i*4+2]);
            names.push_back(current_line[i*4+3]);
            names.push_back(current_line[i*4+4]);
            p->SetNames(names);
            actual_block->AddImproper(p);
          }
          line_read=true;
          break;
        }
        case 5:{//cmap
          int num_cmaps = (current_line.size()-1) / 8;
          for(int i = 0; i < num_cmaps; ++i){
            InteractionPtr p(new Interaction(CMap));
            std::vector<String> names;
            names.push_back(current_line[i*4+1]);
            names.push_back(current_line[i*4+2]);
            names.push_back(current_line[i*4+3]);
            names.push_back(current_line[i*4+4]);
            names.push_back(current_line[i*4+8]);
            p->SetNames(names);
            actual_block->AddCMap(p);
          }
          line_read=true;
          break;
        }
        case 21:{ //donor
          line_read=true;
          break; //I don't think. that this will ever be read...
        }
        case 22:{//acceptor
          line_read=true;
          break; //I don't think. that this will ever be read...
        }
        case 23:{// ic
          line_read=true;
          break;
        }
        case 24:{ //patching
          line_read=true;
          break;
        }
        case 25:{ //print
          line_read=true;
          break;
        }
        case 26:{ //double =>basically a bond
          int num_bonds = (current_line.size()-1) / 2;
          for(int i = 0; i < num_bonds; ++i){
            InteractionPtr p(new Interaction(HarmonicBond));
            std::vector<String> names;
            names.push_back(current_line[i*2+1]);
            names.push_back(current_line[i*2+2]);
            p->SetNames(names);
            actual_block->AddBond(p);
          }
          line_read=true;
          break;
        }
      }
    }

    if(line_read) continue;

    switch(keyword_id){
      case 9:{
        try {ff_->AddMass(current_line[2],boost::lexical_cast<Real>(current_line[3]));}
        catch(std::exception& e){
          throw ost::io::IOException("Could not cast mass into real type!"); 
        }
        break;
      }
      case 10:{
        break; //we don't read that for now....
      }
      case 11:{//we don't allow patches for now
        break;
      }
      case 12:{
        break; //we don't read that for now
      }
      case 13:{
        //check, whether there is an "open block"
        if(in_block){
          ff_->AddBuildingBlock(actual_block_name,actual_block);
          std::cerr<<"added building block: "<<actual_block_name<<std::endl;
        }
        actual_block_name = current_line[1];
        actual_block = BuildingBlockPtr(new BuildingBlock);
        in_block = true;
        std::cerr<<actual_block_name<<std::endl;
        break;
      }
      case 14:{
        //check, whether there is an "open block"
        if(in_block){
          ff_->AddBuildingBlock(actual_block_name,actual_block);
          in_block = false;
        }
        break;
      }


    }

  }

  //add last building block if necessary
  if(in_block) ff_->AddBuildingBlock(actual_block_name,actual_block);

}

}}}
