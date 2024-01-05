#include <ost/io/mol/chemdict_parser.hh>

namespace ost { namespace io {

using namespace ost::conop;
  
bool ChemdictParser::OnBeginData(const StringRef& data_name) 
{
  compound_.reset(new Compound(data_name.str()));
  compound_->SetDialect(dialect_);
  if (last_!=data_name[0]) {
    last_=data_name[0];
    std::cout << last_ << std::flush;
  }
  atom_map_.clear();
  return true;
}

bool ChemdictParser::OnBeginLoop(const StarLoopDesc& header)
{ 
  if (header.GetCategory()=="chem_comp_atom") {
    loop_type_=ATOM_SPEC;
    indices_[ATOM_NAME]=header.GetIndex("atom_id");
    indices_[ALT_ATOM_NAME]=header.GetIndex("alt_atom_id");
    indices_[ELE]=header.GetIndex("type_symbol");      
    indices_[IS_LEAVING]=header.GetIndex("pdbx_leaving_atom_flag");
    indices_[IS_AROMATIC]=header.GetIndex("pdbx_aromatic_flag");
    indices_[ORDINAL]=header.GetIndex("pdbx_ordinal");
    indices_[CHARGE]=header.GetIndex("charge");
    return true;
  } else if (header.GetCategory()=="chem_comp_bond") {
    loop_type_=BOND_SPEC;
    indices_[ATOM_ID1]=header.GetIndex("atom_id_1");
    indices_[ATOM_ID2]=header.GetIndex("atom_id_2");
    indices_[BOND_ORDER]=header.GetIndex("value_order");
    return true;
  } else if (header.GetCategory()=="pdbx_chem_comp_descriptor") {
    loop_type_=DESC_SPEC;
    indices_[DESC_TYPE]=header.GetIndex("type");
    indices_[DESC]=header.GetIndex("descriptor");
    indices_[PROGRAM]=header.GetIndex("program");
    return true;
  }
  loop_type_=DONT_KNOW;
  return false;
}

void ChemdictParser::OnDataRow(const StarLoopDesc& header, 
                               const std::vector<StringRef>& columns) 
{                        
  if (loop_type_==ATOM_SPEC) {
    AtomSpec atom;
    atom.is_leaving=columns[indices_[IS_LEAVING]][0]=='Y';
    atom.name=columns[indices_[ATOM_NAME]].str();
    atom.alt_name=columns[indices_[ALT_ATOM_NAME]].str();
    atom.ordinal=columns[indices_[ORDINAL]].to_int().second-1;
    atom.element=columns[indices_[ELE]].str();
    atom.is_aromatic=columns[indices_[IS_AROMATIC]][0]=='Y';
    atom.charge=columns[indices_[CHARGE]].to_int().second;
    compound_->AddAtom(atom);
    atom_map_[atom.name]=atom.ordinal;
  } else if (loop_type_==BOND_SPEC) {
    BondSpec bond;
    String atom_one=columns[indices_[ATOM_ID1]].str();
    String atom_two=columns[indices_[ATOM_ID2]].str();      
    bond.atom_one=atom_map_[atom_one];
    bond.atom_two=atom_map_[atom_two];
    char b=columns[indices_[BOND_ORDER]][0];
    bond.order=b=='D' ? 2 : (b=='S' ? 1 : 3);
    if (bond.atom_one>bond.atom_two) {
      std::swap(bond.atom_one, bond.atom_two);
    }
    compound_->AddBond(bond);
  } else if (loop_type_==DESC_SPEC) {
    if (columns[indices_[DESC_TYPE]] == StringRef("InChI", 5)) {
      // for type InChi check prefix 'InChI='
      if (columns[indices_[DESC]].substr(0, 6) != StringRef("InChI=", 6)) {
        std::cout << "InChI problem: compound " << compound_->GetID()
                  << " has an InChI descriptor not starting with the "
                     "'InChI=' prefix." << std::endl;
        return;
      }
      compound_->SetInchi(columns[indices_[DESC]].str());
    } else if (columns[indices_[DESC_TYPE]] == StringRef("InChIKey", 8)) {
      compound_->SetInchiKey(columns[indices_[DESC]].str());
    } else if (columns[indices_[DESC_TYPE]] == StringRef("SMILES_CANONICAL", 16) &&
              columns[indices_[PROGRAM]] == StringRef("OpenEye OEToolkits", 18)) {
      compound_->SetSMILES(columns[indices_[DESC]].str());
    }
  }
}

void ChemdictParser::OnDataItem(const StarDataItem& item)
{
  if (item.GetCategory()==StringRef("chem_comp", 9)) {
    if (item.GetName()==StringRef("type", 4)) {
      // convert type to uppercase
      String type=item.GetValue().str();
      for (String::iterator i=type.begin(), e=type.end(); i!=e; ++i) {
        *i=toupper(*i);
      }
      // The type of water is set to "?". let's change it to water...
      if (compound_->GetID()=="HOH") {
        compound_->SetChemClass(mol::ChemClass(mol::ChemClass::WATER));
        compound_->SetOneLetterCode('.');
      } else {
        std::map<String, mol::ChemClass>::iterator i=tm_.find(type);
        if (i!=tm_.end()) {
          compound_->SetChemClass(i->second);
        } else {
          std::cout << "unknown type '" << type << "' for compound "
                    << compound_->GetID() << std::endl;
        }
      }
    } else if (item.GetName()==StringRef("pdbx_type", 9)) {
      String type=item.GetValue().str();
      for (String::iterator i=type.begin(), e=type.end(); i!=e; ++i) {
        *i=toupper(*i);
      }
      std::map<String, mol::ChemType>::iterator i=xtm_.find(type);
      if (i!=xtm_.end()) {
        compound_->SetChemType(i->second);
      } else {
        std::cout << "unknown pdbx_type '" << type << "' for compound "
                  << compound_->GetID() << std::endl;
      }
    } else if (item.GetName()==StringRef("name", 4)) {
      compound_->SetName(item.GetValue().str());
      if (compound_->GetName()==""){
        std::cout << "unknown compound name, chemcomp.name field empty for compound: "
                  << compound_->GetID() << std::endl;
      }
    } else if (item.GetName()==StringRef("formula", 7)) {
      compound_->SetFormula(item.GetValue().str());
      if (compound_->GetFormula()=="H2 O") {
        compound_->SetChemClass(mol::ChemClass(mol::ChemClass::WATER));
        compound_->SetOneLetterCode('.');
      }
    } else if (item.GetName()==StringRef("pdbx_release_status", 19)) {
      String release_status = item.GetValue().str();
      if (release_status == "OBS") {
        compound_->SetObsolete(true);
      }
    } else if (item.GetName()==StringRef("pdbx_replaced_by", 16)) {
      String replaced_by = item.GetValue().str();
      if (replaced_by != "?") {
        compound_->SetReplacedBy(replaced_by);
      }
    }  else if (item.GetName()==StringRef("one_letter_code", 15)) {
      if (item.GetValue().length()==1) {
        compound_->SetOneLetterCode(item.GetValue()[0]);   
      }
    } else if (item.GetName()==StringRef("pdbx_initial_date", 17)) {        
      compound_->SetCreationDate(Date::FromString(item.GetValue()));
    } else if (item.GetName()==StringRef("pdbx_modified_date", 18)) {
      compound_->SetModificationDate(Date::FromString(item.GetValue()));
    }
  } else if (item.GetName()==StringRef("atom_id", 7)) {
    atom_.name=item.GetValue().str();
  } else if (item.GetName()==StringRef("alt_atom_id", 11)) {
    if (compound_->GetID()=="ILE" && item.GetValue()==StringRef("CD1", 3)) {
      atom_.alt_name="CD";
    } else {
      atom_.alt_name=item.GetValue().str();       
    }

  } else if (item.GetName()==StringRef("type_symbol", 11)) {
    atom_.element=item.GetValue().str();
  } else if (item.GetName()==StringRef("pdbx_ordinal", 12)) {
    atom_.ordinal=item.GetValue().to_int().second-1;
  }
}

void ChemdictParser::OnEndData()
{
  if (compound_)
  {
    if (compound_->GetID() != "UNL" &&
        ! (ignore_reserved_ && IsNameReserved(compound_->GetID())) &&
        ! (ignore_obsolete_ && compound_->GetObsolete()))
    {
      if (compound_->GetAtomSpecs().empty()) {
        compound_->AddAtom(atom_);
      }
      lib_->AddCompound(compound_);
    }
  }
}

std::map<String, mol::ChemClass> ChemdictParser::tm_=std::map<String, mol::ChemClass>();
std::map<String, mol::ChemType> ChemdictParser::xtm_=std::map<String, mol::ChemType>();

void ChemdictParser::InitTypeMap()
{
  // This list is manually kept up to date to include all types appearing in PDB
  // See doc for how to create a compound library from scratch to test this
  // https://openstructure.org/docs/conop/compoundlib/#creating-a-compound-library
  if (!tm_.empty())
    return;
  tm_["L-PEPTIDE COOH CARBOXY TERMINUS"]=mol::ChemClass(mol::ChemClass::L_PEPTIDE_LINKING);
  tm_["L-PEPTIDE NH3 AMINO TERMINUS"]=mol::ChemClass(mol::ChemClass::L_PEPTIDE_LINKING);
  tm_["L-GAMMA-PEPTIDE, C-DELTA LINKING"]=mol::ChemClass(mol::ChemClass::L_PEPTIDE_LINKING);
  tm_["L-BETA-PEPTIDE, C-GAMMA LINKING"]=mol::ChemClass(mol::ChemClass::L_PEPTIDE_LINKING);
  tm_["D-PEPTIDE COOH CARBOXY TERMINUS"]=mol::ChemClass(mol::ChemClass::D_PEPTIDE_LINKING);
  tm_["D-PEPTIDE NH3 AMINO TERMINUS"]=mol::ChemClass(mol::ChemClass::D_PEPTIDE_LINKING);
  tm_["D-BETA-PEPTIDE, C-GAMMA LINKING"]=mol::ChemClass(mol::ChemClass::D_PEPTIDE_LINKING);
  tm_["D-GAMMA-PEPTIDE, C-DELTA LINKING"]=mol::ChemClass(mol::ChemClass::D_PEPTIDE_LINKING);
  tm_["L-SACCHARIDE, ALPHA LINKING"]=mol::ChemClass(mol::ChemClass::L_SACCHARIDE);
  tm_["L-SACCHARIDE, BETA LINKING"]=mol::ChemClass(mol::ChemClass::L_SACCHARIDE);
  tm_["L-SACCHARIDE 1,4 AND 1,4 LINKING"]=mol::ChemClass(mol::ChemClass::L_SACCHARIDE);
  tm_["L-SACCHARIDE 1,4 AND 1,6 LINKING"]=mol::ChemClass(mol::ChemClass::L_SACCHARIDE);
  tm_["D-SACCHARIDE 1,4 AND 1,4 LINKING"]=mol::ChemClass(mol::ChemClass::D_SACCHARIDE);
  tm_["L-SACCHARIDE"]=mol::ChemClass(mol::ChemClass::L_SACCHARIDE);
  tm_["D-SACCHARIDE"]=mol::ChemClass(mol::ChemClass::D_SACCHARIDE);
  tm_["D-SACCHARIDE, BETA LINKING"]=mol::ChemClass(mol::ChemClass::D_SACCHARIDE);
  tm_["D-SACCHARIDE, ALPHA LINKING"]=mol::ChemClass(mol::ChemClass::D_SACCHARIDE);
  tm_["SACCHARIDE"]=mol::ChemClass(mol::ChemClass::SACCHARIDE);
  tm_["D-PEPTIDE LINKING"]=mol::ChemClass(mol::ChemClass::D_PEPTIDE_LINKING);
  tm_["L-PEPTIDE LINKING"]=mol::ChemClass(mol::ChemClass::L_PEPTIDE_LINKING);
  tm_["L-PEPTIDE-LINKING"]=mol::ChemClass(mol::ChemClass::L_PEPTIDE_LINKING);
  tm_["DNA LINKING"]=mol::ChemClass(mol::ChemClass::DNA_LINKING);
  tm_["RNA LINKING"]=mol::ChemClass(mol::ChemClass::RNA_LINKING);
  tm_["L-DNA LINKING"]=mol::ChemClass(mol::ChemClass::DNA_LINKING);
  tm_["L-RNA LINKING"]=mol::ChemClass(mol::ChemClass::RNA_LINKING);
  tm_["R-DNA LINKING"]=mol::ChemClass(mol::ChemClass::DNA_LINKING);
  tm_["R-RNA LINKING"]=mol::ChemClass(mol::ChemClass::RNA_LINKING);
  tm_["DNA OH 3 PRIME TERMINUS"]=mol::ChemClass(mol::ChemClass::DNA_LINKING);
  tm_["DNA OH 5 PRIME TERMINUS"]=mol::ChemClass(mol::ChemClass::DNA_LINKING);
  tm_["PEPTIDE-LIKE"]=mol::ChemClass(mol::ChemClass::PEPTIDE_LINKING);
  tm_["PEPTIDE LINKING"]=mol::ChemClass(mol::ChemClass::PEPTIDE_LINKING);
  tm_["PEPTIDE-LINKING"]=mol::ChemClass(mol::ChemClass::PEPTIDE_LINKING);
  tm_["NON-POLYMER"]=mol::ChemClass(mol::ChemClass::NON_POLYMER);
  tm_["RNA OH 3 PRIME TERMINUS"]=mol::ChemClass(mol::ChemClass::RNA_LINKING);
  tm_["RNA OH 5 PRIME TERMINUS"]=mol::ChemClass(mol::ChemClass::RNA_LINKING);
  tm_["?"]=mol::ChemClass(mol::ChemClass::UNKNOWN);
  tm_["WATER"]=mol::ChemClass(mol::ChemClass::WATER);
  tm_["OTHER"]=mol::ChemClass(mol::ChemClass::UNKNOWN);
}

void ChemdictParser::InitPDBXTypeMap()
{
  if (!xtm_.empty())
    return;
  xtm_["HETAI"]=mol::ChemType(mol::ChemType::IONS);
  xtm_["HETAIN"]=mol::ChemType(mol::ChemType::NONCANONICALMOLS);
  xtm_["ATOMS"]=mol::ChemType(mol::ChemType::SACCHARIDES);
  xtm_["ATOMN"]=mol::ChemType(mol::ChemType::NUCLEOTIDES);
  xtm_["ATOMP"]=mol::ChemType(mol::ChemType::AMINOACIDS);
  xtm_["HETAC"]=mol::ChemType(mol::ChemType::COENZYMES);
  xtm_["HETIC"]=mol::ChemType(mol::ChemType::WATERCOORDIONS);
  xtm_["HETAD"]=mol::ChemType(mol::ChemType::DRUGS);
  xtm_["HETAS"]=mol::ChemType(mol::ChemType::WATERS);
  xtm_["?"]=mol::ChemType(mol::ChemType::UNKNOWN);
}

bool ChemdictParser::IsNameReserved(const String& data_name)
{
  // This checks if the compound name is one of the reserved name that will
  // never be used in the PDB. See
  // https://www.rcsb.org/news/feature/630fee4cebdf34532a949c34
  if ( // Check if in range 01-99
       (data_name.length() == 2 &&
         data_name[0] >= '0' && data_name[0] <= '9' &&
         data_name[1] >= '0' && data_name[1] <= '9' &&
         data_name != "00") ||
       data_name == "DRG" ||
       data_name == "INH" ||
       data_name == "LIG"
     )
  {
    std::cout << "Ignoring reserved compound " << data_name << std::endl;
    return true;
  }
  return false;
}

}}
