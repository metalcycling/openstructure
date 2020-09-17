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

#include <ctype.h>

#include <ost/profile.hh>
#include <ost/log.hh>
#include <ost/dyn_cast.hh>
#include <ost/mol/xcs_editor.hh>
#include <ost/conop/conop.hh>

#include <ost/io/mol/mmcif_reader.hh>

namespace ost { namespace io {


bool is_undef(StringRef value)
{
	return value.size()==1 && (value[0]=='?' || value[0]=='.');
}

MMCifReader::MMCifReader(std::istream& stream, mol::EntityHandle& ent_handle,
                         const IOProfile& profile):
  StarParser(stream, true), profile_(profile), ent_handle_(ent_handle)
{
  this->Init();
}

MMCifReader::MMCifReader(const String& filename, mol::EntityHandle& ent_handle,
                         const IOProfile& profile):
  StarParser(filename, true), profile_(profile), ent_handle_(ent_handle)
{
  this->Init();
}

void MMCifReader::Init()
{
  warned_name_mismatch_ = false;
  category_             = DONT_KNOW;
  memset(category_counts_, 0, DONT_KNOW * sizeof(int));
  chain_count_          = 0;
  atom_count_           = 0;
  residue_count_        = 0;
  auth_chain_id_        = false;
  seqres_can_           = false;
  has_model_            = false;
  restrict_chains_      = "";
  subst_res_id_         = "";
  curr_chain_           = mol::ChainHandle();
  curr_residue_         = mol::ResidueHandle();
  seqres_               = seq::CreateSequenceList();
  read_seqres_          = false;
  warned_rule_based_    = false;
  info_                 = MMCifInfo();
}

void MMCifReader::ClearState()
{
  curr_chain_           = mol::ChainHandle();
  curr_residue_         = mol::ResidueHandle();
  chain_count_          = 0;
  residue_count_        = 0;
  atom_count_           = 0;
  category_             = DONT_KNOW;
  warned_name_mismatch_ = false;
  seqres_               = seq::CreateSequenceList();
  info_                 = MMCifInfo();
  entity_desc_map_.clear();
  authors_map_.clear();
  bu_origin_map_.clear();
  bu_assemblies_.clear();
  helix_list_.clear();
  strand_list_.clear();
  revisions_.clear();
  revision_types_.clear();
  database_PDB_rev_added_ = false;
  entity_branch_link_map_.clear();
}

void MMCifReader::SetRestrictChains(const String& restrict_chains)
{
  restrict_chains_ = restrict_chains;
}

bool MMCifReader::OnBeginData(const StringRef& data_name) 
{
  LOG_DEBUG("MCIFFReader: " << profile_);
  Profile profile_import("MMCifReader::OnBeginData");

  // IDs in mmCIF files can be any string, so no restrictions here

  this->ClearState();

  return true;
}

bool MMCifReader::OnBeginLoop(const StarLoopDesc& header)
{
  bool cat_available = false;
  category_ = DONT_KNOW;
  // set whole array to -1
  memset(indices_, -1, MAX_ITEMS_IN_ROW * sizeof(int));
  // walk through possible categories
  if (header.GetCategory() == "atom_site") {
    category_ = ATOM_SITE;
    // mandatory items
    this->TryStoreIdx(AUTH_ASYM_ID,    "auth_asym_id",    header);
    this->TryStoreIdx(AS_ID,           "id",              header);
    this->TryStoreIdx(LABEL_ALT_ID,    "label_alt_id",    header);
    this->TryStoreIdx(LABEL_ASYM_ID,   "label_asym_id",   header);
    this->TryStoreIdx(LABEL_ATOM_ID,   "label_atom_id",   header);
    this->TryStoreIdx(LABEL_COMP_ID,   "label_comp_id",   header);
    this->TryStoreIdx(LABEL_ENTITY_ID, "label_entity_id", header);
    this->TryStoreIdx(LABEL_SEQ_ID,    "label_seq_id",    header);
    this->TryStoreIdx(TYPE_SYMBOL,     "type_symbol",     header);
    // assumed mandatory
    this->TryStoreIdx(CARTN_X, "Cartn_x", header);
    this->TryStoreIdx(CARTN_Y, "Cartn_y", header);
    this->TryStoreIdx(CARTN_Z, "Cartn_z", header);
    // optional
    indices_[OCCUPANCY]          = header.GetIndex("occupancy");
    indices_[B_ISO_OR_EQUIV]     = header.GetIndex("B_iso_or_equiv");
    indices_[GROUP_PDB]          = header.GetIndex("group_PDB");
    indices_[AUTH_SEQ_ID]        = header.GetIndex("auth_seq_id");
    indices_[PDBX_PDB_INS_CODE]  = header.GetIndex("pdbx_PDB_ins_code");
    indices_[PDBX_PDB_MODEL_NUM] = header.GetIndex("pdbx_PDB_model_num");

    // post processing
    if (category_counts_[category_] > 0) {
      if ((has_model_ && (indices_[PDBX_PDB_MODEL_NUM] == -1))||
          (!has_model_ && (indices_[PDBX_PDB_MODEL_NUM] != -1))) {
        throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                              "Not all atom_site entries carry a model number.",
                                                 this->GetCurrentLinenum()));
      }
    }
    cat_available = true;
  } else if (header.GetCategory() == "entity") {
    category_ = ENTITY;
    // mandatory items
    this->TryStoreIdx(E_ID, "id", header);
    // optional
    indices_[E_TYPE]  = header.GetIndex("type");
    indices_[PDBX_DESCRIPTION] = header.GetIndex("pdbx_description");
    cat_available = true;
  } else if (header.GetCategory() == "entity_poly") {
    category_ = ENTITY_POLY;
    // mandatory
    this->TryStoreIdx(ENTITY_ID, "entity_id", header);
    // optional
    indices_[EP_TYPE]  = header.GetIndex("type");
    indices_[PDBX_SEQ_ONE_LETTER_CODE] =
      header.GetIndex("pdbx_seq_one_letter_code");
    indices_[PDBX_SEQ_ONE_LETTER_CODE_CAN] =
      header.GetIndex("pdbx_seq_one_letter_code_can");
    cat_available = true;
  } else if (header.GetCategory() == "citation") {
    category_ = CITATION;
    // mandatory items
    this->TryStoreIdx(CITATION_ID, "id", header);
    // optional
    indices_[ABSTRACT_ID_CAS]         = header.GetIndex("abstract_id_CAS");
    indices_[BOOK_ID_ISBN]            = header.GetIndex("book_id_ISBN");
    indices_[BOOK_TITLE]              = header.GetIndex("book_title");
    indices_[BOOK_PUBLISHER]          = header.GetIndex("book_publisher");
    indices_[BOOK_PUBLISHER_CITY]     = header.GetIndex("book_publisher_city");
    indices_[JOURNAL_ABBREV]          = header.GetIndex("journal_abbrev");
    indices_[YEAR]                    = header.GetIndex("year");
    indices_[TITLE]                   = header.GetIndex("title");
    indices_[JOURNAL_VOLUME]          = header.GetIndex("journal_volume");
    indices_[PAGE_FIRST]              = header.GetIndex("page_first");
    indices_[PAGE_LAST]               = header.GetIndex("page_last");
    indices_[PDBX_DATABASE_ID_DOI]    = header.GetIndex("pdbx_database_id_DOI");
    indices_[PDBX_DATABASE_ID_PUBMED] =
      header.GetIndex("pdbx_database_id_PubMed");
    cat_available = true;
  } else if (header.GetCategory()=="citation_author") {
    category_ = CITATION_AUTHOR;
    // mandatory items
    this->TryStoreIdx(AUTHOR_CITATION_ID, "citation_id", header);
    this->TryStoreIdx(AUTHOR_NAME,        "name", header);
    this->TryStoreIdx(ORDINAL,            "ordinal", header);
    cat_available = true;
  } else if (header.GetCategory() == "exptl") {
    category_ = EXPTL;
    // mandatory items
    this->TryStoreIdx(CITATION_ID, "entry_id", header);
    this->TryStoreIdx(METHOD,      "method", header);
    cat_available = true;
  } else if (header.GetCategory() == "refine") {
    category_ = REFINE;
    // mandatory items
    this->TryStoreIdx(LS_D_RES_HIGH, "ls_d_res_high", header);
    // optional items
    indices_[REFINE_ENTRY_ID] = header.GetIndex("entry_id");
    indices_[LS_D_RES_LOW] = header.GetIndex("ls_d_res_low");
    indices_[LS_R_FACTOR_R_WORK] = header.GetIndex("ls_R_factor_R_work");
    indices_[LS_R_FACTOR_R_FREE] = header.GetIndex("ls_R_factor_R_free");
    cat_available = true;
  } else if (header.GetCategory() == "pdbx_struct_assembly") {
    category_ = PDBX_STRUCT_ASSEMBLY;
    // mandatory items
    this->TryStoreIdx(PSA_ID, "id", header);
    // optional
    indices_[PSA_DETAILS] = header.GetIndex("details");
    indices_[METHOD_DETAILS] = header.GetIndex("method_details");
    cat_available = true;
  } else if (header.GetCategory() == "pdbx_struct_assembly_gen") {
    category_ = PDBX_STRUCT_ASSEMBLY_GEN;
    // mandatory items
    this->TryStoreIdx(ASSEMBLY_ID,     "assembly_id", header);
    this->TryStoreIdx(ASYM_ID_LIST,    "asym_id_list", header);
    this->TryStoreIdx(OPER_EXPRESSION, "oper_expression", header);
    cat_available = true;
  } else if (header.GetCategory() == "pdbx_struct_oper_list") {
    category_ = PDBX_STRUCT_OPER_LIST;
    // mandatory items
    this->TryStoreIdx(PSOL_ID,   "id",   header);
    this->TryStoreIdx(PSOL_TYPE, "type", header);
    // optional items
    indices_[VECTOR_1]   = header.GetIndex("vector[1]");
    indices_[VECTOR_2]   = header.GetIndex("vector[2]");
    indices_[VECTOR_3]   = header.GetIndex("vector[3]");
    indices_[MATRIX_1_1] = header.GetIndex("matrix[1][1]");
    indices_[MATRIX_1_2] = header.GetIndex("matrix[1][2]");
    indices_[MATRIX_1_3] = header.GetIndex("matrix[1][3]");
    indices_[MATRIX_2_1] = header.GetIndex("matrix[2][1]");
    indices_[MATRIX_2_2] = header.GetIndex("matrix[2][2]");
    indices_[MATRIX_2_3] = header.GetIndex("matrix[2][3]");
    indices_[MATRIX_3_1] = header.GetIndex("matrix[3][1]");
    indices_[MATRIX_3_2] = header.GetIndex("matrix[3][2]");
    indices_[MATRIX_3_3] = header.GetIndex("matrix[3][3]");
    cat_available = true;
  } else if (header.GetCategory() == "struct") {
    category_ = STRUCT;
    // mandatory items
    this->TryStoreIdx(STRUCT_ENTRY_ID, "entry_id", header);
    // optional items
    indices_[PDBX_CASP_FLAG]             = header.GetIndex("pdbx_CASP_flag");
    indices_[PDBX_DESCRIPTOR]            = header.GetIndex("pdbx_descriptor");
    indices_[PDBX_FORMULA_WEIGHT]      = header.GetIndex("pdbx_formula_weight");
    indices_[PDBX_FORMULA_WEIGHT_METHOD]
       = header.GetIndex("pdbx_formula_weight_method");
    indices_[PDBX_MODEL_DETAILS]        = header.GetIndex("pdbx_model_details");
    indices_[PDBX_MODEL_TYPE_DETAILS]
       = header.GetIndex("pdbx_model_type_details");
    indices_[STRUCT_TITLE]               = header.GetIndex("title");
    cat_available = true;
  } else if (header.GetCategory() == "struct_conf") {
    category_ = STRUCT_CONF;
    // mandatory items
    this->TryStoreIdx(SC_BEG_LABEL_ASYM_ID, "beg_label_asym_id", header);
    this->TryStoreIdx(SC_BEG_LABEL_COMP_ID, "beg_label_comp_id", header);
    this->TryStoreIdx(SC_BEG_LABEL_SEQ_ID,  "beg_label_seq_id", header);
    this->TryStoreIdx(SC_CONF_TYPE_ID,      "conf_type_id", header);
    this->TryStoreIdx(SC_END_LABEL_ASYM_ID, "end_label_asym_id", header);
    this->TryStoreIdx(SC_END_LABEL_COMP_ID, "end_label_comp_id", header);
    this->TryStoreIdx(SC_END_LABEL_SEQ_ID,  "end_label_seq_id", header);
    this->TryStoreIdx(SC_ID,                "id", header);
    // optional items
    indices_[SC_BEG_AUTH_ASYM_ID] = header.GetIndex("beg_auth_asym_id");
    indices_[SC_END_AUTH_ASYM_ID] = header.GetIndex("end_auth_asym_id");
    cat_available = true;
  } else if (header.GetCategory() == "struct_sheet_range") {
    category_ = STRUCT_SHEET_RANGE;
    // mandatory items
    this->TryStoreIdx(SSR_BEG_LABEL_ASYM_ID,     "beg_label_asym_id", header);
    this->TryStoreIdx(SSR_BEG_LABEL_COMP_ID,     "beg_label_comp_id", header);
    this->TryStoreIdx(SSR_BEG_LABEL_SEQ_ID,      "beg_label_seq_id", header);
    this->TryStoreIdx(SSR_END_LABEL_ASYM_ID,     "end_label_asym_id", header);
    this->TryStoreIdx(SSR_END_LABEL_COMP_ID,     "end_label_comp_id", header);
    this->TryStoreIdx(SSR_END_LABEL_SEQ_ID,      "end_label_seq_id", header);
    this->TryStoreIdx(SSR_SHEET_ID,              "sheet_id", header);
    this->TryStoreIdx(SSR_ID,                    "id", header);
    // optional items
    indices_[SSR_BEG_AUTH_ASYM_ID] = header.GetIndex("beg_auth_asym_id");
    indices_[SSR_END_AUTH_ASYM_ID] = header.GetIndex("end_auth_asym_id");
    cat_available = true;
  } else if (header.GetCategory() == "pdbx_database_PDB_obs_spr") {
    category_ = PDBX_DATABASE_PDB_OBS_SPR;
    // mandatory items
    this->TryStoreIdx(DATE,           "date", header);
    this->TryStoreIdx(PDPOS_ID,       "id", header);
    this->TryStoreIdx(PDB_ID,         "pdb_id", header);
    this->TryStoreIdx(REPLACE_PDB_ID, "replace_pdb_id", header);
    cat_available = true;
  } else if (header.GetCategory() == "struct_ref") {
    category_ = STRUCT_REF;
    this->TryStoreIdx(SR_ENTITY_ID, "entity_id", header);
    this->TryStoreIdx(SR_ID, "id", header);
    this->TryStoreIdx(SR_DB_NAME, "db_name", header);
    this->TryStoreIdx(SR_DB_CODE, "db_code", header);
    indices_[SR_DB_ACCESS]=header.GetIndex("pdbx_db_accession");
    cat_available = true;
  } else if (header.GetCategory() == "struct_ref_seq") {
    category_ = STRUCT_REF_SEQ;	
    this->TryStoreIdx(SRS_ALIGN_ID, "align_id", header);
    this->TryStoreIdx(SRS_STRUCT_REF_ID, "ref_id", header);
    this->TryStoreIdx(SRS_ENT_ALIGN_BEG, "seq_align_beg", header);
    this->TryStoreIdx(SRS_ENT_ALIGN_END, "seq_align_end", header);
    this->TryStoreIdx(SRS_DB_ALIGN_BEG, "db_align_beg", header);
    this->TryStoreIdx(SRS_DB_ALIGN_END, "db_align_end", header);
    indices_[SRS_PDBX_STRAND_ID]=header.GetIndex("pdbx_strand_id");
    cat_available = true;
  } else if (header.GetCategory()=="struct_ref_seq_dif") {
    category_ = STRUCT_REF_SEQ_DIF;
    this->TryStoreIdx(SRSD_ALIGN_ID, "align_id", header);
    this->TryStoreIdx(SRSD_SEQ_RNUM, "seq_num", header);
    this->TryStoreIdx(SRSD_DB_RNUM, "pdbx_seq_db_seq_num", header);
    indices_[SRSD_DETAILS]=header.GetIndex("details");
    cat_available = true;
  } else if (header.GetCategory()=="database_PDB_rev") {
    // THIS IS FOR mmCIF versions < 5
    category_ = DATABASE_PDB_REV;
    // mandatory items
    this->TryStoreIdx(DPI_NUM, "num", header);
    // optional items
    indices_[DPI_DATE] = header.GetIndex("date");
    indices_[DPI_DATE_ORIGINAL] = header.GetIndex("date_original");
    indices_[DPI_STATUS] = header.GetIndex("status");
    cat_available = true;
  } else if (header.GetCategory()=="pdbx_audit_revision_history") {
    // THIS IS FOR mmCIF versions >= 5
    category_ = PDBX_AUDIT_REVISION_HISTORY;
    // mandatory items
    this->TryStoreIdx(PARH_ORDINAL, "ordinal", header);
    this->TryStoreIdx(PARH_MAJOR, "major_revision", header);
    this->TryStoreIdx(PARH_MINOR, "minor_revision", header);
    this->TryStoreIdx(PARH_REVISION_DATE, "revision_date", header);
    cat_available = true;
  } else if (header.GetCategory()=="pdbx_audit_revision_details") {
    // THIS IS FOR mmCIF versions >= 5
    category_ = PDBX_AUDIT_REVISION_DETAILS;
    // mandatory items
    this->TryStoreIdx(PARD_REVISION_ORDINAL, "revision_ordinal", header);
    // optional items
    indices_[PARD_TYPE] = header.GetIndex("type");
    if (indices_[PARD_TYPE] == -1) {
      LOG_WARNING("No 'pdbx_audit_revision_details.type' items "
                  "found! The revision history will not have status entries!");
    }
    cat_available = true;
  } else if (header.GetCategory()=="pdbx_database_status") {
    // THIS IS FOR mmCIF versions >= 5
    category_ = PDBX_DATABASE_STATUS;
    // optional items
    indices_[PDS_RECVD_INITIAL_DEPOSITION_DATE]
     = header.GetIndex("recvd_initial_deposition_date");
    cat_available = true;
  } else if (header.GetCategory() == "pdbx_entity_branch") {
    category_ = PDBX_ENTITY_BRANCH;
    // mandatory
    this->TryStoreIdx(BR_ENTITY_ID, "entity_id", header);
    this->TryStoreIdx(BR_ENTITY_TYPE, "type", header); 
    cat_available = true;
 } else if (header.GetCategory() == "pdbx_entity_branch_link") {
    category_ = PDBX_ENTITY_BRANCH_LINK;
    // mandatory
    this->TryStoreIdx(BL_ENTITY_ID, "entity_id", header);
    this->TryStoreIdx(BL_ATOM_ID_1, "atom_id_1", header);
    this->TryStoreIdx(BL_ATOM_ID_2, "atom_id_2", header);
    this->TryStoreIdx(BL_COMP_ID_1, "comp_id_1", header);
    this->TryStoreIdx(BL_COMP_ID_2, "comp_id_2", header);
    this->TryStoreIdx(BL_ENTITY_BRANCH_LIST_NUM_1, "entity_branch_list_num_1",
                      header);
    this->TryStoreIdx(BL_ENTITY_BRANCH_LIST_NUM_2, "entity_branch_list_num_2",
                      header);
    // optional items
    indices_[BL_ATOM_STEREO_CONFIG_1] = header.GetIndex("atom_stereo_config_1");
    indices_[BL_ATOM_STEREO_CONFIG_2] = header.GetIndex("atom_stereo_config_2");
    indices_[BL_VALUE_ORDER] = header.GetIndex("value_order");
    cat_available = true;
 }
  category_counts_[category_]++;
  return cat_available;
}

mol::ResNum to_res_num(int num, char ins_code)
{
  return mol::ResNum(num, ins_code==' ' ? '\0' : ins_code);
}

bool MMCifReader::ParseAtomIdent(const std::vector<StringRef>& columns,
                                 String& auth_chain_name,
                                 String& cif_chain_name,
                                 StringRef& res_name,
                                 mol::ResNum& resnum,
                                 bool& valid_res_num,
                                 StringRef& atom_name,
                                 char& alt_loc)
{
  // ATOM name
  atom_name = columns[indices_[LABEL_ATOM_ID]];
  if (profile_.calpha_only) {
    if (atom_name != StringRef("CA", 2)) {
      return false;
    }
  }
  // CHAIN ID
  auth_chain_name = columns[indices_[AUTH_ASYM_ID]].str();
  cif_chain_name = columns[indices_[LABEL_ASYM_ID]].str();

  if (! restrict_chains_.empty() &&
      restrict_chains_.find(cif_chain_name) == String::npos) {
    return false;
  } 

  this->TryGetInt(columns[indices_[AS_ID]],
                  "atom_site.id",
                  profile_.fault_tolerant); // unit test

  alt_loc = columns[indices_[LABEL_ALT_ID]][0];
  res_name = columns[indices_[LABEL_COMP_ID]];
  std::pair<bool, int> res_num;
  if (columns[indices_[LABEL_SEQ_ID]][0] != '.') {
    res_num =this->TryGetInt(columns[indices_[LABEL_SEQ_ID]],
                             "atom_site.label_seq_id",
                             profile_.fault_tolerant); // unit test
    if (!res_num.first) { // unit test
      if (profile_.fault_tolerant) {
        return false;
      }
    }
    valid_res_num = true;
  } else {
    valid_res_num = false;
    return true;
  }

  resnum = to_res_num(res_num.second, ' ');

  return true;
}

void MMCifReader::ParseAndAddAtom(const std::vector<StringRef>& columns)
{
  mol::XCSEditor editor=ent_handle_.EditXCS(mol::BUFFERED_EDIT); // potbl
  char alt_loc=0;
  String auth_chain_name;
  String cif_chain_name;
  StringRef res_name, atom_name;
  mol::ResNum res_num(0);
  bool valid_res_num = false;
  if (indices_[PDBX_PDB_MODEL_NUM] != -1) {
    if (has_model_) {
      if (curr_model_ != TryGetInt(columns[indices_[PDBX_PDB_MODEL_NUM]],
                                   "atom_site.pdbx_PDB_model_num")) {
        return;
      }
    } else {
      has_model_ = true;
      curr_model_ = TryGetInt(columns[indices_[PDBX_PDB_MODEL_NUM]],
      "atom_site.pdbx_PDB_model_num");
    }
  }

  if (!this->ParseAtomIdent(columns,
                            auth_chain_name,
                            cif_chain_name,
                            res_name,
                            res_num,
                            valid_res_num,
                            atom_name,
                            alt_loc)) {// unit test
    return;                            
  }
  Real occ = 1.00f, temp = 0;
  geom::Vec3 apos;
  
  for (int i = CARTN_X; i <= CARTN_Z; ++i) {
    std::pair<bool, float> result = this->TryGetFloat(columns[indices_[i]],
                                                      "atom_site.cartn_[xyz]",
                                                      profile_.fault_tolerant);
    if (!result.first) { // unit test
      if (profile_.fault_tolerant) { // unit test
        return;
      }
    }
    apos[i - CARTN_X] = result.second;
  }

  if (indices_[OCCUPANCY] != -1) { // unit test
    occ = this->TryGetReal(columns[indices_[OCCUPANCY]], "atom_site.occupancy");
  }
  if (indices_[B_ISO_OR_EQUIV] != -1) { // unit test
    temp = this->TryGetReal(columns[indices_[B_ISO_OR_EQUIV]],
                            "atom_site.B_iso_or_equiv");
  }

  // determine element
  String s_ele(columns[indices_[TYPE_SYMBOL]].str());

  String aname(atom_name.str());  
  // some postprocessing
  LOG_TRACE( "s_chain: [" << cif_chain_name << "]" );

  // determine chain and residue update
  bool update_chain = false;
  bool update_residue = false;
  if(!curr_chain_) { // unit test
      update_chain=true;
      update_residue=true;
  } else if(curr_chain_.GetName() != cif_chain_name) { // unit test
    update_chain=true;
    update_residue=true;
  }

  if(!curr_residue_) { // unit test
    update_residue=true;
    subst_res_id_ = cif_chain_name +
                    columns[indices_[AUTH_SEQ_ID]].str() +
                    columns[indices_[PDBX_PDB_INS_CODE]].str();
  } else if (!valid_res_num) {
    if (indices_[AUTH_SEQ_ID] != -1 &&
        indices_[PDBX_PDB_INS_CODE] != -1) {
      if (subst_res_id_ !=
          cif_chain_name +
          columns[indices_[AUTH_SEQ_ID]].str() +
          columns[indices_[PDBX_PDB_INS_CODE]].str()) {
        update_residue=true;

        subst_res_id_ = cif_chain_name +
                        columns[indices_[AUTH_SEQ_ID]].str() +
                        columns[indices_[PDBX_PDB_INS_CODE]].str();
      }
    } else {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                           "Missing residue number information",
                                               this->GetCurrentLinenum()));
    }
  } else if(curr_residue_.GetNumber() != res_num) { // unit test
    update_residue=true;
  }

  if(update_chain) { // unit test
    curr_chain_ = ent_handle_.FindChain(cif_chain_name);
    if(!curr_chain_.IsValid()) { // unit test
      LOG_DEBUG("new chain " << cif_chain_name);
      curr_chain_=editor.InsertChain(cif_chain_name);
      curr_chain_.SetStringProp("pdb_auth_chain_name", auth_chain_name);
      ++chain_count_;
      // store entity id
      String ent_id = columns[indices_[LABEL_ENTITY_ID]].str();
      chain_id_pairs_.push_back(std::pair<mol::ChainHandle,String>(curr_chain_,
                                                                   ent_id));
      info_.AddMMCifEntityIdTr(cif_chain_name, ent_id);
    }
    assert(curr_chain_.IsValid());
  } else if (chain_id_pairs_.back().second != // unit test
             columns[indices_[LABEL_ENTITY_ID]].str()) {
    // check that label_entity_id stays the same
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
        "Change of 'atom_site.label_entity_id' item for chain " +
        curr_chain_.GetName() + "! Expected: " + chain_id_pairs_.back().second +
        ", found: " + columns[indices_[LABEL_ENTITY_ID]].str() + ".",
                                             this->GetCurrentLinenum()));
  }

  if(update_residue) { // unit test
    curr_residue_=mol::ResidueHandle();
    if (valid_res_num && profile_.join_spread_atom_records) { // unit test
      curr_residue_=curr_chain_.FindResidue(res_num);
    }
    if (!curr_residue_.IsValid()) { // unit test
      LOG_DEBUG("new residue " << res_name << " " << res_num);
      if (valid_res_num) {
        curr_residue_ = editor.AppendResidue(curr_chain_,
                                             res_name.str(),
                                             res_num);
      } else {
        curr_residue_ = editor.AppendResidue(curr_chain_, res_name.str());
      }
      warned_name_mismatch_=false;
      ++residue_count_; 
    }
    assert(curr_residue_.IsValid());
  }

  // finally add atom
  LOG_DEBUG("adding atom " << aname << " (" << s_ele << ") @" << apos);
  mol::AtomHandle ah;
  if (curr_residue_.GetName()!=res_name.str()) { // unit test
    if (!profile_.fault_tolerant && alt_loc=='.') { // unit test
      std::stringstream ss;
      ss << "Residue with number " << res_num << " has more than one name.";
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                               ss.str(),
                                               this->GetCurrentLinenum()));
    }
    if (!warned_name_mismatch_) { // unit test
      if (alt_loc=='.') { // unit test
        LOG_WARNING("Residue with number " << res_num << " has more than one "
                    "name. Ignoring atoms for everything but the first");
      } else {
        LOG_WARNING("Residue with number " << res_num 
                    << " contains a microheterogeneity. Everything but atoms "
                    "for the residue '" << curr_residue_.GetName() 
                    << "' will be ignored");
      }
    }
    warned_name_mismatch_=true;
    return;
  }
  if (alt_loc!='.') { // unit test
    // Check if there is already a atom with the same name.
    mol::AtomHandle me=curr_residue_.FindAtom(aname);
    if (me.IsValid()) { // unit test
      try {
        editor.AddAltAtomPos(String(1, alt_loc), me, apos);
      } catch (Error&) {
        LOG_INFO("Ignoring atom alt location since there is already an atom "
                 "with name " << aname << ", but without an alt loc");
        return;
      }
      return;
    } else {
      ah = editor.InsertAltAtom(curr_residue_, aname,
                                String(1, alt_loc), apos, s_ele);
      ++atom_count_;
      }
  } else {
    mol::AtomHandle atom=curr_residue_.FindAtom(aname);
    if (atom.IsValid() && !profile_.quack_mode) { // unit test
      if (profile_.fault_tolerant) { // unit test
        LOG_WARNING("duplicate atom '" << aname << "' in residue " 
                    << curr_residue_);
        return;
      }
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                               "Duplicate atom '"+aname+
                                               "' in residue "+
                                               curr_residue_.GetQualifiedName(),
                                               this->GetCurrentLinenum()));
    }
    ah = editor.InsertAtom(curr_residue_, aname, apos, s_ele);
    ++atom_count_;
  }
  ah.SetBFactor(temp);

  ah.SetOccupancy(occ);

  // record type
  ah.SetHetAtom(indices_[GROUP_PDB] == -1 ? false :  
                columns[indices_[GROUP_PDB]][0]=='H');
}

MMCifReader::MMCifEntityDescMap::iterator MMCifReader::GetEntityDescMapIterator(
  const String& entity_id)
{
  MMCifEntityDescMap::iterator edm_it = entity_desc_map_.find(entity_id);
  // if the entity ID is not already stored, insert it with empty values
  if (edm_it == entity_desc_map_.end()) {
    MMCifEntityDesc desc = {.type=mol::CHAINTYPE_N_CHAINTYPES,
                            .details="",
                            .seqres=""};
    edm_it = entity_desc_map_.insert(entity_desc_map_.begin(),
                                     MMCifEntityDescMap::value_type(entity_id,
                                                                    desc));
  }
  return edm_it;
}

void MMCifReader::ParseEntity(const std::vector<StringRef>& columns)
{
  MMCifEntityDescMap::iterator edm_it =
    GetEntityDescMapIterator(columns[indices_[E_ID]].str());

  // type
  if (indices_[E_TYPE] != -1) {
    // only use the entity type if no other is set, entity_poly type is
    // more precise, so if that was set before just leave it in
    if (edm_it->second.type == mol::CHAINTYPE_N_CHAINTYPES) {
      edm_it->second.type = mol::ChainTypeFromString(columns[indices_[E_TYPE]]);
    }
  } else {
    // don't deal with entities without type
    entity_desc_map_.erase(edm_it);
    return;
  }

  // description
  if (indices_[PDBX_DESCRIPTION] != -1) {
    edm_it->second.details = columns[indices_[PDBX_DESCRIPTION]].str();
  }
}

void MMCifReader::ParseEntityPoly(const std::vector<StringRef>& columns)
{
  MMCifEntityDescMap::iterator edm_it =
    GetEntityDescMapIterator(columns[indices_[ENTITY_ID]].str());

  // store type
  if (indices_[EP_TYPE] != -1) {
    edm_it->second.type = mol::ChainTypeFromString(columns[indices_[EP_TYPE]]);
  }

  // store seqres
  if (edm_it->second.seqres.length() > 0) {
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
     "entity_poly.pdbx_seq_one_letter_code[_can] clash: sequence for entry '" +
                                            columns[indices_[ENTITY_ID]].str() +
                                             "' is already set to '" +
                                             edm_it->second.seqres + "'.",
                                             this->GetCurrentLinenum()));
  }
  if (read_seqres_) {
    StringRef seqres;
    if (seqres_can_) {
      if (indices_[PDBX_SEQ_ONE_LETTER_CODE_CAN] != -1) {
        seqres=columns[indices_[PDBX_SEQ_ONE_LETTER_CODE_CAN]];
        edm_it->second.seqres = seqres.str_no_whitespace();        
      } else {
        throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                   "'entity_poly.pdbx_seq_one_letter_code_can' not available.'",
                                                 this->GetCurrentLinenum()));
      }
    } else if (indices_[PDBX_SEQ_ONE_LETTER_CODE] != -1) {
      seqres=columns[indices_[PDBX_SEQ_ONE_LETTER_CODE]];

      conop::CompoundLibPtr comp_lib=conop::Conopology::Instance()
                                            .GetDefaultLib();
      if (!comp_lib) {
        if (!warned_rule_based_) {
          LOG_WARNING("SEQRES import requires a compound library. "
                       "Ignoring SEQRES records");      
        }
        warned_rule_based_=true;
        return;
      }
      edm_it->second.seqres = this->ConvertSEQRES(seqres.str_no_whitespace(),
                                                  comp_lib);
    } else {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                       "'entity_poly.pdbx_seq_one_letter_code' not available.'",
                                               this->GetCurrentLinenum()));
    }
  }
}

String MMCifReader::ConvertSEQRES(const String& seqres, 
                                  conop::CompoundLibPtr comp_lib)
{
  String can_seqres;
  for (String::const_iterator i=seqres.begin(), e=seqres.end(); i!=e; ++i) {
    if (*i=='(') {
      bool found_end_paren=false;
      String tlc;
      tlc.reserve(3);
      while ((++i)!=seqres.end()) {
        if (*i==')') {
          found_end_paren=true;
          break;
        }
        tlc.push_back(*i);
      }
      if (!found_end_paren) {
        throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                          "'entity_poly.pdbx_seq_one_letter_code' contains "
                          "unmatched '('", this->GetCurrentLinenum()));
      }
      conop::CompoundPtr compound=comp_lib->FindCompound(tlc, 
                                                         conop::Compound::PDB);
      if (!compound) {
        if (tlc!="UNK") {

          LOG_WARNING("unknown residue '" << tlc << "' in SEQRES record. "
                      "Setting one-letter-code to 'X'");
        }
        can_seqres.push_back('X');
        continue;
      }
      if (compound->GetOneLetterCode()=='?') {
        can_seqres.push_back('X');
      } else {
        can_seqres.push_back(compound->GetOneLetterCode());
      }

    } else {
      can_seqres.push_back(*i);
    }
  }
  return can_seqres;
}

void MMCifReader::ParseCitation(const std::vector<StringRef>& columns)
{
  // create citation object
  MMCifInfoCitation cit = MMCifInfoCitation();
  // just add info
  cit.SetID(columns[indices_[CITATION_ID]].str());
  if (indices_[ABSTRACT_ID_CAS] != -1) {
    if (columns[indices_[ABSTRACT_ID_CAS]][0]!='?') {
      cit.SetCAS(columns[indices_[ABSTRACT_ID_CAS]].str());
    }
  }
  if (indices_[BOOK_ID_ISBN] != -1) {
    if (columns[indices_[BOOK_ID_ISBN]][0]!='?') {
      cit.SetISBN(columns[indices_[BOOK_ID_ISBN]].str());
    }
  }
  if (indices_[JOURNAL_ABBREV] != -1) {
    if ((columns[indices_[JOURNAL_ABBREV]] != StringRef(".", 1)) &&
        (columns[indices_[JOURNAL_ABBREV]][0] != '?')) {
          cit.SetPublishedIn(columns[indices_[JOURNAL_ABBREV]].str());
          cit.SetCitationTypeJournal();
        }
  }
  if (indices_[BOOK_TITLE] != -1) {
    // this is only set in few PDB entries and RCSB overrides it with
    // the journal_abbrev for their citations
    // -> as of August 1, 2017, 5 entries known: 5b1j, 5b1k, 5fax, 5fbz, 5ffn
    //    -> all those have journal_abbrev set
    if ((columns[indices_[BOOK_TITLE]] != StringRef(".", 1)) &&
        (columns[indices_[BOOK_TITLE]][0] != '?')) {
      // This will override published_in if already set by journal_abbrev. We
      // consider this OK for now since usually the book title is copied to
      // the journal_abbrev attribute.
      cit.SetPublishedIn(columns[indices_[BOOK_TITLE]].str());
      cit.SetCitationTypeBook();
      
      // In theory, book_publisher and book_publisher_city are only set for
      // books and book chapters, so we only try to fetch them if the citation
      // type points to book.
      if (indices_[BOOK_PUBLISHER] != -1) {
        cit.SetBookPublisher(columns[indices_[BOOK_PUBLISHER]].str());
      }
      if (indices_[BOOK_PUBLISHER_CITY] != -1) {
        cit.SetBookPublisherCity(columns[indices_[BOOK_PUBLISHER_CITY]].str());
      }
    }
  }
  if (indices_[JOURNAL_VOLUME] != -1) {
    if (columns[indices_[JOURNAL_VOLUME]][0]!='?') {
      cit.SetVolume(columns[indices_[JOURNAL_VOLUME]].str());
    }
  }
  if (indices_[PAGE_FIRST] != -1) {
    if (columns[indices_[PAGE_FIRST]][0]!='?') {
      cit.SetPageFirst(columns[indices_[PAGE_FIRST]].str());
    }
  }
  if (indices_[PAGE_LAST] != -1) {
    if (columns[indices_[PAGE_LAST]][0]!='?') {
      cit.SetPageLast(columns[indices_[PAGE_LAST]].str());
    }
  }
  if (indices_[PDBX_DATABASE_ID_DOI] != -1) {
    if (columns[indices_[PDBX_DATABASE_ID_DOI]][0]!='?') {
      cit.SetDOI(columns[indices_[PDBX_DATABASE_ID_DOI]].str());
    }
  }
  if (indices_[PDBX_DATABASE_ID_PUBMED] != -1) {
    if (!is_undef(columns[indices_[PDBX_DATABASE_ID_PUBMED]])) {
      cit.SetPubMed(this->TryGetInt(columns[indices_[PDBX_DATABASE_ID_PUBMED]],
                                    "citation.pdbx_database_id_PubMed"));
    }
  }
  if (indices_[YEAR] != -1) {
    if (!is_undef(columns[indices_[YEAR]])) {
      cit.SetYear(this->TryGetInt(columns[indices_[YEAR]], "citation.year"));
    }
  }
  if (indices_[TITLE] != -1) {
    cit.SetTitle(columns[indices_[TITLE]].str());
  }

  // store citation (wo author, yet)
  info_.AddCitation(cit);
}

void MMCifReader::ParseCitationAuthor(const std::vector<StringRef>& columns)
{
  // get/ pack values
  MMCifCitationAuthorMap::iterator atm_it;
  std::vector<String> at_vec;
  std::vector<int> pos_vec;
  atm_it = authors_map_.find(columns[indices_[AUTHOR_CITATION_ID]].str());
  if (atm_it != authors_map_.end()) {
    at_vec = atm_it->second.second;
    pos_vec = atm_it->second.first;
  }
  at_vec.push_back(columns[indices_[AUTHOR_NAME]].str());
  pos_vec.push_back(this->TryGetInt(columns[indices_[ORDINAL]],
                            "citation_author.ordinal"));

  // sort new author into right position
  std::vector<int>::iterator pos_it;
  std::vector<String>::iterator atv_it;
  int ti;
  String ts; 
  pos_it = pos_vec.end();
  atv_it = at_vec.end();
  --pos_it;
  --atv_it;
  for (; pos_it != pos_vec.begin(); --pos_it, --atv_it) {
    if (*pos_it < *(pos_it-1)) {
      ti = *pos_it;
      *pos_it = *(pos_it-1);
      *(pos_it-1) = ti;
      ts = *atv_it;
      *atv_it = *(atv_it-1);
      *(atv_it-1) = ts;
    }
    else {
      break;
    }
  }

  // store new values in map
  if (atm_it != authors_map_.end()) {
    atm_it->second.second = at_vec;
    atm_it->second.first  = pos_vec;
  } else {
    authors_map_.insert(MMCifCitationAuthorMap::value_type(
                               columns[indices_[AUTHOR_CITATION_ID]].str(),
              std::pair<std::vector<int>, std::vector<String> >(pos_vec, at_vec)
                               ));
  }
}

void MMCifReader::ParseExptl(const std::vector<StringRef>& columns)
{
  info_.SetMethod(columns[indices_[METHOD]].str());
}

void MMCifReader::ParseRefine(const std::vector<StringRef>& columns)
{
  StringRef col = columns[indices_[LS_D_RES_HIGH]];
  if (col.size()!=1 || (col[0]!='?' && col[0]!='.')) {
    info_.SetResolution(this->TryGetReal(col, "refine.ls_d_res_high"));
  }
  if (indices_[LS_R_FACTOR_R_WORK] != -1) {
    col = columns[indices_[LS_R_FACTOR_R_WORK]];
    if (col.size()!=1 || (col[0]!='?' && col[0]!='.')) {
      info_.SetRWork(this->TryGetReal(col, "refine.ls_R_factor_R_work"));
    }
  }
  if (indices_[LS_R_FACTOR_R_FREE] != -1) {
    col = columns[indices_[LS_R_FACTOR_R_FREE]];
    if (col.size()!=1 || (col[0]!='?' && col[0]!='.')) {
      info_.SetRFree(this->TryGetReal(col, "refine.ls_R_factor_R_free"));
    }
  }
}

void MMCifReader::ParsePdbxStructAssembly(const std::vector<StringRef>& columns)
{
  MMCifPSAEntry psa;

  if (indices_[PSA_DETAILS] != -1) {
    psa.details = columns[indices_[PSA_DETAILS]].str();
  } else {
    psa.details = "?";
  }

  if (indices_[METHOD_DETAILS] != -1) {
    psa.method_details = columns[indices_[METHOD_DETAILS]].str();
  } else {
    psa.method_details = "?";
  }

  bu_origin_map_.insert(std::pair<String,
                         MMCifPSAEntry>(columns[indices_[PSA_ID]].str(), psa));
}

void MMCifReader::StoreExpression(const char* l, const char* s,
                                  bool& is_range, int lborder,
                                  std::vector<String>& single_block)
{
  std::stringstream ss;
  int rborder;

  if (l != s) {
    if (is_range) {
      is_range = false;
      rborder = this->TryGetInt(StringRef(l, s-l),
                                "pdbx_struct_assembly_gen.oper_expression");
      for (lborder += 1; lborder < rborder; lborder++) {
        ss << lborder;
        single_block.push_back(ss.str());
        ss.str("");
      }
    }
    single_block.push_back(String(l, s-l));
  }
}

void MMCifReader::StoreRange(const char*& l, const char* s, bool& is_range,
                             int& lborder, std::vector<String>& single_block)
{
  if (is_range) {
    throw IOException(this->FormatDiagnostic(STAR_DIAG_WARNING,
                                             "pdbx_struct_assembly_gen.oper_expression is missing a right border for a range expression.",
                                             this->GetCurrentLinenum()));
  }
  is_range = true;
  if (l != s) {
    lborder = this->TryGetInt(StringRef(l, s-l),
                              "pdbx_struct_assembly_gen.oper_expression");
    single_block.push_back(String(l, s-l));
  }
  l = s+1;
}

std::vector<std::vector<String> > MMCifReader::UnPackOperExperession(StringRef expression)
{
  std::vector<std::vector<String> > unpacked;
  std::vector<String> single_block;
  int lborder;
  bool is_range = false;
  std::stringstream ss;
  const char* s = expression.begin();
  const char* e = expression.end();
  const char* l = expression.begin();

  if (*s == '(') {
    ++s;
    ++l;
    // multiple blocks
    while (s != e) {
      if (*s == ',') {
        StoreExpression(l, s, is_range, lborder, single_block);
        l = s+1;
      } else if (*s == '-') {
        StoreRange(l, s, is_range, lborder, single_block);
      } else if (*s == '(') {
        ++l;
      } else if (*s == ')') {
        StoreExpression(l, s, is_range, lborder, single_block);
        l = s+1;
        if (! single_block.empty()) {
          unpacked.push_back(single_block);
        }
        single_block.clear();
      }
      ++s;
    }
  } else {
    // single block
    while (s != e) {
      if (*s == ',') {
        StoreExpression(l, s, is_range, lborder, single_block);
        l = s+1;
      } else if (*s == '-') {
        StoreRange(l, s, is_range, lborder, single_block);
      }
      ++s;
    }
    StoreExpression(l, s, is_range, lborder, single_block);

    if (is_range) {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_WARNING,
                                               "pdbx_struct_assembly_gen.oper_expression is missing a right border for a range expression.",
                                               this->GetCurrentLinenum()));
    }
    unpacked.push_back(single_block);
  }

  return unpacked;
}

void MMCifReader::ParsePdbxStructAssemblyGen(const std::vector<StringRef>& columns)
{
  MMCifBioUAssembly assembly;

  assembly.biounit_id = columns[indices_[ASSEMBLY_ID]].str();

  std::vector<StringRef> tmp_chains=columns[indices_[ASYM_ID_LIST]].split(',');
  std::vector<StringRef>::const_iterator tc_it;
  for (tc_it = tmp_chains.begin(); tc_it != tmp_chains.end(); ++tc_it) {
    assembly.chains.push_back(tc_it->str());
  }

  assembly.operations =
    this->UnPackOperExperession(columns[indices_[OPER_EXPRESSION]]);

  bu_assemblies_.push_back(assembly);
}

void MMCifReader::ParsePdbxStructOperList(const std::vector<StringRef>& columns)
{
  MMCifInfoTransOpPtr op(new MMCifInfoTransOp);

  op->SetID(columns[indices_[PSOL_ID]].str());
  op->SetType(columns[indices_[PSOL_TYPE]].str());

  if ((indices_[VECTOR_1] != -1)&&
      (indices_[VECTOR_2] != -1)&&
      (indices_[VECTOR_3] != -1)) {
    op->SetVector(this->TryGetReal(columns[indices_[VECTOR_1]],
                                   "pdbx_struct_oper_list.vector[1]"),
                  this->TryGetReal(columns[indices_[VECTOR_2]],
                                   "pdbx_struct_oper_list.vector[2]"),
                  this->TryGetReal(columns[indices_[VECTOR_3]],
                                   "pdbx_struct_oper_list.vector[3]"));
  }

  if ((indices_[MATRIX_1_1] != -1)&&
      (indices_[MATRIX_1_2] != -1)&&
      (indices_[MATRIX_1_3] != -1)&&
      (indices_[MATRIX_2_1] != -1)&&
      (indices_[MATRIX_2_2] != -1)&&
      (indices_[MATRIX_2_3] != -1)&&
      (indices_[MATRIX_3_1] != -1)&&
      (indices_[MATRIX_3_2] != -1)&&
      (indices_[MATRIX_3_3] != -1)) {
    op->SetMatrix(this->TryGetReal(columns[indices_[MATRIX_1_1]],
                                   "pdbx_struct_oper_list.matrix[1][1]"),
                  this->TryGetReal(columns[indices_[MATRIX_1_2]],
                                   "pdbx_struct_oper_list.matrix[1][2]"),
                  this->TryGetReal(columns[indices_[MATRIX_1_3]],
                                   "pdbx_struct_oper_list.matrix[1][3]"),
                  this->TryGetReal(columns[indices_[MATRIX_2_1]],
                                   "pdbx_struct_oper_list.matrix[2][1]"),
                  this->TryGetReal(columns[indices_[MATRIX_2_2]],
                                   "pdbx_struct_oper_list.matrix[2][2]"),
                  this->TryGetReal(columns[indices_[MATRIX_2_3]],
                                   "pdbx_struct_oper_list.matrix[2][3]"),
                  this->TryGetReal(columns[indices_[MATRIX_3_1]],
                                   "pdbx_struct_oper_list.matrix[3][1]"),
                  this->TryGetReal(columns[indices_[MATRIX_3_2]],
                                   "pdbx_struct_oper_list.matrix[3][2]"),
                  this->TryGetReal(columns[indices_[MATRIX_3_3]],
                                   "pdbx_struct_oper_list.matrix[3][3]"));
  }

  info_.AddOperation(op);
}

void MMCifReader::ParseStruct(const std::vector<StringRef>& columns)
{
  MMCifInfoStructDetails details = MMCifInfoStructDetails();

  details.SetEntryID(columns[indices_[STRUCT_ENTRY_ID]].str());

  if (indices_[STRUCT_TITLE] != -1) {
    details.SetTitle(columns[indices_[STRUCT_TITLE]].str());
  }

  if ((indices_[PDBX_CASP_FLAG] != -1) &&
      (columns[indices_[PDBX_CASP_FLAG]][0] != '?')) {
    details.SetCASPFlag(columns[indices_[PDBX_CASP_FLAG]][0]);
  }

  if (indices_[PDBX_DESCRIPTOR] != -1) {
    details.SetDescriptor(columns[indices_[PDBX_DESCRIPTOR]].str());
  }

  if (indices_[PDBX_FORMULA_WEIGHT] != -1) {
    if (!is_undef(columns[indices_[PDBX_FORMULA_WEIGHT]])) {
      details.SetMass(this->TryGetReal(columns[indices_[PDBX_FORMULA_WEIGHT]],
                                       "struct.pdbx_formula_weight"));
    }
  }

  if (indices_[PDBX_FORMULA_WEIGHT_METHOD] != -1) {
    details.SetMassMethod(columns[indices_[PDBX_FORMULA_WEIGHT_METHOD]].str());
  }

  if ((indices_[PDBX_MODEL_DETAILS] != -1) &&
      (columns[indices_[PDBX_MODEL_DETAILS]][0] != '?')) {
    details.SetModelDetails(columns[indices_[PDBX_MODEL_DETAILS]].str());
  }

  if ((indices_[PDBX_MODEL_TYPE_DETAILS] != -1) &&
      (columns[indices_[PDBX_MODEL_TYPE_DETAILS]][0] != '?')) {
    details.SetModelTypeDetails(
                              columns[indices_[PDBX_MODEL_TYPE_DETAILS]].str());
  }

  info_.SetStructDetails(details);
}

MMCifReader::MMCifSecStructElement MMCifReader::DetermineSecStructType(
                                                    const StringRef& type) const
{
  if (type == StringRef("HELX_P", 6)) {
    return MMCIF_HELIX;
  } else if (type == StringRef("HELX_OT_P", 9)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_P", 9)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_OT_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_AL_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_GA_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_OM_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_PI_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_27_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_3T_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_PP_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_P",     9)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_OT_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_AL_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_GA_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_OM_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_PI_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_27_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_3T_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_PP_P", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_N", 6)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_OT_N", 9)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_N", 9)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_OT_N", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_A_N", 11)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_B_N", 11)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_RH_Z_N", 11)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_N", 9)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_OT_N", 12)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_A_N", 11)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_B_N", 11)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("HELX_LH_Z_N", 11)) {
    return MMCIF_HELIX;
  }
  else if (type == StringRef("TURN_P", 6)) {
    return MMCIF_TURN;
  }
  else if (type == StringRef("TURN_OT_P", 9)) {
    return MMCIF_TURN;
  }
  else if (type == StringRef("TURN_TY1_P", 10)) {
    return MMCIF_TURN;
  }
  else if (type == StringRef("TURN_TY1P_P", 11)) {
    return MMCIF_TURN;
  }
  else if (type == StringRef("TURN_TY2_P", 10)) {
    return MMCIF_TURN;
  }
  else if (type == StringRef("TURN_TY2P_P", 11)) {
    return MMCIF_TURN;
  }
  else if (type == StringRef("TURN_TY3_P", 10)) {
    return MMCIF_TURN;
  }
  else if (type == StringRef("TURN_TY3P_P", 11)) {
    return MMCIF_TURN;
  }
  else if (type == StringRef("STRN", 4)) {
    return MMCIF_STRAND;
  }

  throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                    "Unknown secondary structure class found: "+
                                           type.str(),
                                           this->GetCurrentLinenum()));
}

void MMCifReader::ParseStructConf(const std::vector<StringRef>& columns)
{
  StringRef chain_name;
  int s_res_num;
  int e_res_num;

  // fetch chain name, first
  if(auth_chain_id_) {
    if (indices_[SC_BEG_AUTH_ASYM_ID] != -1) {
      chain_name = columns[indices_[SC_BEG_AUTH_ASYM_ID]];
    } else {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
"Chain name by author requested but 'struct_conf.beg_auth_asym_id' is not set.",
                                               this->GetCurrentLinenum()));
    }
  } else {
    chain_name = columns[indices_[SC_BEG_LABEL_ASYM_ID]];
  }

  if (restrict_chains_.size() == 0 ||
      restrict_chains_.find(chain_name.str()) != String::npos) {
    // fetch start and end
    std::pair<bool,int> s_beg = this->TryGetInt(columns[indices_[SC_BEG_LABEL_SEQ_ID]], // s_res_num
                                                "struct_conf.beg_label_seq_id",
                                                profile_.fault_tolerant);
    std::pair<bool,int> s_end = this->TryGetInt(columns[indices_[SC_END_LABEL_SEQ_ID]], // e_res_num
                                                "struct_conf.end_label_seq_id",
                                                profile_.fault_tolerant);
    if (!(s_beg.first && s_end.first)) {
      return;
    }
    s_res_num = s_beg.second;
    e_res_num = s_end.second;

    MMCifHSEntry hse = {to_res_num(s_res_num, ' '),
                        to_res_num(e_res_num, ' '),
                        chain_name.str()};
    
    MMCifSecStructElement type =
      DetermineSecStructType(columns[indices_[SC_CONF_TYPE_ID]]);
    if (type == MMCIF_HELIX) {
      helix_list_.push_back(hse);
    } else if (type == MMCIF_STRAND) {
      strand_list_.push_back(hse);
    }
  }
}

void MMCifReader::ParseStructSheetRange(const std::vector<StringRef>& columns)
{
  StringRef chain_name;
  int s_res_num;
  int e_res_num;

  if(auth_chain_id_) {
    if (indices_[SSR_BEG_AUTH_ASYM_ID] != -1) {
      chain_name = columns[indices_[SSR_BEG_AUTH_ASYM_ID]];
    } else {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
"Chain name by author requested but 'struct_sheet_range.beg_auth_asym_id' is not set.",
                                               this->GetCurrentLinenum()));
    }
  } else {
    chain_name = columns[indices_[SSR_BEG_LABEL_ASYM_ID]];
  }

  // restrict_chains feature not unit tested, since its about to be changed in
  // the future
  if (restrict_chains_.size() == 0 ||
      restrict_chains_.find(chain_name.str()) != String::npos) {

    std::pair<bool,int> s_beg = this->TryGetInt(columns[indices_[SSR_BEG_LABEL_SEQ_ID]], // s_res_num
                                                "struct_sheet_range.beg_label_seq_id",
                                                profile_.fault_tolerant);
    std::pair<bool,int> s_end = this->TryGetInt(columns[indices_[SSR_END_LABEL_SEQ_ID]], // e_res_num
                                                "struct_sheet_range.end_label_seq_id",
                                                profile_.fault_tolerant);
    if (!(s_beg.first && s_end.first)) {
      return;
    }
    s_res_num = s_beg.second;
    e_res_num = s_end.second;

    MMCifHSEntry hse = {to_res_num(s_res_num, ' '),
                        to_res_num(e_res_num, ' '),
                        chain_name.str()};
    strand_list_.push_back(hse);
  }
}

void MMCifReader::ParsePdbxDatabasePdbObsSpr(const std::vector<StringRef>&
                                             columns)
{
  MMCifInfoObsolete obs_data = MMCifInfoObsolete();

  obs_data.SetDate(columns[indices_[DATE]].str());
  obs_data.SetID(columns[indices_[PDPOS_ID]]);
  obs_data.SetPDBID(columns[indices_[PDB_ID]].str());
  obs_data.SetReplacedPDBID(columns[indices_[REPLACE_PDB_ID]].str());

  info_.SetObsoleteInfo(obs_data);
}

void MMCifReader::ParseDatabasePDBRev(const std::vector<StringRef>& columns)
{
  int num;
  StringRef date;
  StringRef status;

  num = this->TryGetInt(columns[indices_[DPI_NUM]], "database_PDB_rev.num");
  if (indices_[DPI_DATE] != -1) {
    date = columns[indices_[DPI_DATE]];
  } else {
    date = StringRef("", 0);
  }
  if (indices_[DPI_DATE_ORIGINAL] != -1) {
    info_.SetRevisionsDateOriginal(columns[indices_[DPI_DATE_ORIGINAL]].str());
  }
  if (indices_[DPI_STATUS] != -1) {
    status = columns[indices_[DPI_STATUS]];
  } else {
    status = StringRef("", 0);
  }
  info_.AddRevision(num, date.str(), status.str());
  database_PDB_rev_added_ = true;
}

void MMCifReader::ParsePdbxAuditRevisionHistory(
                                       const std::vector<StringRef>& columns) {
  // get ordinal and date
  int num = this->TryGetInt(columns[indices_[PARH_ORDINAL]],
                            "pdbx_audit_revision_history.ordinal");
  StringRef date = columns[indices_[PARH_REVISION_DATE]];
  int major = this->TryGetInt(columns[indices_[PARH_MAJOR]],
                              "pdbx_audit_revision_history.major_revision");
  int minor = this->TryGetInt(columns[indices_[PARH_MINOR]],
                              "pdbx_audit_revision_history.minor_revision");
  // add to map
  revisions_.push_back(MMCifRevisionDesc(num, date.str(), major, minor));
}

void MMCifReader::ParsePdbxAuditRevisionDetails(
                                       const std::vector<StringRef>& columns) {
  // get ordinal
  int num = this->TryGetInt(columns[indices_[PARD_REVISION_ORDINAL]],
                            "pdbx_audit_revision_details.revision_ordinal");
  // add type to map if available
  if (indices_[PARD_TYPE] != -1) {
    StringRef type = columns[indices_[PARD_TYPE]];
    revision_types_[num] = type.str();
  }
}

void MMCifReader::ParsePdbxDatabaseStatus(
                                       const std::vector<StringRef>& columns) {
  const int idx = indices_[PDS_RECVD_INITIAL_DEPOSITION_DATE];
  if (idx != -1) {
    info_.SetRevisionsDateOriginal(columns[idx].str());
  }
}

void MMCifReader::OnDataRow(const StarLoopDesc& header, 
                            const std::vector<StringRef>& columns)
{
  switch(category_) {
  case ATOM_SITE:
    LOG_TRACE("processing atom_site entry");
    this->ParseAndAddAtom(columns);
    break;
  case ENTITY:
    LOG_TRACE("processing entity entry");
    this->ParseEntity(columns);
    break;
  case ENTITY_POLY:
    LOG_TRACE("processing entity_poly entry");
    this->ParseEntityPoly(columns);
    break;
  case CITATION:
    LOG_TRACE("processing citation entry");
    this->ParseCitation(columns);
    break;
  case CITATION_AUTHOR:
    LOG_TRACE("processing citation_author entry")
    this->ParseCitationAuthor(columns);
    break;
  case EXPTL:
    LOG_TRACE("processing exptl entry")
    this->ParseExptl(columns);
    break;
  case REFINE:
    LOG_TRACE("processing refine entry")
    this->ParseRefine(columns);
    break;
  case PDBX_STRUCT_ASSEMBLY:
    LOG_TRACE("processing pdbx_struct_assembly entry")
    this->ParsePdbxStructAssembly(columns);
    break;
  case PDBX_STRUCT_ASSEMBLY_GEN:
    LOG_TRACE("processing pdbx_struct_assembly_gen entry")
    this->ParsePdbxStructAssemblyGen(columns);
    break;
  case PDBX_STRUCT_OPER_LIST:
    LOG_TRACE("processing pdbx_struct_oper_list entry")
    this->ParsePdbxStructOperList(columns);
    break;
  case STRUCT:
    LOG_TRACE("processing struct entry")
    this->ParseStruct(columns);
    break;
  case STRUCT_CONF:
    LOG_TRACE("processing struct_conf entry")
    this->ParseStructConf(columns);
    break;
  case STRUCT_SHEET_RANGE:
    LOG_TRACE("processing struct_sheet_range entry")
    this->ParseStructSheetRange(columns);
    break;
  case PDBX_DATABASE_PDB_OBS_SPR:
    LOG_TRACE("processing pdbx_database_PDB_obs_spr entry")
    this->ParsePdbxDatabasePdbObsSpr(columns);
    break;
  case STRUCT_REF:
    LOG_TRACE("processing struct_ref entry");
    this->ParseStructRef(columns);
    break;
  case STRUCT_REF_SEQ:
    LOG_TRACE("processing struct_ref entry");
    this->ParseStructRefSeq(columns);
    break;
  case STRUCT_REF_SEQ_DIF:
    LOG_TRACE("processing struct_ref entry");
    this->ParseStructRefSeqDif(columns);
    break;
  case DATABASE_PDB_REV:
    LOG_TRACE("processing database_PDB_rev entry");
    this->ParseDatabasePDBRev(columns);
    break;
  case PDBX_AUDIT_REVISION_HISTORY:
    LOG_TRACE("processing pdbx_audit_revision_history entry");
    this->ParsePdbxAuditRevisionHistory(columns);
    break;
  case PDBX_AUDIT_REVISION_DETAILS:
    LOG_TRACE("processing pdbx_audit_revision_details entry");
    this->ParsePdbxAuditRevisionDetails(columns);
    break;
  case PDBX_DATABASE_STATUS:
    LOG_TRACE("processing pdbx_database_status entry");
    this->ParsePdbxDatabaseStatus(columns);
    break;
  case PDBX_ENTITY_BRANCH:
    LOG_TRACE("processing pdbx_entity_branch entry");
    this->ParsePdbxEntityBranch(columns);
    break;
  case PDBX_ENTITY_BRANCH_LINK:
    LOG_TRACE("processing pdbx_entity_branch_link entry");
    this->ParsePdbxEntityBranchLink(columns);
    break;
  default:
    throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                       "Uncatched category '"+ header.GetCategory() +"' found.",
                                             this->GetCurrentLinenum()));
    return;
  }
}

void MMCifReader::AssignSecStructure(mol::EntityHandle ent)
{
  for (MMCifHSVector::const_iterator i=helix_list_.begin(), e=helix_list_.end();
       i!=e; ++i) {
    mol::ChainHandle chain = ent.FindChain(i->chain_name);
    if (!chain.IsValid()) {
      LOG_INFO("ignoring helix record for unknown chain " + i->chain_name);
      continue;
    }
    mol::SecStructure alpha(mol::SecStructure::ALPHA_HELIX);
    chain.AssignSecondaryStructure(alpha, i->start, i->end);
  }

  for (MMCifHSVector::const_iterator i=strand_list_.begin(),
         e=strand_list_.end();
       i!=e; ++i) {
    mol::ChainHandle chain=ent.FindChain(i->chain_name);
    if (!chain.IsValid()) {
      LOG_INFO("ignoring strand record for unknown chain " + i->chain_name);
      continue;
    }
    mol::SecStructure extended(mol::SecStructure::EXTENDED);
    chain.AssignSecondaryStructure(extended, i->start, i->end);
  }
}

void MMCifReader::ParseStructRef(const std::vector<StringRef>& columns)
{
	String ent_id=columns[indices_[SR_ENTITY_ID]].str();
	String db_name=columns[indices_[SR_DB_NAME]].str();
	String db_code=columns[indices_[SR_DB_CODE]].str();
	String id=columns[indices_[SR_ID]].str();
	String db_access;
	if (indices_[SR_DB_ACCESS]!=-1) {
		db_access=columns[indices_[SR_DB_ACCESS]].str();
	}
	MMCifInfoStructRefPtr sr(new MMCifInfoStructRef(id, ent_id, db_name, 
				                                          db_code, db_access));
	struct_refs_.push_back(sr);
}

void MMCifReader::ParseStructRefSeq(const std::vector<StringRef>& columns)
{
 String aln_id=columns[indices_[SRS_ALIGN_ID]].str();
 String sr_id=columns[indices_[SRS_STRUCT_REF_ID]].str();
 String chain_name;
 if (indices_[SRS_PDBX_STRAND_ID]!=-1) {
 	 chain_name=columns[indices_[SRS_PDBX_STRAND_ID]].str();
 }
 std::pair<bool,int> dbbeg=this->TryGetInt(columns[indices_[SRS_DB_ALIGN_BEG]], 
 		                                        "_struct_ref_seq.db_align_beg",
 		                                        profile_.fault_tolerant);
 std::pair<bool,int> dbend=this->TryGetInt(columns[indices_[SRS_DB_ALIGN_END]], 
 		                                       "_struct_ref_seq.db_align_end",
 		                                       profile_.fault_tolerant);
 std::pair<bool,int> entbeg=this->TryGetInt(columns[indices_[SRS_ENT_ALIGN_BEG]], 
 		                                        "_struct_ref_seq.seq_align_beg",
 		                                        profile_.fault_tolerant);
 std::pair<bool,int> entend=this->TryGetInt(columns[indices_[SRS_ENT_ALIGN_END]], 
 		                                        "_struct_ref_seq.seq_align_END",
 		                                        profile_.fault_tolerant);
 if (!(dbbeg.first && dbend.first && entbeg.first && entend.first)) {
 	 return;
 }
 bool found=false;
 for (MMCifInfoStructRefs::iterator i=struct_refs_.begin(), 
 		  e=struct_refs_.end(); i!=e; ++i) { 
 	 if ((*i)->GetID()==sr_id) {
		 (*i)->AddAlignedSeq(aln_id, chain_name, entbeg.second, entend.second, 
		 		                 dbbeg.second, dbend.second);
		 found=true;
 	 	 break;
 	 }
 }
 if (!found) {
 	 if (profile_.fault_tolerant) {
 	 	 LOG_ERROR("struct_ref_seq.ref_id points to inexistent struct_ref '"
 	 	 		       << sr_id <<  "'");
 	 	 return;
 	 }
	 std::stringstream ss;
	 ss << "struct_ref_seq.ref_id points to inexistent struct_ref '";
	 ss << sr_id << "'";
	 throw IOException(ss.str());
 }
}

void MMCifReader::ParseStructRefSeqDif(const std::vector<StringRef>& columns)
{
  String aln_id=columns[indices_[SRSD_ALIGN_ID]].str();
  String db_rnum;
  if (!is_undef(columns[indices_[SRSD_DB_RNUM]])) {
    db_rnum=columns[indices_[SRSD_DB_RNUM]].str();
  }
  std::pair<bool,int> seq_rnum(true, -1);
  if (!is_undef(columns[indices_[SRSD_SEQ_RNUM]])) {
    seq_rnum=this->TryGetInt(columns[indices_[SRSD_SEQ_RNUM]],
                             "_struct_ref_seq_dif.seq_num",
                             profile_.fault_tolerant);
    
  }
  if (!seq_rnum.first) {
    return;
  }
  String details;
  if (indices_[SRSD_DETAILS]!=-1) {
	  details=columns[indices_[SRSD_DETAILS]].str();
  }
  bool found=false;
  for (MMCifInfoStructRefs::iterator i=struct_refs_.begin(), 
 		  e=struct_refs_.end(); i!=e; ++i) { 
 	 if (MMCifInfoStructRefSeqPtr s=(*i)->GetAlignedSeq(aln_id)) {
		 s->AddDif(seq_rnum.second, db_rnum, details); 
		 found=true;
 	 	 break;
 	 }
 }
 if (!found) {
 	 if (profile_.fault_tolerant) {
 	 	 LOG_ERROR("struct_ref_seq_dif.align_id points to inexistent "
 	 	 		       "struct_ref_seq '" << aln_id <<  "'");
 	 	 return;
 	 }
	 std::stringstream ss;
	 ss << "struct_ref_seq.ref_id points to inexistent struct_ref '";
	 ss << aln_id << "'";
	 throw IOException(ss.str());
 }
}

void MMCifReader::ParsePdbxEntityBranch(const std::vector<StringRef>& columns)
{
  // get entity/ descreption entry
  MMCifEntityDescMap::iterator edm_it =
    GetEntityDescMapIterator(columns[indices_[BR_ENTITY_ID]].str());

  // store type
  if (indices_[BR_ENTITY_TYPE] != -1) {
    edm_it->second.type = mol::ChainTypeFromString(columns[indices_[EP_TYPE]]);
  }
}

void MMCifReader::ParsePdbxEntityBranchLink(const std::vector<StringRef>& columns)
{
  MMCifPdbxEntityBranchLink link_pair;

  String entity_id(columns[indices_[BL_ENTITY_ID]].str());

  // list of entities -> pairs of info for link
  link_pair.res_num_1 =
    this->TryGetInt(columns[indices_[BL_ENTITY_BRANCH_LIST_NUM_1]],
                    "pdbx_entity_branch_link.entity_branch_list_num_1");
  link_pair.cmp_1 = columns[indices_[BL_COMP_ID_1]].str();
  link_pair.atm_nm_1 = columns[indices_[BL_ATOM_ID_1]].str();
  link_pair.res_num_2 =
    this->TryGetInt(columns[indices_[BL_ENTITY_BRANCH_LIST_NUM_2]],
                    "pdbx_entity_branch_link.entity_branch_list_num_2");
  link_pair.cmp_2 = columns[indices_[BL_COMP_ID_2]].str();
  link_pair.atm_nm_2 = columns[indices_[BL_ATOM_ID_2]].str();

  /*if (indices_[BL_ATOM_STEREO_CONFIG_1] != -1) {
    char A = *columns[indices_[BL_ATOM_STEREO_CONFIG_1]].begin();
  }*/
  // check stereo values to be N S R
  /*if (indices_[BL_ATOM_STEREO_CONFIG_2] != -1) {
  }*/
  // check value order
  if (indices_[BL_VALUE_ORDER] != -1) {
    link_pair.bond_order = MMCifValueOrderToOSTBondOrder(
                                              columns[indices_[BL_VALUE_ORDER]]);
  } else {
    link_pair.bond_order = 1;
  }

  std::pair<MMCifPdbxEntityBranchLinkMap::iterator, bool> rit;

  // check if element already exists
  MMCifPdbxEntityBranchLinkMap::iterator blm_it =
    entity_branch_link_map_.find(entity_id);

  // if the entity was not seen before, create it in the map
  if (blm_it == entity_branch_link_map_.end()) {
    rit = entity_branch_link_map_.insert(
                   MMCifPdbxEntityBranchLinkMap::value_type(entity_id,
                                      std::vector<MMCifPdbxEntityBranchLink>()));
    blm_it = rit.first;
  }

  // add the link pair
  blm_it->second.push_back(link_pair);
}

void MMCifReader::OnEndData()
{
  mol::XCSEditor editor=ent_handle_.EditXCS(mol::BUFFERED_EDIT);

  // process chain types
  std::vector<std::pair<mol::ChainHandle, String> >::const_iterator css;
  MMCifEntityDescMap::const_iterator edm_it;
  MMCifPdbxEntityBranchLinkMap::const_iterator blm_it;
  std::vector<MMCifPdbxEntityBranchLink>::const_iterator bl_it;
  String pdb_auth_chain_name;
  for (css = chain_id_pairs_.begin(); css != chain_id_pairs_.end(); ++css) {
    // chain description
    edm_it = entity_desc_map_.find(css->second);
    if (edm_it != entity_desc_map_.end()) {
      editor.SetChainType(css->first, edm_it->second.type);
      editor.SetChainDescription(css->first, edm_it->second.details);
      if (edm_it->second.seqres.length() > 0) {
        seqres_.AddSequence(seq::CreateSequence(css->first.GetName(),
                                                edm_it->second.seqres));
        pdb_auth_chain_name = css->first.GetStringProp("pdb_auth_chain_name");
        info_.AddMMCifPDBChainTr(css->first.GetName(), pdb_auth_chain_name);
        info_.AddPDBMMCifChainTr(pdb_auth_chain_name, css->first.GetName());
      } else if (edm_it->second.type!=mol::CHAINTYPE_WATER) {
        // mark everything that doesn't have SEQRES and isn't of type
        // water as ligand
        mol::ChainHandle chain=css->first;
        mol::ResidueHandleList residues=chain.GetResidueList();
        for (mol::ResidueHandleList::iterator 
             i=residues.begin(), e=residues.end(); i!=e; ++i) {
          (*i).SetIsLigand(true);
        }
      }
    } else {
      LOG_WARNING("No entity description found for atom_site.label_entity_id '"
                  << css->second << "'");
    }
    // find
    blm_it = entity_branch_link_map_.find(css->second);
    // store linker pair
    if (blm_it != entity_branch_link_map_.end()) {
      for (bl_it = blm_it->second.begin(); bl_it != blm_it->second.end();
           ++bl_it) {
        mol::ResidueHandle res1 = css->first.FindResidue(to_res_num(
                                                         bl_it->res_num_1, ' '));
        mol::ResidueHandle res2 = css->first.FindResidue(to_res_num(
                                                         bl_it->res_num_2, ' '));
        info_.AddEntityBranchLink(css->first.GetName(),
                                  res1.FindAtom(bl_it->atm_nm_1),
                                  res2.FindAtom(bl_it->atm_nm_2),
                                  bl_it->bond_order);
      }
    }
  }

  // process citations (couple with authors
  // iterate citations
  MMCifCitationAuthorMap::const_iterator atm_it;
  for (atm_it = authors_map_.begin(); atm_it != authors_map_.end(); ++atm_it) {
    info_.AddAuthorsToCitation(StringRef(atm_it->first.c_str(),
                                         atm_it->first.length()),
                               atm_it->second.second);
  }

  bool found;
  MMCifBioUAssemblyVector::iterator bua_it;
  std::vector<std::vector<String> >::const_iterator aol_it;
  std::vector<String>::const_iterator aob_it;
  std::vector<MMCifInfoTransOpPtr> operation_list;
  std::map<String, MMCifPSAEntry>::const_iterator buom_it;
  std::vector<MMCifInfoTransOpPtr> operations = info_.GetOperations();
  info_.SetStructRefs(struct_refs_);
  std::vector<MMCifInfoTransOpPtr>::const_iterator buop_it;
  MMCifInfoBioUnit biounit;
  for (bua_it = bu_assemblies_.begin();
       bua_it != bu_assemblies_.end();
       ++bua_it) {
    biounit = MMCifInfoBioUnit();
    // pair with pdbx_struct_assembly entry
    buom_it = bu_origin_map_.find(bua_it->biounit_id);
    if (buom_it == bu_origin_map_.end()) {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                               "No pdbx_struct_assembly.id '"+
                                               bua_it->biounit_id +
                         "' found as requested by pdbx_struct_assembly_gen.")); 
    }
    biounit.SetDetails(buom_it->second.details);
    biounit.SetMethodDetails(buom_it->second.method_details);
    biounit.SetID(buom_it->first);
    biounit.SetChainList(bua_it->chains);

    // pair with pdbx_struct_oper_list
    for (aol_it = bua_it->operations.begin();
         aol_it != bua_it->operations.end();
         ++aol_it) {
      operation_list.clear();
      for (aob_it = aol_it->begin(); aob_it != aol_it->end(); aob_it++) {
        found = false;
        for (buop_it = operations.begin();
             buop_it != operations.end();
             ++buop_it) {
          if ((*buop_it)->GetID() == *aob_it) {
            operation_list.push_back(*buop_it);
            found = true;
            break;
          }
        }
        if (!found) {
          throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                                "No pdbx_struct_oper_list.id '"+
                                                   *aob_it +
                          "' found as requested by pdbx_struct_assembly_gen."));
        }
      }
      biounit.AddOperations(operation_list);
    }
    info_.AddBioUnit(biounit);
  }
  bu_assemblies_.clear();

  // create secondary structure from struct_conf info
  this->AssignSecStructure(ent_handle_);

  // add revision history for new style mmCIFs (only if no old data there)
  if (!database_PDB_rev_added_) {
    std::vector<MMCifRevisionDesc>::const_iterator r_it;
    for (r_it = revisions_.begin(); r_it != revisions_.end(); ++r_it) {
      // look for status
      const int num = r_it->num;
      std::map<int, String>::const_iterator rt_it = revision_types_.find(num);
      if (rt_it != revision_types_.end()) {
        info_.AddRevision(num, r_it->date, rt_it->second, r_it->major,
                          r_it->minor);
      } else {
        info_.AddRevision(num, r_it->date, "?", r_it->major, r_it->minor);
      }
    }
  }

  LOG_INFO("imported "
           << chain_count_ << " chains, "
           << residue_count_ << " residues, "
           << atom_count_ << " atoms;"
           << helix_list_.size() << " helices and "
           << strand_list_.size() << " strands");
}

unsigned char MMCifValueOrderToOSTBondOrder(const StringRef value_order)
{
  if (value_order == StringRef("sing", 4)) {
      return 1;
  }
  if (value_order == StringRef("doub", 4)) {
      return 2;
  }
  if (value_order == StringRef("trip", 4)) {
      return 3;
  }
  LOG_WARNING("Non-covered bond order found: '" << value_order << "'");
  if (value_order == StringRef("arom", 4)) {
      return 4;
  }
  if (value_order == StringRef("delo", 4)) {
      return 5;
  }
  if (value_order == StringRef("pi", 2)) {
      return 6;
  }
  if (value_order == StringRef("poly", 4)) {
      return 7;
  }
  if (value_order == StringRef("quad", 4)) {
      return 8;
  }
  return 1;
}

String OSTBondOrderToMMCifValueOrder(const unsigned char bond_order)
{
  if (bond_order == 1) {
    return String("sing");
  }
  if (bond_order == 2) {
    return String("doub");
  }
  if (bond_order == 3) {
    return String("trip");
  }
  if (bond_order == 4) {
    return String("arom");
  }
  if (bond_order == 5) {
      return String("delo");
  }
  if (bond_order == 6) {
      return String("pi");
  }
  if (bond_order == 7) {
      return String("poly");
  }
  if (bond_order == 8) {
      return String("quad");
  }
  LOG_WARNING("Unknow bond order found: '" << (int)bond_order << "'");
  return String("");
}

}}
