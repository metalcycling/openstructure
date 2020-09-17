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
#ifndef OST_MMCIF_READER_HH
#define OST_MMCIF_READER_HH

#include <map>

#include <ost/geom/geom.hh>
#include <ost/seq/sequence_list.hh>
#include <ost/mol/residue_handle.hh>
#include <ost/mol/chain_type.hh>
#include <ost/conop/compound_lib.hh>
#include <ost/io/mol/io_profile.hh>
#include <ost/io/io_exception.hh>
#include <ost/io/mol/star_parser.hh>
#include <ost/io/mol/mmcif_info.hh>

namespace ost { namespace io {

/// \brief reader for the mmcif file format
/// 
/// \section mmcif_format mmcif format description/ coverage
/// 
/// mmcif is an instance of the \link StarParser STAR format\endlink to store
/// entries of the PDB. The following data categories should be covered by this
/// reader:
/// 
/// \li atom_site
/// \li entity
/// \li entity_poly
/// \li citation
/// \li citation_author
/// \li exptl
/// \li refine
/// \li pdbx_struct_assembly
/// \li pdbx_struct_assembly_gen
/// \li pdbx_struct_oper_list
/// \li struct
/// \li struct_conf
/// \li struct_sheet_range
/// \li pdbx_database_PDB_obs_spr
/// \li database_PDB_rev
/// \li pdbx_entity_branch
/// \li pdbx_entity_branch_link
class DLLEXPORT_OST_IO MMCifReader : public StarParser  {
public:
  /// \brief create a MMCifReader
  ///
  /// \param stream input stream
  /// \param ent_handle
  /// \param profile

  MMCifReader(std::istream& stream, mol::EntityHandle& ent_handle,
              const IOProfile& profile);

  /// \brief create a MMCifReader
  ///
  /// \param filename input file
  /// \param ent_handle
  /// \param profile
  MMCifReader(const String& filename, mol::EntityHandle& ent_handle,
              const IOProfile& profile);

  /// \brief Initialise the reader.
  ///
  void Init();

  /// \brief Set up a fresh instance
  void ClearState();

  /// \brief Set names of restricted chains for the reader.
  ///
  /// \param restrict_chains chain name
  void SetRestrictChains(const String& restrict_chains);

  /// \brief Toggle reading of canonical sequence residues
  ///        (entity_poly.pdbx_seq_one_letter_code_can instead of
  ///        entity_poly.pdbx_seq_one_letter_code). This flag is exclusive.
  ///
  /// \param flag True for reading canonical sequences.
  void SetReadCanonicalSeqRes(bool flag)
  {
    seqres_can_ = flag;
  }

  const String& GetRestrictChains() const
  {
    return restrict_chains_;
  }

  /// \brief Enable or disable reading of auth_chain_id instead aof label_chain
  /// id (default)
  ///
  /// \param id enable (true) or disable (false) reading of auth_chain_id.
  void SetAuthChainID(bool id)
  {
    auth_chain_id_ = id;
  }

  /// \brief check mmcif input to be read. Substitutional function for
  /// \link StarParser StarParser\endlink.
  ///
  /// \param data_name value of the data_ tag
  ///
  /// \return true, if the blockcode (PDB id) is valid, false otherwise
  virtual bool OnBeginData(const StringRef& data_name);

  /// \brief check if a current loop is to be parsed
  ///
  /// \param header categories of the upcoming loop block
  ///
  /// \return bool
  virtual bool OnBeginLoop(const StarLoopDesc& header); // tested

  /// \brief read a row of data
  ///
  /// \param header categories and items
  /// \param columns data
  virtual void OnDataRow(const StarLoopDesc& header, 
                         const std::vector<StringRef>& columns);

  /// \brief Finalise parsing.
  virtual void OnEndData();

  /// \brief Return sequences
  ///
  /// \return List of sequences
  seq::SequenceList GetSeqRes() const {
    return seqres_;
  }

  /// \brief Toggle reading of SEQRES
  ///
  /// \param flag True enables, False disables reading SEQRES
  void SetReadSeqRes(bool flag)
  {
    read_seqres_ = flag;
  }

  /// \brief Check if reading of SEQRES is enabled
  ///
  /// \return True if reading of SEQRES is enabled
  bool GetReadSeqRes() const
  {
    return read_seqres_;
  }

  /// \brief Get additional information of the mmCIF file.
  ///
  /// \return MMCitfInfo object
  const MMCifInfo& GetInfo() { return info_; }

protected:
  /// \brief Store an item index from loop header in preparation for reading a 
  ///        row. Throws an exception if the item does not exist.
  ///
  /// \param mapping position the item index is stored at
  /// \param item exact item name to fetch
  /// \param header loop header to pull index from
  void TryStoreIdx(const int mapping,
                   const String& item,
                   const StarLoopDesc& header)
  {
    indices_[mapping] = header.GetIndex(item);
    
    if (indices_[mapping] == -1) {
      throw IOException(this->FormatDiagnostic(STAR_DIAG_ERROR,
                                               "No item '" + item +
                                               "' found in '" +
                                               header.GetCategory()+
                                               "' header",
                                               this->GetCurrentLinenum()));
     }
  } // tested

  /// \brief fetch values identifying atoms
  ///
  /// \param[in]  columns data row
  /// \param[out] auth_chain_name takes atom_site.label_asym_id or, if
  ///             auth_chain_id_ is set, atom_site.auth_asym_id as a chain name
  /// \param[out] cif_chain_name
  /// \param[out] res_name fetches atom_site.label_comp_id
  /// \param[out] resnum gets atom_site.label_seq_id if available, consecutive
  ///             numbers, otherwise
  /// \param[out] valid_res_num shows if we have a valid residue number or if
  ///             we have to invent our own
  /// \param[out] atom_name corresponds to label_atom_id
  /// \param[out] alt_loc gets first letter of atom_site.label_alt_id
  bool ParseAtomIdent(const std::vector<StringRef>& columns,
                      String& auth_chain_name,
                      String& cif_chain_name,
                      StringRef& res_name,
                      mol::ResNum& resnum,
                      bool& valid_res_num,
                      StringRef& atom_name,
                      char& alt_loc);

  /// \brief Fetch atom information and store it.
  ///
  /// \param columns data row
  void ParseAndAddAtom(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF entity information
  ///
  /// \param columns data row
  void ParseEntity(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF entity_poly information
  ///
  /// \param columns data row
  void ParseEntityPoly(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF citation information
  ///
  /// \param columns data row
  void ParseCitation(const std::vector<StringRef>& columns);

	const MMCifInfoStructRefs& GetStructRefs() const { return struct_refs_; }
  /// \brief convert the seqres data item to canonical form. 
  /// 
  /// The seqres sequence lists non-standard residues in parenthesis. For 
  /// proper handling of our sequence classes, these need to be converted to 
  /// one-letter-codes. Ideally, we would use the canonical SEQRES. This is 
  /// not possible, however, since the PDB assigns multiple one letter codes 
  /// to some of the residues. To be consistent, we have to do the conversion
  /// on our own.
  String ConvertSEQRES(const String& seqres, conop::CompoundLibPtr compound_lib);
  /// \brief Fetch mmCIF citation_author information
  ///
  /// \param columns data row
  void ParseCitationAuthor(const std::vector<StringRef>& columns);
  
  /// \ brief parse a row in the struct_ref category
  void ParseStructRef(const std::vector<StringRef>& columns);
  
  /// \brief parse row in the struct_ref_seq category
  void ParseStructRefSeq(const std::vector<StringRef>& columns);

  ///  \brief parse row in the struct_ref_seq_dif category
  void ParseStructRefSeqDif(const std::vector<StringRef>& columns);
  /// \brief Fetch mmCIF exptl information
  ///
  /// \param columns data row
  void ParseExptl(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF refine information
  ///
  /// \param columns data row
  void ParseRefine(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF pdbx_struct_assembly information
  ///
  /// \param columns data row
  void ParsePdbxStructAssembly(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF pdbx_struct_assembly_gen information
  ///
  /// \param columns data row
  void ParsePdbxStructAssemblyGen(const std::vector<StringRef>& columns);

  std::vector<std::vector<String> > UnPackOperExperession(StringRef expression);

  void StoreExpression(const char* l, const char* s,
                       bool& is_range, int lborder,
                       std::vector<String>& single_block);

  void StoreRange(const char*& l, const char* s, bool& is_range, int& lborder,
                  std::vector<String>& single_block);

  /// \brief Fetch mmCIF pdbx_struct_oper_list information
  ///
  /// \param columns data row
  void ParsePdbxStructOperList(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF database_PDB_rev information
  ///
  /// \param columns data row
  void ParseDatabasePDBRev(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF pdbx_audit_revision_history information
  ///
  /// \param columns data row
  void ParsePdbxAuditRevisionHistory(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF pdbx_audit_revision_details information
  ///
  /// \param columns data row
  void ParsePdbxAuditRevisionDetails(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF pdbx_database_status information
  ///
  /// \param columns data row
  void ParsePdbxDatabaseStatus(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF struct information
  ///
  /// \param columns data row
  void ParseStruct(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF struct_conf (secondary structure) information
  ///
  /// \param columns data row
  void ParseStructConf(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF struct_sheet_range (beta sheets) information
  ///
  /// \param columns data row
  void ParseStructSheetRange(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF pdbx_database_PDB_obs_spr information
  ///
  /// \param columns data row
  void ParsePdbxDatabasePdbObsSpr(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF pdbx_entity_branch information
  ///
  /// \param columns data row
  void ParsePdbxEntityBranch(const std::vector<StringRef>& columns);

  /// \brief Fetch mmCIF pdbx_entity_branch_link information
  ///
  /// \param columns data row
  void ParsePdbxEntityBranchLink(const std::vector<StringRef>& columns);

  /// \struct types of secondary structure
  typedef enum {
    MMCIF_HELIX,
    MMCIF_STRAND,
    MMCIF_TURN
  } MMCifSecStructElement;

  /// \brief Check whether an element was classified sheet or helix
  ///
  /// \param type Type to be classified
  MMCifSecStructElement DetermineSecStructType(const StringRef& type) const;

  /// \brief Transform data from struct_conf entry into secondary structure
  ///
  /// \param ent Entity to assign secondary structure to
  void AssignSecStructure(mol::EntityHandle ent);

private:
  /// \enum magic numbers of this class
  typedef enum {
    MAX_ITEMS_IN_ROW=18 ///< count for possible items in a loop row
  } MMCifMagicNos;

  /// \enum items of the atom_site category
  typedef enum {
    AUTH_ASYM_ID,      ///< chain name by author as in PDB
    AS_ID,             ///< atom serial id
    LABEL_ALT_ID,      ///< AltLoc
    LABEL_ASYM_ID,     ///< chain name by PDB
    LABEL_ATOM_ID,
    LABEL_COMP_ID,
    LABEL_ENTITY_ID,   ///< link to category entity
    LABEL_SEQ_ID,      ///< residue no.
    AUTH_SEQ_ID,       ///< residue no. by author
    TYPE_SYMBOL,       ///< chemical element
    CARTN_X,           ///< Coordinates ||IMPORTANT: This 3 entries have to stay
    CARTN_Y,           ///< Coordinates ||together for the reader to work!
    CARTN_Z,           ///< Coordinates ||
    OCCUPANCY,
    B_ISO_OR_EQUIV,
    PDBX_PDB_INS_CODE,
    GROUP_PDB,         ///< record name
    PDBX_PDB_MODEL_NUM ///< model no. (especially NMR structures)
  } AtomSiteItems;

  /// \enum items of the entity category
  typedef enum {
    E_ID,              ///< unique identifier
    E_TYPE,            ///< polymer, non-polymer or water
    PDBX_DESCRIPTION   ///< special aspects of the entity
  } EntityItems;

  /// \enum items of the entity_poly category
  typedef enum {
    ENTITY_ID,                    ///< pointer to entity.id
    EP_TYPE,                      ///< type of polymer
    PDBX_SEQ_ONE_LETTER_CODE,     ///< sequence, 1-letter code
    PDBX_SEQ_ONE_LETTER_CODE_CAN  ///< canonical sequence, 1-letter code
  } EntityPolyItems;

  /// \enum items of the citation category
  typedef enum {
    CITATION_ID,                  ///< unique identifier
    ABSTRACT_ID_CAS,              ///< CAS identifier
    BOOK_ID_ISBN,                 ///< ISBN code assigned, if book cited
    BOOK_TITLE,                   ///< title of book storing the citation
    BOOK_PUBLISHER,               ///< name of publisher f a book
    BOOK_PUBLISHER_CITY,          ///< location of a publisher of a book
    JOURNAL_ABBREV,               ///< abbreviated journal title for articles
    JOURNAL_VOLUME,               ///< volume of cited journal
    PAGE_FIRST,                   ///< first page of citation
    PAGE_LAST,                    ///< last page of citation
    PDBX_DATABASE_ID_DOI,         ///< Document Object Identifier of doi.org
    PDBX_DATABASE_ID_PUBMED,      ///< Ascession number of PubMed
    YEAR,                         ///< year of the citation
    TITLE                         ///< title of the citation
  } CitationItems;

  /// \enum items of the citation_author category
  typedef enum {
    AUTHOR_CITATION_ID,           ///< link to CITATION_ID
    AUTHOR_NAME,                  ///< name of an author
    ORDINAL                       ///< position in author list
  } CitationAuthorItems;

  /// \enum items of the exptl category
  typedef enum {
    EXPTL_ENTRY_ID,               ///< identifier
    METHOD                        ///< method of the experiment
  } ExptlItems;

  /// \enum items of the refine category
  typedef enum {
    REFINE_ENTRY_ID,              ///< id
    LS_D_RES_HIGH,                ///< crystal resolution
    LS_D_RES_LOW,
    LS_R_FACTOR_R_WORK,           ///< R-work value
    LS_R_FACTOR_R_FREE            ///< R-free value
  } RefineItems;

  /// \enum items of the pdbx_struct_assembly category
  typedef enum {
    PSA_DETAILS,                  ///< special aspects of the assembly
    PSA_ID,                       ///< unique identifier
    METHOD_DETAILS                ///< details about assembly computation
  } PdbxStructAssemblyItems;

  // \enum items of the struct_ref category
  typedef enum {
  	SR_ENTITY_ID,
  	SR_ID,
  	SR_DB_CODE,
  	SR_DB_NAME,
  	SR_DB_ACCESS
	} StructRefItems;
	
  /// \enum items of the struct_ref_seq category
  typedef enum {
    SRS_ALIGN_ID,
    SRS_STRUCT_REF_ID,
    SRS_PDBX_STRAND_ID,
    SRS_DB_ALIGN_BEG,
    SRS_DB_ALIGN_END,
    SRS_ENT_ALIGN_BEG,
    SRS_ENT_ALIGN_END
  } StructRefSeqItems;
  
  /// \enum items of the struct_ref_seq_dif category
  typedef enum {
    SRSD_ALIGN_ID,
    SRSD_SEQ_RNUM,
    SRSD_DB_RNUM,
    SRSD_DETAILS
  } StructRefSeqDifItems;

  /// \enum items of the pdbx_struct_assembly_gen category
  typedef enum {
    ASSEMBLY_ID,                  ///< link to pdbx_struct_assembly.id
    ASYM_ID_LIST,                 ///< list of chains
    OPER_EXPRESSION               ///< list of pdbx_struct_oper_list.ids
  } PdbxStructAssemblyGenItems;

  /// \enum items of the pdbx_struct_oper_list category
  typedef enum {
    PSOL_ID,    ///< unique identifier
    PSOL_TYPE,  ///< type of operation
    VECTOR_1,   ///< vector component
    VECTOR_2,   ///< vector component
    VECTOR_3,   ///< vector component
    MATRIX_1_1, ///< matrix component
    MATRIX_1_2, ///< matrix component
    MATRIX_1_3, ///< matrix component
    MATRIX_2_1, ///< matrix component
    MATRIX_2_2, ///< matrix component
    MATRIX_2_3, ///< matrix component
    MATRIX_3_1, ///< matrix component
    MATRIX_3_2, ///< matrix component
    MATRIX_3_3  ///< matrix component
  } PdbxStructOperListItems;

  /// \enum items of the struct category
  typedef enum {
     STRUCT_ENTRY_ID,     ///< name of the structure
     PDBX_CASP_FLAG,      ///< CASP/ CASD-NMR/ other contest target
     PDBX_DESCRIPTOR,     ///< descriptor for NDB structure/ PDB COMPND record
     PDBX_FORMULA_WEIGHT, ///< mass in daltons
     PDBX_FORMULA_WEIGHT_METHOD, ///< how mass was determined
     PDBX_MODEL_DETAILS,  ///< methodology
     PDBX_MODEL_TYPE_DETAILS, ///< type of structure model
     STRUCT_TITLE         ///< title for the data block
  } StructItems;

  /// \enum items of the struct_conf category
  typedef enum {
    SC_BEG_AUTH_ASYM_ID,  ///< Starting residue (atom_site.auth_asym_id)
    SC_BEG_LABEL_ASYM_ID, ///< Starting residue (atom_site.label_asym_id)
    SC_BEG_LABEL_COMP_ID, ///< Starting residue (atom_site.label_comp_id)
    SC_BEG_LABEL_SEQ_ID,  ///< Starting residue (atom_site.label_seq_id)
    SC_CONF_TYPE_ID,      ///< Pointer to struct_conf_type.id
    SC_END_AUTH_ASYM_ID,  ///< Ending residue, points to atom_site.auth_asym_id
    SC_END_LABEL_ASYM_ID, ///< Ending residue, points to atom_site.label_asym_id
    SC_END_LABEL_COMP_ID, ///< Ending residue, points to atom_site.label_comp_id
    SC_END_LABEL_SEQ_ID,  ///< Ending residue, points to atom_site.label_seq_id
    SC_ID,                ///< Unique identifier
  } StructConfItems;

  /// \enum items of the struct_sheet_range category
  typedef enum {
    SSR_BEG_LABEL_ASYM_ID,     ///< start, chain name (atom_site.label_asym_id)
    SSR_BEG_LABEL_COMP_ID,     ///< start, atom_site.label_comp_id
    SSR_BEG_LABEL_SEQ_ID,    ///< start, residue number (atom_site.label_seq_id)
    SSR_END_LABEL_ASYM_ID,     ///< end, chain name (atom_site.label_asym_id)
    SSR_END_LABEL_COMP_ID,     ///< end, atom_site.label_comp_id
    SSR_END_LABEL_SEQ_ID,      ///< end, residue number (atom_site.label_seq_id)
    SSR_SHEET_ID,              ///< unique identifier
    SSR_ID,                    ///< link to struct_sheet.id
    SSR_BEG_AUTH_ASYM_ID,      ///< alternative start, (atom_site.auth_asym_id)
    SSR_END_AUTH_ASYM_ID,      ///< alternative end, (atom_site.auth_asym_id)
    } StructSheetRangeItems;

  /// \enum items of the pdbx_database_PDB_obs_spr category
  typedef enum {
    DATE,           ///< date of replacement
    PDPOS_ID,       ///< type of obsolete of this entry
    PDB_ID,         ///< NEW PDB ID
    REPLACE_PDB_ID, ///< OLD PDB ID
  } PdbxDatabasePDBObsSpr;

  /// \enum categories of the database_PDB_rev category
  typedef enum {
    DPI_NUM,           ///< unique identifier
    DPI_DATE,          ///< revision date
    DPI_DATE_ORIGINAL, ///< date of first sight
    DPI_STATUS,        ///< status of a revision
  } DatabasePDBRevItems;

  /// \enum categories of the pdbx_audit_revision_history category
  typedef enum {
    PARH_ORDINAL,         ///< unique identifier
    PARH_REVISION_DATE,   ///< revision date
    PARH_MAJOR,           ///< major version
    PARH_MINOR,           ///< minor version
  } PdbxAuditRevisionHistoryItems;

  /// \enum categories of the pdbx_audit_revision_details category
  typedef enum {
    PARD_REVISION_ORDINAL, ///< link to pdbx_audit_revision_history.ordinal
    PARD_TYPE,             ///< type classification of the revision
  } PdbxAuditRevisionDetailsItems;

  /// \enum categories of the pdbx_database_status category
  typedef enum {
    PDS_RECVD_INITIAL_DEPOSITION_DATE, ///< date of initial deposition
  } PdbxDatabaseStatusItems;

  /// \enum items of the pdbx_entity_branch category (pendant to entity_poly)
  typedef enum {
    BR_ENTITY_ID,  ///< pointer to entity.id
    BR_ENTITY_TYPE ///< type of branched molecular entity
  } EntityBranchItems;

  /// \enum items of the pdbx_entity_branch_link category
  typedef enum {
    BL_ENTITY_ID,                /// < pointer to entity.id
    BL_ATOM_ID_1,                /// < atom identifier (element + number)
    BL_ATOM_ID_2,                /// < atom identifier (element + number)
    BL_COMP_ID_1,                /// < tlc of component
    BL_COMP_ID_2,                /// < tlc of component
    BL_ENTITY_BRANCH_LIST_NUM_1, /// < res. no. (pdbx_entity_branch_list.num)
    BL_ENTITY_BRANCH_LIST_NUM_2, /// < res. no. (pdbx_entity_branch_list.num)
    BL_ATOM_STEREO_CONFIG_1,     /// < chiral configuration (R/ S)
    BL_ATOM_STEREO_CONFIG_2,     /// < chiral configuration (R/ S)
    BL_VALUE_ORDER               /// < bond order
  } EntityBranchLinkItems;

  /// \enum categories of the mmcif format
  typedef enum {
    ATOM_SITE,
    ENTITY,
    ENTITY_POLY,
    CITATION,
    CITATION_AUTHOR,
    EXPTL,
    REFINE,
    PDBX_STRUCT_ASSEMBLY,
    PDBX_STRUCT_ASSEMBLY_GEN,
    PDBX_STRUCT_OPER_LIST,
    STRUCT,
    STRUCT_CONF,
    STRUCT_SHEET_RANGE,
    PDBX_DATABASE_PDB_OBS_SPR,
    STRUCT_REF,
    STRUCT_REF_SEQ,
    STRUCT_REF_SEQ_DIF,
    DATABASE_PDB_REV,
    PDBX_AUDIT_REVISION_HISTORY,
    PDBX_AUDIT_REVISION_DETAILS,
    PDBX_DATABASE_STATUS,
    PDBX_ENTITY_BRANCH,
    PDBX_ENTITY_BRANCH_LINK,
    DONT_KNOW
  } MMCifCategory;

  /// \struct keeping track of entity information
  typedef struct {
    mol::ChainType type; ///< characterise entity
    String details;      ///< description of this entity
    String seqres;       ///< chain of monomers
  } MMCifEntityDesc;
  typedef std::map<String, MMCifEntityDesc> MMCifEntityDescMap;

  /// \brief Get an iterator for MMCifEntityDescMap by finding an element or
  ///        inserting a new one into the map.
  /// \param entity_id ID of the entity to talk to
  MMCifEntityDescMap::iterator GetEntityDescMapIterator(const String& entity_id);

  /// \struct assembly information
  typedef struct {
    String biounit_id;                              ///< identifier for the bu
    std::vector<String> chains;                     ///< chains affected by
                                                    /// this operations
    std::vector<std::vector<String> > operations;   ///< list of links to
                                                    /// MMCifBioUOperation
  } MMCifBioUAssembly;
  typedef std::vector<MMCifBioUAssembly> MMCifBioUAssemblyVector;

  typedef std::map<String, std::pair<std::vector<int>, std::vector<String> > >
    MMCifCitationAuthorMap;

  /// \struct store struct_conf info (secondary structure)
  typedef struct {
    mol::ResNum start;
    mol::ResNum end;
    String chain_name;
  } MMCifHSEntry;
  typedef std::vector<MMCifHSEntry> MMCifHSVector;

  /// \struct store pdbx_struct_assembly info
  typedef struct {
    String details;
    String method_details;
  } MMCifPSAEntry;
  typedef std::map<String, MMCifPSAEntry> MMCifPSAMap;

  /// \struct keeping track of revision information
  struct MMCifRevisionDesc {
    // silly GCC note: major() & minor() exist as macros...facepalm
    MMCifRevisionDesc(int _num, const String& _date, int _major, int _minor)
                     : date(_date) {
      num = _num;
      major = _major;
      minor = _minor;
    }
    int num;      ///< unique identifier
    String date;  ///< revision date
    int major;    ///< major version
    int minor;    ///< minor version
  };

  /// \struct to keep entity_branch linker pairs while parsing
  typedef struct {
    int res_num_1;    ///< index of the linked residue relative to its chain
    String cmp_1;     ///< tlc of residue (sugar)
    String atm_nm_1;  ///< (IUPAC) name of the linking atom
    int res_num_2;    ///< index of the linked residue relative to its chain
    String cmp_2;     ///< tlc of residue (sugar)
    String atm_nm_2;  ///< index of the linked residue relative to its chain
    unsigned char bond_order; ///< ID of value_order as OST bond_order
  } MMCifPdbxEntityBranchLink;
  typedef std::map<String, std::vector<MMCifPdbxEntityBranchLink> >
    MMCifPdbxEntityBranchLinkMap;

  // members
  MMCifCategory category_;
  int category_counts_[DONT_KNOW+1]; ///< overall no. of atom_site loops
  int indices_[MAX_ITEMS_IN_ROW]; ///< map items to values in loops
  const IOProfile& profile_;
  mol::EntityHandle& ent_handle_;
  String restrict_chains_;
  bool auth_chain_id_;       ///< use chain IDs given by authors rather than pdb
  bool seqres_can_;          ///< read canonical 1-letter residues?
  mol::ChainHandle curr_chain_;
  mol::ResidueHandle curr_residue_;
  int chain_count_;
  int residue_count_;
  int atom_count_;
  bool warned_name_mismatch_;
  bool warned_rule_based_;
  String subst_res_id_; ///< work around for missing label_seq_id's
  bool has_model_;      ///< keep track of models through different atom_sites
  int curr_model_;      ///< if we have pdbx_PDB_model_num, store no.
  std::vector<std::pair<mol::ChainHandle, String> > chain_id_pairs_;
  ///< chain and label_entity_id
  MMCifEntityDescMap entity_desc_map_; ///< stores entity items
  seq::SequenceList seqres_;
  bool read_seqres_;
  MMCifInfo info_;      ///< info container
  MMCifCitationAuthorMap authors_map_;
  MMCifBioUAssemblyVector bu_assemblies_;
  MMCifPSAMap bu_origin_map_; ///< pdbx_struct_assembly
  MMCifHSVector helix_list_; ///< for storing struct_conf sec.struct. data
  MMCifHSVector strand_list_; ///< for storing struct_conf sec.struct. data
  MMCifInfoStructRefs struct_refs_;
  // for storing revisions 
  std::vector<MMCifRevisionDesc> revisions_;
  std::map<int, String> revision_types_;
  bool database_PDB_rev_added_;
  // for entity_branch connections
  MMCifPdbxEntityBranchLinkMap entity_branch_link_map_;
};

/// \brief Translate mmCIF info on bond type (e.g.
///        pdbx_entity_branch_link.value_order) to OST bond_order
///
/// \param value_order abbreviation detemrining the bond order
DLLEXPORT_OST_IO unsigned char MMCifValueOrderToOSTBondOrder(
                                                    const StringRef value_order);

/// \brief Translate an OST bond_order to mmCIF value_order
///
/// \param bond_order OST bond order
DLLEXPORT_OST_IO String OSTBondOrderToMMCifValueOrder(
                                                 const unsigned char bond_order);
}}

#endif

//  LocalWords:  MMCifEntityDescMap
