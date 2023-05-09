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

#include "USalign/TMalign.h" // include for the external TMalign
#include "USalign/MMalign.h"

#include <ost/mol/atom_view.hh>
#include <ost/message.hh>
#include <ost/bindings/wrap_tmalign.hh>

namespace ost{ namespace bindings{


// hacked version of MMalign_final which does not print anything to stdout
// but returns an MMAlignResult object

MMAlignResult MMalign_final_hacked(
    const string xname, const string yname,
    const vector<string> chainID_list1, const vector<string> chainID_list2,
    string fname_super, string fname_lign, string fname_matrix,
    const vector<vector<vector<double> > >&xa_vec,
    const vector<vector<vector<double> > >&ya_vec,
    const vector<vector<char> >&seqx_vec, const vector<vector<char> >&seqy_vec,
    const vector<vector<char> >&secx_vec, const vector<vector<char> >&secy_vec,
    const vector<int> &mol_vec1, const vector<int> &mol_vec2,
    const vector<int> &xlen_vec, const vector<int> &ylen_vec,
    double **xa, double **ya, char *seqx, char *seqy, char *secx, char *secy,
    int len_aa, int len_na, int chain1_num, int chain2_num,
    double **TMave_mat,
    vector<vector<string> >&seqxA_mat, vector<vector<string> >&seqM_mat,
    vector<vector<string> >&seqyA_mat, int *assign1_list, int *assign2_list,
    vector<string>&sequence, const double d0_scale, const bool m_opt,
    const int o_opt, const int outfmt_opt, const int ter_opt,
    const int split_opt, const bool a_opt, const bool d_opt,
    const bool fast_opt, const bool full_opt, const int mirror_opt,
    const vector<string>&resi_vec1, const vector<string>&resi_vec2)
{
    int i,j;
    int xlen=0;
    int ylen=0;
    for (i=0;i<chain1_num;i++) xlen+=xlen_vec[i];
    for (j=0;j<chain2_num;j++) ylen+=ylen_vec[j];
    //if (xlen<=3 || ylen<=3) return;
    if (xlen<=3 || ylen<=3) return MMAlignResult();

    seqx = new char[xlen+1];
    secx = new char[xlen+1];
    NewArray(&xa, xlen, 3);
    seqy = new char[ylen+1];
    secy = new char[ylen+1];
    NewArray(&ya, ylen, 3);

    int mol_type=copy_chain_pair_data(xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy, chain1_num, chain2_num,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence);

    /* declare variable specific to this pair of TMalign */
    double t0[3], u0[3][3];
    double TM1, TM2;
    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
    double d0_0, TM_0;
    double d0A, d0B, d0u, d0a;
    double d0_out=5.0;
    string seqM, seqxA, seqyA;// for output alignment
    double rmsd0 = 0.0;
    int L_ali;                // Aligned length in standard_TMscore
    double Liden=0;
    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
    int n_ali=0;
    int n_ali8=0;

    double Lnorm_ass=len_aa+len_na;

    /* entry function for structure alignment */
    TMalign_main(xa, ya, seqx, seqy, secx, secy,
        t0, u0, TM1, TM2, TM3, TM4, TM5,
        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
        xlen, ylen, sequence, Lnorm_ass, d0_scale,
        3, a_opt, false, d_opt, fast_opt, mol_type, -1);

    // this whole stuff is not needed in the hacked version

    /* prepare full complex alignment */
    //string chainID1="";
    //string chainID2="";
    //sequence.clear();
    //sequence.push_back(""); // seqxA
    //sequence.push_back(""); // seqyA
    //sequence.push_back(""); // seqM
    //int aln_start=0;
    //int aln_end=0;
    //for (i=0;i<chain1_num;i++)
    //{
    //    j=assign1_list[i];
    //    if (j<0) continue;
    //    chainID1+=chainID_list1[i];
    //    chainID2+=chainID_list2[j];
    //    sequence[0]+=seqxA_mat[i][j]+'*';
    //    sequence[1]+=seqyA_mat[i][j]+'*';

    //    aln_end+=seqxA_mat[i][j].size();
    //    seqM_mat[i][j]=seqM.substr(aln_start,aln_end-aln_start);
    //    sequence[2]+=seqM_mat[i][j]+'*';
    //    aln_start=aln_end;
    //}

    /* prepare unaligned region */
    //for (i=0;i<chain1_num;i++)
    //{
    //    if (assign1_list[i]>=0) continue;
    //    chainID1+=chainID_list1[i];
    //    chainID2+=':';
    //    string s(seqx_vec[i].begin(),seqx_vec[i].end());
    //    sequence[0]+=s.substr(0,xlen_vec[i])+'*';
    //    sequence[1]+=string(xlen_vec[i],'-')+'*';
    //    s.clear();
    //    sequence[2]+=string(xlen_vec[i],' ')+'*';
    //}
    //for (j=0;j<chain2_num;j++)
    //{
    //    if (assign2_list[j]>=0) continue;
    //    chainID1+=':';
    //    chainID2+=chainID_list2[j];
    //    string s(seqy_vec[j].begin(),seqy_vec[j].end());
    //    sequence[0]+=string(ylen_vec[j],'-')+'*';
    //    sequence[1]+=s.substr(0,ylen_vec[j])+'*';
    //    s.clear();
    //    sequence[2]+=string(ylen_vec[j],' ')+'*';
    //}

    /* print alignment */
    //output_results(xname, yname, chainID1.c_str(), chainID2.c_str(),
    //    xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
    //    sequence[2].c_str(), sequence[0].c_str(), sequence[1].c_str(),
    //    Liden, n_ali8, L_ali, TM_ali, rmsd_ali,
    //    TM_0, d0_0, d0A, d0B, 0, d0_scale, d0a, d0u, 
    //    (m_opt?fname_matrix:"").c_str(), outfmt_opt, ter_opt, true,
    //    split_opt, o_opt, fname_super,
    //    false, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);


    MMAlignResult res;
    res.rmsd = rmsd0;
    res.tm_score = TM1;
    res.transform = geom::Mat4(u0[0][0], u0[0][1], u0[0][2], t0[0],
                               u0[1][0], u0[1][1], u0[1][2], t0[1],
                               u0[2][0], u0[2][1], u0[2][2], t0[2],
                               0.0, 0.0, 0.0, 1.0);

    for (i=0;i<chain1_num;i++)
    {
        j=assign1_list[i];
        if (j<0) continue;
        res.ent1_mapped_chains.push_back(chainID_list1[i]);
        res.ent2_mapped_chains.push_back(chainID_list2[j]);
        ost::seq::AlignmentHandle aln = ost::seq::CreateAlignment();
        aln.AddSequence(ost::seq::CreateSequence(chainID_list1[i],
                                                 seqxA_mat[i][j]));
        aln.AddSequence(ost::seq::CreateSequence(chainID_list2[j],
                                                 seqyA_mat[i][j]));
        res.alignments.push_back(aln);
    }

    /* clean up */
    seqM.clear();
    seqxA.clear();
    seqyA.clear();
    delete [] seqx;
    delete [] seqy;
    delete [] secx;
    delete [] secy;
    DeleteArray(&xa,xlen);
    DeleteArray(&ya,ylen);
    // sequence gets never filled in hacked version
    //sequence[0].clear();
    //sequence[1].clear();
    //sequence[2].clear();

    return res;

    // deleted remainder that is only active if full_opt is enabled in original
    // code

}

void parse_chain_list_hacked(const std::vector<geom::Vec3List>& pos,
                             const ost::seq::SequenceList& seq,
                             const std::vector<bool>& rna,
                             const std::vector<string>&chain_list,
                             std::vector<std::vector<std::vector<double> > >&a_vec,
                             std::vector<std::vector<char> >&seq_vec,
                             std::vector<std::vector<char> >&sec_vec,
                             std::vector<int>&mol_vec, vector<int>&len_vec,
                             std::vector<std::string>&chainID_list,
                             const int ter_opt, const int split_opt,
                             const string mol_opt, const int infmt_opt,
                             const string atom_opt,
                             const int mirror_opt, const int het_opt,
                             int &len_aa, int &len_na,  
                             const int o_opt, std::vector<std::string>&resi_vec) {

  // variables used for injection
  // - pos => position data for each chain
  // - seq => sequence data for each chain, sequence name is chain name
  // - rna => whether were dealing with amini acid or nucleotide chains

  // variables that need assignment:
  // - a_vec => position data with 3 dimensions. 1: chain, 2: residue, 3:xyz
  // - seq_vec => sequence data
  // - sec_vec => secondary structure data
  // - mol_vec => one entry per chain, 0 if peptide chain, > 0 if nucleotides
  // - len_vec => length of chains
  // - chainID_list => chain names
  // - len_aa => total number of peptide residues
  // - len_na => total number of nucleotide residues
  // - resi_vec => leave empty... this is only used if byresi_opt is True or
  //               for generating output which is disabled anyways.

  // all other variables are simply ignored but stay here to keep the interface
  // as similar as possible

  int n_chains = pos.size();
  a_vec.resize(n_chains);
  seq_vec.resize(n_chains);
  sec_vec.resize(n_chains);
  mol_vec.resize(n_chains);
  len_vec.resize(n_chains);
  chainID_list.resize(n_chains);
  len_aa = 0;
  len_na = 0;
  for(int i = 0; i < n_chains; ++i) {
    int n_residues = pos[i].size();

    // assign a_vec
    a_vec[i].resize(n_residues, std::vector<double>(3));
    for(int j = 0; j < n_residues; ++j) {
      a_vec[i][j][0] = pos[i][j][0];
      a_vec[i][j][1] = pos[i][j][1];
      a_vec[i][j][2] = pos[i][j][2];
    }

    // assign seq_vec
    const std::string& s = seq[i].GetString();
    seq_vec[i].assign(s.begin(), s.end());

    // assign sec_vec
    double **xa;
    NewArray(&xa, n_residues, 3);
    // yet another conversion needed...
    for(int j = 0; j < n_residues; ++j) {
      xa[j][0] = a_vec[i][j][0];
      xa[j][1] = a_vec[i][j][1];
      xa[j][2] = a_vec[i][j][2];
    }
    char* sec = new char[n_residues + 1];
    if(rna[i]) {
      // make a const cast here... USalign doesn't do anything to that value
      // if it does in the future, thats a recipe for disaster
      make_sec(const_cast<char*>(s.c_str()), xa, n_residues, sec, " C3'");
    } else {
      make_sec(xa, n_residues, sec);
    }
    sec_vec[i].assign(sec, sec + n_residues);
    DeleteArray(&xa, n_residues);
    delete [] sec;

    // assign mol_vec
    if(rna[i]) {
      mol_vec[i] = n_residues;
    } else {
      mol_vec[i] = (-1)*n_residues;
    }

    // assign len_vec
    len_vec[i] = n_residues;

    // assign chainID_list
    // chainID_list[i] = ":1," + seq[i].GetName();
    chainID_list[i] = seq[i].GetName();

    // update length variables
    if(rna[i]) {
      len_na += n_residues;
    } else {
      len_aa += n_residues;
    }
  }
}


TMAlignResult WrappedTMAlign(const geom::Vec3List& pos_one, 
                             const geom::Vec3List& pos_two, 
                             const ost::seq::SequenceHandle& seq1,
                             const ost::seq::SequenceHandle& seq2,
                             bool fast,
                             bool rna) {

  int xlen = pos_one.size();
  int ylen = pos_two.size();  

  if(xlen <= 5 || ylen <= 5) {
    throw ost::Error("Input sequence too short!");
  }

  if(xlen != seq1.GetLength() || ylen != seq2.GetLength()) {
    throw ost::Error("Positions and Sequence must have consistent size "
                     "to run TMAlign!");
  }

  // squeeze input into right format
  char* seqx = new char[xlen+1];
  char* seqy = new char[ylen+1];
  seqx[xlen] = '\0';
  seqy[ylen] = '\0';
  char* secx = new char[xlen];
  char* secy = new char[ylen];

  // use TMalign functionality to generate position arrays
  double** xa;
  double** ya;
  NewArray(&xa, xlen, 3);
  NewArray(&ya, ylen, 3);

  for(int i = 0; i < xlen; ++i) {
    xa[i][0] = pos_one[i][0];
    xa[i][1] = pos_one[i][1];
    xa[i][2] = pos_one[i][2];
    seqx[i] = seq1[i];
  }

  for(int i = 0; i < ylen; ++i) {
    ya[i][0] = pos_two[i][0];
    ya[i][1] = pos_two[i][1];
    ya[i][2] = pos_two[i][2];
    seqy[i] = seq2[i];
  }

  if(rna) {
    make_sec(seqx, xa, xlen, secx, " C3'");
    make_sec(seqy, ya, ylen, secy, " C3'");
  } else {
    make_sec(xa, xlen, secx);
    make_sec(ya, ylen, secy);
  }

  // these variables are chosen such that running TMalign_main is the same as 
  // you would call the executable without any additional parameters
  double Lnorm_ass = -1.0;
  double d0_scale = -1.0;
  bool u_opt = false;
  bool d_opt = false;
  int a_opt = 0; 
  std::vector<String> sequence; 
  bool i_opt = false;
  double TMcut = -1; 

  // following variables are copied from the TMAlign source code
  double t0[3], u0[3][3];
  double TM1, TM2;
  double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
  double d0_0, TM_0;
  double d0A, d0B, d0u, d0a;
  double d0_out=5.0;
  String seqM, seqxA, seqyA;// for output alignment
  double rmsd0 = 0.0;
  int L_ali;                // Aligned length in standard_TMscore
  double Liden=0;
  double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
  int n_ali=0;
  int n_ali8=0;
  int mol_type=static_cast<int>(rna); // Treated as RNA if mol_type > 0

  TMalign_main(xa, ya, seqx, seqy, secx, secy, t0, u0, TM1, TM2, TM3, TM4, TM5,
               d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out, seqM, seqxA, seqyA,
               rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8, xlen, ylen, 
               sequence, Lnorm_ass, d0_scale, i_opt, a_opt, u_opt, d_opt, 
               fast, mol_type, TMcut);

  // cleanup
  DeleteArray(&xa, xlen);
  DeleteArray(&ya, ylen);
  delete [] seqx;
  delete [] seqy;
  delete [] secx;
  delete [] secy;

  // collect results and return
  TMAlignResult res;
  res.tm_score = TM1;
  res.rmsd = rmsd0;
  res.aligned_length = n_ali8;
  res.transform = geom::Mat4(u0[0][0], u0[0][1], u0[0][2], t0[0],
                             u0[1][0], u0[1][1], u0[1][2], t0[1],
                             u0[2][0], u0[2][1], u0[2][2], t0[2],
                             0.0, 0.0, 0.0, 1.0);
  res.alignment = ost::seq::CreateAlignment();
  ost::seq::SequenceHandle aligned_seq1 = 
  ost::seq::CreateSequence(seq1.GetName(), seqxA);
  ost::seq::SequenceHandle aligned_seq2 = 
  ost::seq::CreateSequence(seq2.GetName(), seqyA);
  res.alignment.AddSequence(aligned_seq1);
  res.alignment.AddSequence(aligned_seq2);
  return res;
}

MMAlignResult WrappedMMAlign(const std::vector<geom::Vec3List>& pos_one,
                             const std::vector<geom::Vec3List>& pos_two,
                             const ost::seq::SequenceList& seq1,
                             const ost::seq::SequenceList& seq2,
                             const std::vector<bool>& rna1,
                             const std::vector<bool>& rna2,
                             bool fast) {

  // input checks
  if(pos_one.empty() || pos_two.empty()) {
    throw ost::Error("Cannot compute MMAlign on empty chains!");
  }

  if(static_cast<int>(pos_one.size()) != seq1.GetCount() ||
     pos_one.size() != rna1.size()) {
    throw ost::Error("Inconsistent input sizes in WrappedMMAlign");
  }

  if(static_cast<int>(pos_two.size()) != seq2.GetCount() ||
     pos_two.size() != rna2.size()) {
    throw ost::Error("Inconsistent input sizes in WrappedMMAlign");
  }

  if(pos_one.size() == 1 && pos_two.size() == 1) {
    // just run TMAlign...
    if(rna1[0] != rna2[0]) {
      throw ost::Error("Error in WrappedMMAlign: If both complexes only have "
                       "one chain, they must either be both peptide or both "
                       "RNA.");
    }
    TMAlignResult tm_result = WrappedTMAlign(pos_one[0], pos_two[0], seq1[0],
                                             seq2[0], fast, rna1[0]);
    ost::seq::AlignmentList alns;
    alns.push_back(tm_result.alignment);
    std::vector<String> ent1_mapped_chains;
    std::vector<String> ent2_mapped_chains;
    ent1_mapped_chains.push_back(seq1[0].GetName());
    ent2_mapped_chains.push_back(seq2[0].GetName());
    return MMAlignResult(tm_result.rmsd, tm_result.tm_score,
                         tm_result.transform, tm_result.aligned_length,
                         alns, ent1_mapped_chains,
                         ent2_mapped_chains);
  }

  // the following is a copy of the variable definition section in USalign.cpp
  // main function to get default param
  // changes to the defaults there are explicitely marked
  // most with: silence compiler warning due to unused variables
  std::string xname       = "";
  std::string yname       = "";
  std::string fname_super = ""; // file name for superposed structure
  std::string fname_lign  = ""; // file name for user alignment
  std::string fname_matrix= ""; // file name for output matrix
  vector<std::string> sequence; // get value from alignment file
  //double Lnorm_ass, d0_scale; // silence compiler warning
  double d0_scale = 0.0; // only d0 scale required, directly initialize with
                         // default value to silence compiler warning

  // silence compiler warning
  //bool h_opt = false; // print full help message
  // silence compiler warning
  //bool v_opt = false; // print version
  bool m_opt = false; // flag for -m, output rotation matrix
  int  i_opt = 0;     // 1 for -i, 3 for -I
  int  o_opt = 0;     // 1 for -o, 2 for -rasmol
  int  a_opt = 0;     // flag for -a, do not normalized by average length
  // silence compiler warning
  //bool u_opt = false; // flag for -u, normalized by user specified length
  bool d_opt = false; // flag for -d, user specified d0

  bool   full_opt  = false;// do not show chain level alignment
  double TMcut     =-1;
  // silence compiler warning
  //bool   se_opt    =false;
  int    infmt1_opt=-1;    // PDB or PDBx/mmCIF format for chain_1
  int    infmt2_opt=-1;    // PDB or PDBx/mmCIF format for chain_2
  int    ter_opt   =2;     // END, or different chainID
  int    split_opt =2;     // split each chains
  int    outfmt_opt=0;     // set -outfmt to full output
  bool   fast_opt  =false; // flags for -fast, fTM-align algorithm
  // silence compiler warning
  //int    cp_opt    =0;     // do not check circular permutation
  // silence compiler warning
  //int    closeK_opt=-1;    // number of atoms for SOI initial alignment.
                             // 5 and 0 for -mm 5 and 6
  // silence compiler warning
  // int    hinge_opt =9;     // maximum number of hinge allowed for flexible
  int    mirror_opt=0;     // do not align mirror
  int    het_opt=0;        // do not read HETATM residues
  // silence compiler warning
  // int    mm_opt=0;         // do not perform MM-align
  std::string atom_opt  ="auto";// use C alpha atom for protein and C3' for RNA
  std::string mol_opt   ="auto";// auto-detect the molecule type as protein/RNA
  std::string suffix_opt="";    // set -suffix to empty
  std::string dir_opt   ="";    // set -dir to empty
  std::string dir1_opt  ="";    // set -dir1 to empty
  std::string dir2_opt  ="";    // set -dir2 to empty
  int    byresi_opt=0;     // set -byresi to 0
  vector<std::string> chain1_list; // only when -dir1 is set
  vector<std::string> chain2_list; // only when -dir2 is set


  // The following is pretty much a copy of the MMalign function in USalign.cpp
  // with adaptions to inject our own data structures

  /* declare previously global variables */
  vector<vector<vector<double> > > xa_vec; // structure of complex1
  vector<vector<vector<double> > > ya_vec; // structure of complex2
  vector<vector<char> >seqx_vec; // sequence of complex1
  vector<vector<char> >seqy_vec; // sequence of complex2
  vector<vector<char> >secx_vec; // secondary structure of complex1
  vector<vector<char> >secy_vec; // secondary structure of complex2
  vector<int> mol_vec1;          // molecule type of complex1, RNA if >0
  vector<int> mol_vec2;          // molecule type of complex2, RNA if >0
  vector<string> chainID_list1;  // list of chainID1
  vector<string> chainID_list2;  // list of chainID2
  vector<int> xlen_vec;          // length of complex1
  vector<int> ylen_vec;          // length of complex2
  int    i,j;                    // chain index
  int    xlen, ylen;             // chain length
  double **xa, **ya;             // structure of single chain
  char   *seqx, *seqy;           // for the protein sequence 
  char   *secx, *secy;           // for the secondary structure 
  int    xlen_aa,ylen_aa;        // total length of protein
  int    xlen_na,ylen_na;        // total length of RNA/DNA
  vector<string> resi_vec1;  // residue index for chain1
  vector<string> resi_vec2;  // residue index for chain2


  // COMMENT OUT MMALIGN STRUCTURE PARSING WHICH IS FILE BASED
  ////////////////////////////////////////////////////////////

  ///* parse complex */
  //parse_chain_list(chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
  //    xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
  //    atom_opt, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt, resi_vec1);
  //if (xa_vec.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 1");
  //parse_chain_list(chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
  //    ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
  //    atom_opt, 0, het_opt, ylen_aa, ylen_na, o_opt, resi_vec2);
  //if (ya_vec.size()==0) PrintErrorAndQuit("ERROR! 0 chain in complex 2");

  // INJECT OWN DATA
  //////////////////
  parse_chain_list_hacked(pos_one, seq1, rna1,
      chain1_list, xa_vec, seqx_vec, secx_vec, mol_vec1,
      xlen_vec, chainID_list1, ter_opt, split_opt, mol_opt, infmt1_opt,
      atom_opt, mirror_opt, het_opt, xlen_aa, xlen_na, o_opt, resi_vec1);
  parse_chain_list_hacked(pos_two, seq2, rna2,
      chain2_list, ya_vec, seqy_vec, secy_vec, mol_vec2,
      ylen_vec, chainID_list2, ter_opt, split_opt, mol_opt, infmt2_opt,
      atom_opt, 0, het_opt, ylen_aa, ylen_na, o_opt, resi_vec2);

  int len_aa=getmin(xlen_aa,ylen_aa);
  int len_na=getmin(xlen_na,ylen_na);
  if (a_opt)
  {
      len_aa=(xlen_aa+ylen_aa)/2;
      len_na=(xlen_na+ylen_na)/2;
  }
  if (byresi_opt) i_opt=3; 

  // this is already handled above
  ////////////////////////////////
  
  ///* perform monomer alignment if there is only one chain */
  //if (xa_vec.size()==1 && ya_vec.size()==1)
  //{
  //    xlen = xlen_vec[0];
  //    ylen = ylen_vec[0];
  //    seqx = new char[xlen+1];
  //    seqy = new char[ylen+1];
  //    secx = new char[xlen+1];
  //    secy = new char[ylen+1];
  //    NewArray(&xa, xlen, 3);
  //    NewArray(&ya, ylen, 3);
  //    copy_chain_data(xa_vec[0],seqx_vec[0],secx_vec[0], xlen,xa,seqx,secx);
  //    copy_chain_data(ya_vec[0],seqy_vec[0],secy_vec[0], ylen,ya,seqy,secy);
  //    
  //    /* declare variable specific to this pair of TMalign */
  //    double t0[3], u0[3][3];
  //    double TM1, TM2;
  //    double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
  //    double d0_0, TM_0;
  //    double d0A, d0B, d0u, d0a;
  //    double d0_out=5.0;
  //    string seqM, seqxA, seqyA;// for output alignment
  //    double rmsd0 = 0.0;
  //    int L_ali;                // Aligned length in standard_TMscore
  //    double Liden=0;
  //    double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
  //    int n_ali=0;
  //    int n_ali8=0;
  //    
  //    if (byresi_opt) extract_aln_from_resi(sequence,
  //        seqx,seqy,resi_vec1,resi_vec2,byresi_opt);
  
  //    /* entry function for structure alignment */
  //    TMalign_main(xa, ya, seqx, seqy, secx, secy,
  //        t0, u0, TM1, TM2, TM3, TM4, TM5,
  //        d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
  //        seqM, seqxA, seqyA,
  //        rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
  //        xlen, ylen, sequence, 0, d0_scale,
  //        i_opt, a_opt, false, d_opt, fast_opt,
  //        mol_vec1[0]+mol_vec2[0],TMcut);
  
  //    /* print result */
  //    output_results(
  //        xname.substr(dir1_opt.size()),
  //        yname.substr(dir2_opt.size()),
  //        chainID_list1[0], chainID_list2[0],
  //        xlen, ylen, t0, u0, TM1, TM2, TM3, TM4, TM5, rmsd0, d0_out,
  //        seqM.c_str(), seqxA.c_str(), seqyA.c_str(), Liden,
  //        n_ali8, L_ali, TM_ali, rmsd_ali, TM_0, d0_0, d0A, d0B,
  //        0, d0_scale, d0a, d0u, (m_opt?fname_matrix:"").c_str(),
  //        outfmt_opt, ter_opt, true, split_opt, o_opt, fname_super,
  //        0, a_opt, false, d_opt, mirror_opt, resi_vec1, resi_vec2);
  
  //    /* clean up */
  //    seqM.clear();
  //    seqxA.clear();
  //    seqyA.clear();
  //    delete[]seqx;
  //    delete[]seqy;
  //    delete[]secx;
  //    delete[]secy;
  //    DeleteArray(&xa,xlen);
  //    DeleteArray(&ya,ylen);
  
  //    vector<vector<vector<double> > >().swap(xa_vec); // structure of complex1
  //    vector<vector<vector<double> > >().swap(ya_vec); // structure of complex2
  //    vector<vector<char> >().swap(seqx_vec); // sequence of complex1
  //    vector<vector<char> >().swap(seqy_vec); // sequence of complex2
  //    vector<vector<char> >().swap(secx_vec); // secondary structure of complex1
  //    vector<vector<char> >().swap(secy_vec); // secondary structure of complex2
  //    mol_vec1.clear();       // molecule type of complex1, RNA if >0
  //    mol_vec2.clear();       // molecule type of complex2, RNA if >0
  //    chainID_list1.clear();  // list of chainID1
  //    chainID_list2.clear();  // list of chainID2
  //    xlen_vec.clear();       // length of complex1
  //    ylen_vec.clear();       // length of complex2
  //    return 0;
  //}

    /* declare TM-score tables */
    int chain1_num=xa_vec.size();
    int chain2_num=ya_vec.size();
    vector<string> tmp_str_vec(chain2_num,"");
    double **TMave_mat;
    double **ut_mat; // rotation matrices for all-against-all alignment
    int ui,uj,ut_idx;
    NewArray(&TMave_mat,chain1_num,chain2_num);
    NewArray(&ut_mat,chain1_num*chain2_num,4*3);
    vector<vector<string> >seqxA_mat(chain1_num,tmp_str_vec);
    vector<vector<string> > seqM_mat(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_mat(chain1_num,tmp_str_vec);

    double maxTMmono=-1;
    int maxTMmono_i,maxTMmono_j;

    // assign default values to silence compiler warnings
    maxTMmono_i = -1;
    maxTMmono_j = -1;

    /* get all-against-all alignment */
    if (len_aa+len_na>500) fast_opt=true;
    for (i=0;i<chain1_num;i++)
    {
        xlen=xlen_vec[i];
        if (xlen<3)
        {
            for (j=0;j<chain2_num;j++) TMave_mat[i][j]=-1;
            continue;
        }
        seqx = new char[xlen+1];
        secx = new char[xlen+1];
        NewArray(&xa, xlen, 3);
        copy_chain_data(xa_vec[i],seqx_vec[i],secx_vec[i],
            xlen,xa,seqx,secx);

        for (j=0;j<chain2_num;j++)
        {
            ut_idx=i*chain2_num+j;
            for (ui=0;ui<4;ui++)
                for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=0;
            ut_mat[ut_idx][0]=1;
            ut_mat[ut_idx][4]=1;
            ut_mat[ut_idx][8]=1;

            if (mol_vec1[i]*mol_vec2[j]<0) //no protein-RNA alignment
            {
                TMave_mat[i][j]=-1;
                continue;
            }

            ylen=ylen_vec[j];
            if (ylen<3)
            {
                TMave_mat[i][j]=-1;
                continue;
            }
            seqy = new char[ylen+1];
            secy = new char[ylen+1];
            NewArray(&ya, ylen, 3);
            copy_chain_data(ya_vec[j],seqy_vec[j],secy_vec[j],
                ylen,ya,seqy,secy);

            /* declare variable specific to this pair of TMalign */
            double t0[3], u0[3][3];
            double TM1, TM2;
            double TM3, TM4, TM5;     // for a_opt, u_opt, d_opt
            double d0_0, TM_0;
            double d0A, d0B, d0u, d0a;
            double d0_out=5.0;
            string seqM, seqxA, seqyA;// for output alignment
            double rmsd0 = 0.0;
            int L_ali;                // Aligned length in standard_TMscore
            double Liden=0;
            double TM_ali, rmsd_ali;  // TMscore and rmsd in standard_TMscore
            int n_ali=0;
            int n_ali8=0;

            int Lnorm_tmp=len_aa;
            if (mol_vec1[i]+mol_vec2[j]>0) Lnorm_tmp=len_na;
            
            if (byresi_opt)
            {
                int total_aln=extract_aln_from_resi(sequence,
                    seqx,seqy,resi_vec1,resi_vec2,xlen_vec,ylen_vec, i, j);
                seqxA_mat[i][j]=sequence[0];
                seqyA_mat[i][j]=sequence[1];
                if (total_aln>xlen+ylen-3)
                {
                    for (ui=0;ui<3;ui++) for (uj=0;uj<3;uj++) 
                        ut_mat[ut_idx][ui*3+uj]=(ui==uj)?1:0;
                    for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=0;
                    TMave_mat[i][j]=0;
                    seqM.clear();
                    seqxA.clear();
                    seqyA.clear();

                    delete[]seqy;
                    delete[]secy;
                    DeleteArray(&ya,ylen);
                    continue;
                }
            }

            /* entry function for structure alignment */
            TMalign_main(xa, ya, seqx, seqy, secx, secy,
                t0, u0, TM1, TM2, TM3, TM4, TM5,
                d0_0, TM_0, d0A, d0B, d0u, d0a, d0_out,
                seqM, seqxA, seqyA,
                rmsd0, L_ali, Liden, TM_ali, rmsd_ali, n_ali, n_ali8,
                xlen, ylen, sequence, Lnorm_tmp, d0_scale,
                i_opt, false, true, false, fast_opt,
                mol_vec1[i]+mol_vec2[j],TMcut);

            /* store result */
            for (ui=0;ui<3;ui++)
                for (uj=0;uj<3;uj++) ut_mat[ut_idx][ui*3+uj]=u0[ui][uj];
            for (uj=0;uj<3;uj++) ut_mat[ut_idx][9+uj]=t0[uj];
            seqxA_mat[i][j]=seqxA;
            seqyA_mat[i][j]=seqyA;
            TMave_mat[i][j]=TM4*Lnorm_tmp;
            if (TMave_mat[i][j]>maxTMmono)
            {
                maxTMmono=TMave_mat[i][j];
                maxTMmono_i=i;
                maxTMmono_j=j;
            }

            /* clean up */
            seqM.clear();
            seqxA.clear();
            seqyA.clear();

            delete[]seqy;
            delete[]secy;
            DeleteArray(&ya,ylen);
        }

        delete[]seqx;
        delete[]secx;
        DeleteArray(&xa,xlen);
    }

    /* calculate initial chain-chain assignment */
    int *assign1_list; // value is index of assigned chain2
    int *assign2_list; // value is index of assigned chain1
    assign1_list=new int[chain1_num];
    assign2_list=new int[chain2_num];
    double total_score=enhanced_greedy_search(TMave_mat, assign1_list,
        assign2_list, chain1_num, chain2_num);
    if (total_score<=0) PrintErrorAndQuit("ERROR! No assignable chain");

    /* refine alignment for large oligomers */
    int aln_chain_num=count_assign_pair(assign1_list,chain1_num);
    bool is_oligomer=(aln_chain_num>=3);
    if (aln_chain_num==2) // dimer alignment
    {
        int na_chain_num1,na_chain_num2,aa_chain_num1,aa_chain_num2;
        count_na_aa_chain_num(na_chain_num1,aa_chain_num1,mol_vec1);
        count_na_aa_chain_num(na_chain_num2,aa_chain_num2,mol_vec2);

        /* align protein-RNA hybrid dimer to another hybrid dimer */
        if (na_chain_num1==1 && na_chain_num2==1 && 
            aa_chain_num1==1 && aa_chain_num2==1) is_oligomer=false;
        /* align pure protein dimer or pure RNA dimer */
        else if ((getmin(na_chain_num1,na_chain_num2)==0 && 
                    aa_chain_num1==2 && aa_chain_num2==2) ||
                 (getmin(aa_chain_num1,aa_chain_num2)==0 && 
                    na_chain_num1==2 && na_chain_num2==2))
        {
            adjust_dimer_assignment(xa_vec,ya_vec,xlen_vec,ylen_vec,mol_vec1,
                mol_vec2,assign1_list,assign2_list,seqxA_mat,seqyA_mat);
            is_oligomer=false; // cannot refiner further
        }
        else is_oligomer=true; /* align oligomers to dimer */
    }

    if (aln_chain_num>=3 || is_oligomer) // oligomer alignment
    {
        /* extract centroid coordinates */
        double **xcentroids;
        double **ycentroids;
        NewArray(&xcentroids, chain1_num, 3);
        NewArray(&ycentroids, chain2_num, 3);
        double d0MM=getmin(
            calculate_centroids(xa_vec, chain1_num, xcentroids),
            calculate_centroids(ya_vec, chain2_num, ycentroids));

        /* refine enhanced greedy search with centroid superposition */
        //double het_deg=check_heterooligomer(TMave_mat, chain1_num, chain2_num);
        homo_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na, ut_mat);
        hetero_refined_greedy_search(TMave_mat, assign1_list,
            assign2_list, chain1_num, chain2_num, xcentroids,
            ycentroids, d0MM, len_aa+len_na);
        
        /* clean up */
        DeleteArray(&xcentroids, chain1_num);
        DeleteArray(&ycentroids, chain2_num);
    }

    /* store initial assignment */
    int init_pair_num=count_assign_pair(assign1_list,chain1_num);
    int *assign1_init, *assign2_init;
    assign1_init=new int[chain1_num];
    assign2_init=new int[chain2_num];
    double **TMave_init;
    NewArray(&TMave_init,chain1_num,chain2_num);
    vector<vector<string> >seqxA_init(chain1_num,tmp_str_vec);
    vector<vector<string> >seqyA_init(chain1_num,tmp_str_vec);
    vector<string> sequence_init;
    copy_chain_assign_data(chain1_num, chain2_num, sequence_init,
        seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init);

    /* perform iterative alignment */
    double max_total_score=0; // ignore old total_score because previous
                              // score was from monomeric chain superpositions
    int max_iter=5-(int)((len_aa+len_na)/200);
    if (max_iter<2) max_iter=2;
    if (byresi_opt==0) MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec,
        seqx_vec, seqy_vec, secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec,
        ylen_vec, xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num,
        chain2_num, TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list,
        sequence, d0_scale, fast_opt);

    /* sometime MMalign_iter is even worse than monomer alignment */
    if (byresi_opt==0 && max_total_score<maxTMmono)
    {
        copy_chain_assign_data(chain1_num, chain2_num, sequence,
            seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
            seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat);
        for (i=0;i<chain1_num;i++)
        {
            if (i!=maxTMmono_i) assign1_list[i]=-1;
            else assign1_list[i]=maxTMmono_j;
        }
        for (j=0;j<chain2_num;j++)
        {
            if (j!=maxTMmono_j) assign2_list[j]=-1;
            else assign2_list[j]=maxTMmono_i;
        }
        sequence[0]=seqxA_mat[maxTMmono_i][maxTMmono_j];
        sequence[1]=seqyA_mat[maxTMmono_i][maxTMmono_j];
        max_total_score=maxTMmono;
        MMalign_iter(max_total_score, max_iter, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
            TMave_mat, seqxA_mat, seqyA_mat, assign1_list, assign2_list, sequence,
            d0_scale, fast_opt);
    }

    /* perform cross chain alignment
     * in some cases, this leads to dramatic improvement, esp for homodimer */
    int iter_pair_num=count_assign_pair(assign1_list,chain1_num);
    if (iter_pair_num>=init_pair_num) copy_chain_assign_data(
        chain1_num, chain2_num, sequence_init,
        seqxA_mat, seqyA_mat, assign1_list, assign2_list, TMave_mat,
        seqxA_init, seqyA_init, assign1_init,  assign2_init,  TMave_init);
    double max_total_score_cross=max_total_score;
    if (byresi_opt==0 && len_aa+len_na<10000)
    {
        MMalign_dimer(max_total_score_cross, xa_vec, ya_vec, seqx_vec, seqy_vec,
            secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
            xa, ya, seqx, seqy, secx, secy, len_aa, len_na, chain1_num, chain2_num,
            TMave_init, seqxA_init, seqyA_init, assign1_init, assign2_init,
            sequence_init, d0_scale, fast_opt);
        if (max_total_score_cross>max_total_score) 
        {
            max_total_score=max_total_score_cross;
            copy_chain_assign_data(chain1_num, chain2_num, sequence,
                seqxA_init, seqyA_init, assign1_init, assign2_init, TMave_init,
                seqxA_mat,  seqyA_mat,  assign1_list, assign2_list, TMave_mat);
        }
    } 

    /* final alignment */

    // commented out by Gabriel => avoid include that defines print_version 
    //if (outfmt_opt==0) print_version(); 

    // Call hacked version of MMalign_final that returns MMAlignResult object

    MMAlignResult res = MMalign_final_hacked(xname.substr(dir1_opt.size()), yname.substr(dir2_opt.size()),
        chainID_list1, chainID_list2,
        fname_super, fname_lign, fname_matrix,
        xa_vec, ya_vec, seqx_vec, seqy_vec,
        secx_vec, secy_vec, mol_vec1, mol_vec2, xlen_vec, ylen_vec,
        xa, ya, seqx, seqy, secx, secy, len_aa, len_na,
        chain1_num, chain2_num, TMave_mat,
        seqxA_mat, seqM_mat, seqyA_mat, assign1_list, assign2_list, sequence,
        d0_scale, m_opt, o_opt, outfmt_opt, ter_opt, split_opt,
        a_opt, d_opt, fast_opt, full_opt, mirror_opt, resi_vec1, resi_vec2);

    /* clean up everything */
    delete [] assign1_list;
    delete [] assign2_list;
    DeleteArray(&TMave_mat,chain1_num);
    DeleteArray(&ut_mat,   chain1_num*chain2_num);
    vector<vector<string> >().swap(seqxA_mat);
    vector<vector<string> >().swap(seqM_mat);
    vector<vector<string> >().swap(seqyA_mat);
    vector<string>().swap(tmp_str_vec);

    delete [] assign1_init;
    delete [] assign2_init;
    DeleteArray(&TMave_init,chain1_num);
    vector<vector<string> >().swap(seqxA_init);
    vector<vector<string> >().swap(seqyA_init);

    vector<vector<vector<double> > >().swap(xa_vec); // structure of complex1
    vector<vector<vector<double> > >().swap(ya_vec); // structure of complex2
    vector<vector<char> >().swap(seqx_vec); // sequence of complex1
    vector<vector<char> >().swap(seqy_vec); // sequence of complex2
    vector<vector<char> >().swap(secx_vec); // secondary structure of complex1
    vector<vector<char> >().swap(secy_vec); // secondary structure of complex2
    mol_vec1.clear();       // molecule type of complex1, RNA if >0
    mol_vec2.clear();       // molecule type of complex2, RNA if >0
    vector<string>().swap(chainID_list1);  // list of chainID1
    vector<string>().swap(chainID_list2);  // list of chainID2
    xlen_vec.clear();       // length of complex1
    ylen_vec.clear();       // length of complex2
    vector<string> ().swap(resi_vec1);  // residue index for chain1
    vector<string> ().swap(resi_vec2);  // residue index for chain2


  return res;
}

void ExtractChainInfo(const ost::mol::ChainView& chain, geom::Vec3List& pos,
                      ost::seq::SequenceHandle& s, bool& rna_mode) {

  pos.clear();
  std::vector<char> olcs;
  rna_mode = false;
  ost::mol::ResidueViewList res_list = chain.GetResidueList();

  for(auto it = res_list.begin(); it != res_list.end(); ++it) {
    char olc = it->GetOneLetterCode();
    if(olc == '?') {
      continue;
    }
    if(it->IsPeptideLinking()) {
      ost::mol::AtomView ca = it->FindAtom("CA");
      if(!ca.IsValid()) {
        continue;
      }
      if(rna_mode) {
        std::stringstream ss;
        ss << "Error in WrappedTMAlign: Chains cannot have peptide and RNA ";
        ss << "residues in same chain. Problematic chain: "<<chain.GetName();
        throw ost::Error(ss.str());
      }
      olcs.push_back(olc);
      pos.push_back(ca.GetPos());
    }
    else if(it->IsNucleotideLinking()) {
      ost::mol::AtomView c3 = it->FindAtom("C3'");
      if(!c3.IsValid()) {
        continue;
      }
      if(rna_mode==false && !pos.empty()) {
        std::stringstream ss;
        ss << "Error in WrappedTMAlign: Chains cannot have peptide and RNA ";
        ss << "residues in same chain. Problematic chain: "<<chain.GetName();
        throw ost::Error(ss.str());
      }
      rna_mode = true;
      // for some reason, USalign wants nucleotides to be lower case
      olcs.push_back(tolower(olc));
      pos.push_back(c3.GetPos());
    }
  }
  String str_s = String(olcs.begin(), olcs.end());
  s = ost::seq::CreateSequence(chain.GetName(), str_s);
}


TMAlignResult WrappedTMAlign(const ost::mol::ChainView& chain1,
                             const ost::mol::ChainView& chain2,
                             bool fast) {

  geom::Vec3List pos1;
  ost::seq::SequenceHandle s1;
  bool rna_mode1;
  ExtractChainInfo(chain1, pos1, s1, rna_mode1);

  geom::Vec3List pos2;
  ost::seq::SequenceHandle s2;
  bool rna_mode2;
  ExtractChainInfo(chain2, pos2, s2, rna_mode2);

  if(rna_mode1 != rna_mode2) {
    throw ost::Error("Error in WrappedTMAlign: Cannot compare peptide with "
                     "RNA chains");
  }

  return WrappedTMAlign(pos1, pos2, s1, s2, fast, rna_mode1);
}

MMAlignResult WrappedMMAlign(const ost::mol::EntityView& ent1,
                             const ost::mol::EntityView& ent2,
                             bool fast) {
  ost::mol::ChainViewList chains1 = ent1.GetChainList();
  int n1 = chains1.size();
  std::vector<geom::Vec3List> pos1(n1);
  ost::seq::SequenceList s1 = ost::seq::CreateSequenceList();
  std::vector<bool> rna1(n1);
  for(int i = 0; i < n1; ++i) {
    bool rna;
    ost::seq::SequenceHandle s;
    ExtractChainInfo(chains1[i], pos1[i], s, rna);
    rna1[i] = rna;
    s1.AddSequence(s);
  }

  ost::mol::ChainViewList chains2 = ent2.GetChainList();
  int n2 = chains2.size();
  std::vector<geom::Vec3List> pos2(n2);
  ost::seq::SequenceList s2 = ost::seq::CreateSequenceList();
  std::vector<bool> rna2(n2);
  for(int i = 0; i < n2; ++i) {
    bool rna;
    ost::seq::SequenceHandle s;
    ExtractChainInfo(chains2[i], pos2[i], s, rna);
    rna2[i] = rna;
    s2.AddSequence(s);
  }

  return WrappedMMAlign(pos1, pos2, s1, s2, rna1, rna2, fast);
}

}} //ns
