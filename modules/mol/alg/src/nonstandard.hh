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

#ifndef OST_MOL_ALG_NONSTANDARD_HH
#define OST_MOL_ALG_NONSTANDARD_HH
/*
  Author: Marco Biasini, Juergen Haas
 */
#include "module_config.hh"
#include <ost/conop/compound_lib.hh>



namespace ost { namespace mol { namespace alg {


/// \brief copies all atom of src_res to dst_res, gets compound lib from builder

                             
bool DLLEXPORT_OST_MOL_ALG CopyResidue(ost::mol::ResidueHandle src_res,
                                     ost::mol::ResidueHandle dst_res,
                                     ost::mol::XCSEditor& edi);

/// \brief copies all atom of src_res to dst_res, requires a compound lib 
bool DLLEXPORT_OST_MOL_ALG CopyResidue(ost::mol::ResidueHandle src_res,
                                     ost::mol::ResidueHandle dst_res,
                                     ost::mol::XCSEditor& edi,
                                     ost::conop::CompoundLibPtr lib);

                             
/// \brief copies all atom of src_res to dst_res
/// \param src_res
/// \param dst_res
/// \param edi
/// \param has_cbeta will be set to true if the src_res has a cbeta and the
///      dst_residue is not a glycine                             
bool DLLEXPORT_OST_MOL_ALG CopyIdentical(ost::mol::ResidueHandle src_res,
                                       ost::mol::ResidueHandle dst_res,
                                       ost::mol::XCSEditor& edi,
                                       bool& has_cbeta);

/// \brief copies atoms of src_res to dst_res, gets compound lib from builder
///
/// src_res and dst_res are thought to be conserved, e.g. the parent standard 
/// amino acid of both residues is the same. This includes cases where e.g. the 
/// src_rs is and MSE and the dst_res is a MET. This function automatically 
/// tries to do the right thing an keep as many atoms as possible from src_res


                             
bool DLLEXPORT_OST_MOL_ALG CopyConserved(ost::mol::ResidueHandle src_res, 
                                       ost::mol::ResidueHandle dst_res,
                                       ost::mol::XCSEditor& edi,
                                       bool& has_cbeta,
                                       ost::conop::CompoundLibPtr lib);

/// \brief copies atoms of src_res to dst_res, requires compound lib
///
/// src_res and dst_res are thought to be conserved, e.g. the parent standard 
/// amino acid of both residues is the same. This includes cases where e.g. the 
/// src_rs is and MSE and the dst_res is a MET. This function automatically 
/// tries to do the right thing an keep as many atoms as possible from src_res


                             
bool DLLEXPORT_OST_MOL_ALG CopyConserved(ost::mol::ResidueHandle src_res, 
                                       ost::mol::ResidueHandle dst_res,
                                       ost::mol::XCSEditor& edi,
                                       bool& has_cbeta);


/// \brief construct dst_res in case src_res and dst_res are not conserved.
/// 
/// This essentially copies the backbone of src_res to dst_res. The CB atom is 
/// only copied if dst_res is not equal to glycine.


bool DLLEXPORT_OST_MOL_ALG CopyNonConserved(ost::mol::ResidueHandle src_res, 
                                          ost::mol::ResidueHandle dst_res,
                                          ost::mol::XCSEditor& edi,
                                          bool& has_cbeta);

/// \brief construct dst_res from src_res when src_res is an MSE

bool DLLEXPORT_OST_MOL_ALG CopyMSE(ost::mol::ResidueHandle src_res, 
                                 ost::mol::ResidueHandle dst_res,
                                 ost::mol::XCSEditor& edi,
                                 bool& has_cbeta);
                       
/// \brief construct a dst_res with only atoms matching the standard aminoacid
/// from src_res when src_res is an is modified
                            
bool DLLEXPORT_OST_MOL_ALG CopyModified(ost::mol::ResidueHandle src_res, 
                                      ost::mol::ResidueHandle dst_res,
                                      ost::mol::XCSEditor& edi,
                                      bool& has_cbeta,
                                      ost::conop::CompoundLibPtr lib);



}}}

#endif
