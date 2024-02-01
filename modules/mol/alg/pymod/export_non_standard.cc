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
#include <boost/python.hpp>
#include <ost/mol/mol.hh>
#include <ost/mol/alg/nonstandard.hh>

using namespace boost::python;

using namespace ost::mol::alg;
using namespace ost::mol;
                   
object copy_conserved_handle(ResidueHandle src_res, ResidueHandle dst_res,
                          XCSEditor edi) {
  bool has_cbeta = false;
  bool ret = CopyConserved(src_res, dst_res, edi, has_cbeta);
  return make_tuple(ret, has_cbeta);
}

object copy_non_conserved_handle(ResidueHandle src_res, ResidueHandle dst_res,
                          XCSEditor edi) {
  bool has_cbeta = false;
  bool ret = CopyNonConserved(src_res, dst_res, edi, has_cbeta);
  return make_tuple(ret, has_cbeta);
}

bool (*copy_residue_handle)(ost::mol::ResidueHandle,
                            ost::mol::ResidueHandle,
                            ost::mol::XCSEditor&)=&CopyResidue;


void export_NonStandard()
{
  def("CopyNonConserved",&copy_non_conserved_handle);
  def("CopyConserved", copy_conserved_handle);
  def("CopyResidue", copy_residue_handle);
 }

