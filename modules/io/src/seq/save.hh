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
#ifndef OST_IO_SEQ_SAVE_HH
#define OST_IO_SEQ_SAVE_HH

/*
  Author: Marco Biasini
 */
#include <ost/io/module_config.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/sequence_list.hh>
namespace ost { namespace io {

void DLLEXPORT_OST_IO SaveSequenceList(const seq::ConstSequenceList& seq_list,
                                    const String& filename,
                                    const String& format="auto");
                                      
String DLLEXPORT_OST_IO 
SequenceListToString(const seq::ConstSequenceList& seq_list,
                     const String& format);
                                          
void DLLEXPORT_OST_IO 
SaveSequence(const seq::ConstSequenceHandle& sequence, const String& filename,
             const String& format="auto");
                                      
             
String DLLEXPORT_OST_IO
SequenceToString(const seq::ConstSequenceHandle& sequence,
                 const String& format="auto");

void DLLEXPORT_OST_IO SaveAlignment(const seq::AlignmentHandle& msa,
                                    const String& filename,
                                    const String& format="auto");

String DLLEXPORT_OST_IO AlignmentToString(const seq::AlignmentHandle& msa,
                                          const String& format);

}}

#endif
