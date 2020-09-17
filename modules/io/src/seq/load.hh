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
#ifndef OST_IO_LOAD_ALIGNMENT_HH
#define OST_IO_LOAD_ALIGNMENT_HH

/*
  Author: Marco Biasini
 */

#include <ost/io/module_config.hh>
#include <ost/seq/sequence_list.hh>
#include <ost/seq/alignment_handle.hh>
#include <ost/seq/profile_handle.hh>

namespace ost { namespace io {

seq::AlignmentHandle DLLEXPORT_OST_IO 
LoadAlignment(const String& file_name, const String& format="auto");

seq::AlignmentHandle DLLEXPORT_OST_IO 
AlignmentFromStream(std::istream&  stream, const String& format);

seq::AlignmentHandle DLLEXPORT_OST_IO 
AlignmentFromString(const String& data, const String& format);


seq::SequenceList DLLEXPORT_OST_IO 
LoadSequenceList(const String& file_name, const String& format="auto");

seq::SequenceList DLLEXPORT_OST_IO 
SequenceListFromStream(std::istream&  stream, const String& format);

seq::SequenceList DLLEXPORT_OST_IO 
SequenceListFromString(const String& data, const String& format);

seq::SequenceHandle DLLEXPORT_OST_IO 
LoadSequence(const String& file_name, const String& format="auto");

seq::SequenceHandle DLLEXPORT_OST_IO 
SequenceFromStream(std::istream&  stream, const String& format);

seq::SequenceHandle DLLEXPORT_OST_IO 
SequenceFromString(const String& data, const String& format);

seq::ProfileHandlePtr DLLEXPORT_OST_IO 
LoadSequenceProfile(const String& file_name, const String& format="auto");

seq::ProfileHandlePtr DLLEXPORT_OST_IO
SequenceProfileFromString(const String& data, const String& format);

}}

#endif
