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
#include "load.hh"
#include <ost/io/io_manager.hh>

namespace ost { namespace io {


seq::SequenceList LoadSequenceList(const String& file_name, 
                                   const String& format)
{
  seq::SequenceList a=seq::CreateSequenceList();
  IOManager& m=IOManager::Instance();                               
  SequenceIOHandlerPtr ali_io=m.FindAlignmentImportHandler(file_name, format);

  ali_io->Import(a, file_name);
  return a;
}

seq::SequenceList SequenceListFromStream(std::istream&  stream, 
                                         const String& format)
{
  IOManager& m=IOManager::Instance();                               
  SequenceIOHandlerPtr ali_io=m.FindAlignmentImportHandler("", format);
  if (!ali_io) {
    throw IOException("No IO Handler for '"+format+"' format");
  }
  seq::SequenceList a=seq::CreateSequenceList();
  ali_io->Import(a, stream);
  return a;
}

seq::SequenceList SequenceListFromString(const String& data, 
                                         const String& format)
{
  std::stringstream stream(data, std::ios::in);
  return SequenceListFromStream(stream, format);
} 

seq::AlignmentHandle LoadAlignment(const String& file_name,
                                   const String& format) 
{
  return AlignmentFromSequenceList(LoadSequenceList(file_name, format));
}

seq::AlignmentHandle AlignmentFromStream(std::istream& stream, 
                                         const String& format)
{
  return AlignmentFromSequenceList(SequenceListFromStream(stream, format));
}

seq::AlignmentHandle AlignmentFromString(const String& data, 
                                         const String& format)
{
  return AlignmentFromSequenceList(SequenceListFromString(data, format));
}

seq::SequenceHandle get_one_seq(seq::SequenceList seq_list)
{
  if (seq_list.GetCount()!=1) {
    std::stringstream s;
    s << "Sequence file contains " << seq_list.GetCount() 
      << " sequences." <<  std::endl << "If you want to load multiple "
      << "sequences use io.LoadAlignment or io.LoadSequenceList";
    throw IOException(s.str());
  }
  return seq_list[0];  
}

seq::SequenceHandle LoadSequence(const String& file_name, 
                                 const String& format)
{
  return get_one_seq(LoadSequenceList(file_name, format));
}

seq::SequenceHandle SequenceFromStream(std::istream&  stream, 
                                       const String& format)
{
 return get_one_seq(SequenceListFromStream(stream, format));
}

seq::SequenceHandle SequenceFromString(const String& data, 
                                       const String& format)
{
 return get_one_seq(SequenceListFromString(data, format));
}

seq::ProfileHandlePtr LoadSequenceProfile(const String& file_name,
                                          const String& format)
{
  seq::ProfileHandlePtr prof(new seq::ProfileHandle);
  IOManager& m = IOManager::Instance();
  ProfileIOHandlerPtr prof_io = m.FindProfileImportHandler(file_name, format);
  prof_io->Import(*prof, file_name);
  return prof;
}

seq::ProfileHandlePtr SequenceProfileFromString(const String& data,
                                                const String& format)
{
  seq::ProfileHandlePtr prof(new seq::ProfileHandle);
  String file_name = "fake." + format;
  IOManager& m = IOManager::Instance();
  ProfileIOHandlerPtr prof_io = m.FindProfileImportHandler(file_name, format);
  prof_io->ImportFromString(*prof, data);
  return prof;
}

}}
