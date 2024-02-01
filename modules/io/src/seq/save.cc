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
#include <sstream>


#include <ost/io/io_manager.hh>

#include "save.hh"

/*
  Author: Marco Biasini
 */
namespace ost { namespace io {

void SaveAlignment(const seq::ConstSequenceList& msa, const String& filename,
                   const String& format)
{

}

String AlignmentToString(const seq::ConstSequenceList& ali, 
                         const String& format)
{
  IOManager& m=IOManager::Instance();                               
  SequenceIOHandlerPtr ali_io=m.FindAlignmentExportHandler("", format);
  
  if (!ali_io) {
    throw IOException("Can't export sequences to format '"+format+"'.");
  }
  std::stringstream sstream;
  ali_io->Export(ali, sstream);
  return sstream.str();
}

void DLLEXPORT_OST_IO SaveSequenceList(const seq::ConstSequenceList& seq_list,
                                       const String& filename,
                                       const String& format)
{  
  IOManager& m=IOManager::Instance();                               
  SequenceIOHandlerPtr ali_io=m.FindAlignmentExportHandler(filename, format);

  ali_io->Export(seq_list, filename);  
}

                                      
String DLLEXPORT_OST_IO 
SequenceListToString(const seq::ConstSequenceList& seq_list,
                     const String& format)
{
  IOManager& m=IOManager::Instance();                               
  SequenceIOHandlerPtr ali_io=m.FindAlignmentExportHandler("", format);
  
  if (!ali_io) {
    throw IOException("Can't export sequences to format '"+format+"'.");
  }
  std::stringstream sstream;
  ali_io->Export(seq_list, sstream);
  return sstream.str();
}
                                          
void  SaveSequence(const seq::ConstSequenceHandle& sequence, 
                   const String& filename, const String& format)
{
  seq::ConstSequenceList seq_list=seq::CreateConstSequenceList();
  seq_list.AddSequence(sequence);
  SaveSequenceList(seq_list, filename, format);
}

                                      
             
String SequenceToString(const seq::ConstSequenceHandle& sequence,
                        const String& format)
{
  seq::ConstSequenceList seq_list=seq::CreateConstSequenceList();
  seq_list.AddSequence(sequence);
  return SequenceListToString(seq_list, format);
}

void SaveAlignment(const seq::AlignmentHandle& msa,
                   const String& filename,
                   const String& format)
{
  SaveSequenceList(msa.GetSequences(), filename, format);
}

String AlignmentToString(const seq::AlignmentHandle& msa,
                         const String& format)
{
  return SequenceListToString(msa.GetSequences(), format);
}

}}
