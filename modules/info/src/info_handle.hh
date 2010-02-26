//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2010 by the OpenStructure authors
// Copyright (C) 2003-2010 by the IPLT authors
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


/*
  high level info interface
  
  Author: Ansgar Philippsen
*/

#ifndef OST_DATA_INFO_H
#define OST_DATA_INFO_H

#include <ost/base.hh>

#include "info_path.hh"
#include "info_fw.hh"

namespace ost { namespace info {

class InfoHandle;

//! Create InfoHandle from scratch
DLLEXPORT InfoHandle CreateInfo();
DLLEXPORT InfoHandle CreateInfo(const String& dtdfile);

//! Load InfoHandle from a file
DLLEXPORT InfoHandle LoadInfo(const String& file);


//! main info handle
/*!
  following the handle concept, this class represents a thin wrapper
  to an underlying info class, which is shared among copies of InfoHandles
  unless the Copy() method is used.
*/
class DLLEXPORT InfoHandle {
  friend InfoHandle CreateInfo();
  friend InfoHandle CreateInfo(const String&);
  friend InfoHandle LoadInfo(const String&);

  typedef std::vector<RootPtr> RootPtrList;
public:
  //! empty, ie invalid handle
  InfoHandle();

  //! return clone of the underlying representation
  InfoHandle Copy() const;
  
  //! generate info from file
  void Import(const String& file);

  //! write content to file
  void Export(const String& file) const;

  bool IsValid() const;

  InfoGroup Root() const;

  void AddDefault(const InfoHandle& h);

  bool HasDefaultGroup(const InfoPath& p) const;
  InfoGroup GetDefaultGroup(const InfoPath& p) const;

  bool HasDefaultItem(const InfoPath& p) const;
  InfoItem GetDefaultItem(const InfoPath& p) const;

private:
  InfoHandle(RootPtr impl);

  RootPtr impl_;
};

}} // ns

#endif
