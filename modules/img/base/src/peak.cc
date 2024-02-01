//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2020 by the OpenStructure authors
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
  Author: Ansgar Philippsen
*/

#include <fstream>
#include <sstream>

#include <ost/message.hh>

#include "peak.hh"

namespace ost { namespace img { 

/*PointList PeakList2PointList(const PeakList& plist)
{
  PointList pl;
  for(PeakList::const_iterator it=plist.begin();it!=plist.end();++it) {
    pl.push_back(Point(*it));
  }
  return pl;
}*/

PeakList ImportPeakList(const String& filename)
{
  PeakList plist;

  std::ifstream infile(filename.c_str());
  if (!infile) throw Error(filename+" not found");

  String line;
  while (std::getline(infile,line)) {
    std::istringstream linestream(line);
    int px,py,pz;
    Real vv;
    linestream >> px >> py >> pz >> vv;
    if (linestream.fail()) {}
    plist.push_back(Peak(Point(px,py,pz),vv));
  }

  return plist;
}

void ExportPeakList(const PeakList& plist, const String& filename)
{
  std::ofstream outfile(filename.c_str());
  if (!outfile) throw Error(filename+" cannot be opened");

  for(PeakList::const_iterator it=plist.begin();it!=plist.end();++it) {
    const Peak& pk=*it;
    outfile << pk[0] << " " << pk[1] << " " << pk[2] << " " << pk.GetValue() << std::endl;
  }
}

}} // ns


