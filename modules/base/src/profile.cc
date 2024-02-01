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
#include "profile.hh"

#if !defined(WIN32) && OST_PROFILING_ENABLED

#include <iostream>

#include "log.hh"

namespace ost {

Real delta_time(const timeval& t1, 
                  const timeval& t2)
{
  return 1e-6*static_cast<Real>((t2.tv_sec*1e6+t2.tv_usec)-
                                  (t1.tv_sec*1e6+t1.tv_usec));
}
Profile::Profile(const String& msg):
  msg_(msg)
{
  gettimeofday(&start_, NULL);
}

Profile::~Profile()
{
  struct timeval end;
  gettimeofday(&end,NULL);
  std::cout << msg_ << " finished in " <<delta_time(start_, end) 
            << "s" << std::endl;
}

}

#endif
    
